#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "debug.h"
#include "floating_point_type.h"
#include "parse_csv.h"
#include "utils.h"
#include "water_vapor_continuum.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif


/*Read in the water vapor continuum coefficients.*/
int get_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc,
                                    unsigned int const nws,
                                    int const w0,
                                    double const res,
                                    int put_on_device)
{
    not_null(cc);

    /*Set file names.*/
    char *filepath[NUM_COEFS];
    filepath[MTCKD25_F296] = "INPUT/water_vapor_continuum/296MTCKD25_F.csv";
    filepath[MTCKD25_S296] = "INPUT/water_vapor_continuum/296MTCKD25_S.csv";
    filepath[CKDF] = "INPUT/water_vapor_continuum/CKDF.csv";
    filepath[CKDS] = "INPUT/water_vapor_continuum/CKDS.csv";
    int num_vals[NUM_COEFS];
    num_vals[MTCKD25_F296] = 1;
    num_vals[MTCKD25_S296] = 1;
    num_vals[CKDF] = 3;
    num_vals[CKDS] = 3;

    /*Allocate memory for each of the coefficient pointer.*/
    cc->coefs = NULL;
    check(malloc_ptr((void **)(&(cc->coefs)),
                     sizeof(*(cc->coefs))*NUM_COEFS));

    int i;
    for (i=0;i<NUM_COEFS;++i)
    {
        /*Read in the data.*/
        log_mesg("Reading in water vapor continuum coefficients from"
                     " file %s.",
                 filepath[i]);
        int num_lines;
        int num_cols;
        char **buf;
        check(parse_csv(filepath[i],
                        &num_lines,
                        &num_cols,
                        1,
                        &buf));
        if ((num_vals[i] + 1) != num_cols)
        {
            fatal("The number of columns (%d) in file %s does not match"
                      " the expected number (%d).",
                  num_cols,
                  filepath[i],
                  num_vals[i] + 1);
        }

        /*Convert the data from strings to floating point.*/
        fp_t *fbuf = NULL;
        int data_size = num_lines*num_cols;
        check(malloc_ptr((void **)(&fbuf),
                         sizeof(*fbuf)*data_size));
        int j;
        for (j=0;j<data_size;++j)
        {
            double d;
            check(to_double(buf[j],
                            &d));
            check(to_fp_t(d,
                          &(fbuf[j])));
            free(buf[j]);
        }
        free(buf);

        /*Allocate space for the coefficient values at each wavenumber.*/
        fp_t *c = NULL;
        int num_bytes = sizeof(*c)*nws;
        check(malloc_ptr((void **)&c,
                         num_bytes));
        memset(c,
               0,
               num_bytes);
        
        /*Interpolate to wavenumber grid.*/
        fp_t *x = &(fbuf[0]);
        fp_t *y = &(fbuf[num_lines]);
        for (j=0;j<nws;++j)
        {
            check(linear_interpolation(x,
                                       y,
                                       num_lines,
                                       (fp_t)(w0 + j*res),
                                       &(c[j])));
        }
        free(fbuf);

        if (put_on_device)
        {
            using_gpu();
#ifdef __NVCC__
            HANDLE_ERROR(cudaMalloc(&(cc->coefs[i]),
                                    num_bytes));
            HANDLE_ERROR(cudaMemcpy(cc->coefs[i],
                                    c,
                                    num_bytes,
                                    cudaMemcpyHostToDevice));
#endif
            free(c);
        }
        else
        {
            cc->coefs[i] = c;
        }
    }
    return SUCCESS;
}


int free_water_vapor_continuum_coeffs(WaterVaporContinuumCoefs_t *cc,
                                      int const on_device)
{
    not_null(cc);
    int i;
    for (i=0;i<NUM_COEFS;++i)
    {
        if (on_device)
        {
            using_gpu();
#ifdef __NVCC__
            HANDLE_ERROR(cudaFree(cc->coefs[i]));
#endif
        }
        else
        {
            free(cc->coefs[i]);
        }
    }
    free(cc->coefs);
    cc->coefs = NULL;
    return SUCCESS;
}


#ifdef __NVCC__
__global__
void calc_water_vapor_ctm_optdepth(unsigned int const nF,
                                   int const numLayers,
                                   fp_t * const optdepth,
                                   fp_t const * const CS,
                                   fp_t const * const T,
                                   fp_t const * const Ps,
                                   fp_t const * const N,
                                   fp_t const * const T0,
                                   fp_t const * const CF,
                                   fp_t const * const P,
                                   fp_t const * const T0F)
{
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;

    if (tid < nF)
    {
        fp_t const tref = 296.0;
        int lyr;

#pragma unroll
        for (lyr=0;lyr<numLayers;++lyr)
        {
            optdepth[lyr*nF+tid] += N[lyr]*(tref/T[lyr])*((CS[tid]*Ps[lyr]*
                                    exp(T0[tid]*(tref-T[lyr]))) +
                                    (CF[tid]*(P[lyr]-Ps[lyr])*
                                    exp(T0F[tid]*(tref-T[lyr]))));
        }
    }
    return;
}
#endif


void calc_water_vapor_ctm_optdepth_h(unsigned int const nF,
                                     int const numLayers,
                                     fp_t * const optdepth,
                                     fp_t const * const CS,
                                     fp_t const * const T,
                                     fp_t const * const Ps,
                                     fp_t const * const N,
                                     fp_t const * const T0,
                                     fp_t const * const CF,
                                     fp_t const * const P,
                                     fp_t const * const T0F)
{
    unsigned int tid;
    int lyr;
    fp_t const tref = 296.0;

#pragma omp parallel for collapse(2) \
                         schedule(static) \
                         default(none) \
                         private(lyr,tid)
    for (lyr=0;lyr<numLayers;++lyr)
    {
        for (tid=0;tid<nF;++tid)
        {
            optdepth[lyr*nF+tid] += N[lyr]*(tref/T[lyr])*((CS[tid]*Ps[lyr]*
                                    exp(T0[tid]*(tref-T[lyr]))) +
                                    (CF[tid]*(P[lyr]-Ps[lyr])*
                                    exp(T0F[tid]*(tref-T[lyr]))));
        }
    }
    return;
}
