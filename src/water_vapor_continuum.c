#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __NVCC__
#include "cuda_helpers.cuh"
#endif
#include "debug.h"
#include "floating_point_type.h"
#include "parse_csv.h"
#include "utils.h"
#include "water_vapor_continuum.h"


/*Read in the water vapor continuum coefficients.*/
int get_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc,
                                    char const * const h2o_ctm_dir,
                                    uint64_t const num_wpoints,
                                    double const w0,
                                    double const res)
{
    not_null(cc);
    not_null(h2o_ctm_dir);

    /*Set file names.*/
    char *filepath[NUM_COEFS];
    int num_vals[NUM_COEFS];
    size_t s = strlen(h2o_ctm_dir) + 64;
    int i;
    for (i=0;i<NUM_COEFS;++i)
    {
        check(malloc_ptr((void **)(&(filepath[i])),
                         sizeof(**filepath)*s));
        switch (i)
        {
            case MTCKD25_F296:
                snprintf(filepath[i],
                         s,
                         "%s/296MTCKD25_F.csv",
                         h2o_ctm_dir);
                num_vals[i] = 1;
                break;
            case MTCKD25_S296:
                snprintf(filepath[i],
                         s,
                         "%s/296MTCKD25_S.csv",
                         h2o_ctm_dir);
                num_vals[i] = 1;
                break;
            case CKDF:
                snprintf(filepath[i],
                         s,
                         "%s/CKDF.csv",
                         h2o_ctm_dir);
                num_vals[i] = 3;
                break;
            case CKDS:
                snprintf(filepath[i],
                         s,
                         "%s/CKDS.csv",
                         h2o_ctm_dir);
                num_vals[i] = 3;
                break;
            default:
                sentinel();
        }
    }

    /*Allocate memory for each of the coefficient pointer.*/
    cc->coefs = NULL;
    check(malloc_ptr((void **)(&(cc->coefs)),
                     sizeof(*(cc->coefs))*NUM_COEFS));
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
            fatal(VALUE_ERR,
                  "The number of columns (%d) in file %s does not match"
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
        int num_bytes = sizeof(*c)*num_wpoints;
        check(malloc_ptr((void **)&c,
                         num_bytes));
        memset(c,
               0,
               num_bytes);
        
        /*Interpolate to wavenumber grid.*/
        fp_t *x = &(fbuf[0]);
        fp_t *y = &(fbuf[num_lines]);
        for (j=0;(unsigned int)j<num_wpoints;++j)
        {
            check(linear_interpolation(x,
                                       y,
                                       num_lines,
                                       (fp_t)(w0 + j*res),
                                       &(c[j])));
        }
        free(fbuf);
        cc->coefs[i] = c;
    }
    cc->num_wpoints = num_wpoints;
    for (i=0;i<NUM_COEFS;++i)
    {
        free(filepath[i]);
    }
    return SUCCESS;
}


int free_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc)
{
    not_null(cc);
    int i;
    for (i=0;i<NUM_COEFS;++i)
    {
        free(cc->coefs[i]);
    }
    free(cc->coefs);
    cc->coefs = NULL;
    return SUCCESS;
}


int put_water_vapor_coefs_on_device(WaterVaporContinuumCoefs_t const * const in,
                                    WaterVaporContinuumCoefs_t * const out)
{
    not_null(in);
    not_null(out);
    check(malloc_ptr((void **)(&(out->coefs)),
                     sizeof(*(out->coefs))*NUM_COEFS));
    int i;
    for (i=0;i<NUM_COEFS;++i)
    {
#ifdef __NVCC__
        int num_bytes = sizeof(*(in->coefs[i]))*(in->num_wpoints);
        HANDLE_ERROR(cudaMalloc(&(out->coefs[i]),
                                num_bytes));
        HANDLE_ERROR(cudaMemcpy(out->coefs[i],
                                in->coefs[i],
                                num_bytes,
                                cudaMemcpyHostToDevice));
#endif
    }
    return SUCCESS;
}


int remove_water_vapor_coefs_from_device(WaterVaporContinuumCoefs_t * const in)
{
    not_null(in);
    int i;
    for (i=0;i<NUM_COEFS;++i)
    {
#ifdef __NVCC__
        HANDLE_ERROR(cudaFree(in->coefs[i]));
#endif
    }
    free(in->coefs);
    in->coefs = NULL;
    return SUCCESS;
}


#ifdef __NVCC__
__global__
void calc_water_vapor_ctm_optical_depth(uint64_t const num_wpoints,
                                        int const num_layers,
                                        fp_t * const tau,
                                        fp_t const * const CS,
                                        fp_t const * const T,
                                        fp_t const * const Ps,
                                        fp_t const * const N,
                                        fp_t const * const T0,
                                        fp_t const * const CF,
                                        fp_t const * const P,
                                        fp_t const * const T0F)
{
    unsigned int const tid = blockIdx.x*blockDim.x + threadIdx.x;

    if (tid < num_wpoints)
    {
        fp_t const tref = 296.0;
        int lyr;

#pragma unroll
        for (lyr=0;lyr<num_layers;++lyr)
        {
            tau[lyr*num_wpoints+tid] += N[lyr]*(tref/T[lyr])*((CS[tid]*Ps[lyr]*
                                        exp(T0[tid]*(tref-T[lyr]))) +
                                        (CF[tid]*(P[lyr]-Ps[lyr])*
                                        exp(T0F[tid]*(tref-T[lyr]))));
        }
    }
    return;
}
#endif


void calc_water_vapor_ctm_optical_depth_h(uint64_t const num_wpoints,
                                          int const num_layers,
                                          fp_t * const tau,
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
    for (lyr=0;lyr<num_layers;++lyr)
    {
        for (tid=0;tid<num_wpoints;++tid)
        {
            tau[lyr*num_wpoints+tid] += N[lyr]*(tref/T[lyr])*((CS[tid]*Ps[lyr]*
                                        exp(T0[tid]*(tref-T[lyr]))) +
                                        (CF[tid]*(P[lyr]-Ps[lyr])*
                                        exp(T0F[tid]*(tref-T[lyr]))));
        }
    }
    return;
}
