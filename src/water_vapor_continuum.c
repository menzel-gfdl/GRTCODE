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
                                    double const res,
                                    int const gpu_id)
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
        gmalloc(filepath[i],s,HOST_ONLY);
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
    gmalloc(cc->coefs,NUM_COEFS,HOST_ONLY);
    for (i=0;i<NUM_COEFS;++i)
    {
        /*Read in the data.*/
        log_info("Reading in water vapor continuum coefficients from"
                     " file %s.",
                 filepath[i]);
        int num_lines;
        int num_cols;
        char **buf;
        throw(parse_csv(filepath[i],
                        &num_lines,
                        &num_cols,
                        1,
                        &buf));
        if ((num_vals[i] + 1) != num_cols)
        {
            raise(VALUE_ERR,
                  "The number of columns (%d) in file %s does not match"
                      " the expected number (%d).",
                  num_cols,
                  filepath[i],
                  num_vals[i] + 1);
        }

        /*Convert the data from strings to floating point.*/
        fp_t *fbuf;
        int data_size = num_lines*num_cols;
        gmalloc(fbuf,data_size,HOST_ONLY);
        int j;
        for (j=0;j<data_size;++j)
        {
            double d;
            throw(to_double(buf[j],
                            &d));
            throw(to_fp_t(d,
                          &(fbuf[j])));
            gfree(buf[j],HOST_ONLY);
        }
        gfree(buf,HOST_ONLY);

        /*Allocate space for the coefficient values at each wavenumber.*/
        fp_t *c;
        gmalloc(c,num_wpoints,HOST_ONLY);
        gmemset(c,0,num_wpoints,HOST_ONLY);

        /*Interpolate to wavenumber grid.*/
        fp_t *x = &(fbuf[0]);
        fp_t *y = &(fbuf[num_lines]);
        for (j=0;(unsigned int)j<num_wpoints;++j)
        {
            throw(linear_interpolation(x,
                                       y,
                                       num_lines,
                                       (fp_t)(w0 + j*res),
                                       &(c[j])));
        }
        gfree(fbuf,HOST_ONLY);
        if (gpu_id == HOST_ONLY)
        {
            cc->coefs[i] = c;
        }
        else
        {
            gmalloc(cc->coefs[i],num_wpoints,gpu_id);
            gmemcpy(cc->coefs[i],c,num_wpoints,gpu_id,FROM_HOST);
            gfree(c,HOST_ONLY);
        }
    }
    cc->num_wpoints = num_wpoints;
    cc->gpu_id = gpu_id;
    for (i=0;i<NUM_COEFS;++i)
    {
        gfree(filepath[i],HOST_ONLY);
    }
    return SUCCESS;
}


int free_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc)
{
    not_null(cc);
    int i;
    for (i=0;i<NUM_COEFS;++i)
    {
        gfree(cc->coefs[i],cc->gpu_id);
    }
    gfree(cc->coefs,HOST_ONLY);
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
