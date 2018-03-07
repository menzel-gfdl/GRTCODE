#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "continuum.h"
#include "debug.h"
#include "floating_point_type.h"
#include "omp.h"
#include "utils.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif

static const int MAXCHARSPERLINE=128;


static int parse_CKD(char const * const fname,
                     fp_t *AryPtr,
                     int const maxwavenum,
                     int const minw,
                     fp_t const res)
{
    not_null(fname);
    not_null(AryPtr);

    /*Open the file.*/
    log_mesg("Reading in continuum coefficients from file %s.",
             fname);
    FILE *f = NULL;
    open_file(f,
              fname,
              "r");

    /*Count the number of lines in the file.*/
    unsigned int line_count = 0;
    char line[MAXCHARSPERLINE];
    while (fgets(line,MAXCHARSPERLINE,f) != NULL)
    {
        line_count++;
    }

    /*Read in the data.*/
    int *wavenums = NULL;
    check(malloc_ptr((void **)(&wavenums),
                     sizeof(*wavenums)*line_count));
    fp_t *buf = NULL;
    check(malloc_ptr((void **)(&buf),
                     sizeof(*buf)*line_count));
    int count = 0;
    double v0;
    double v1;
    rewind(f);
    while (fgets(line,MAXCHARSPERLINE,f) != NULL)
    {
        sscanf(line,
               "%lf %lf",
               &v0,
               &v1);
        wavenums[count] = (int)v0;
        buf[count] = (fp_t)v1;
        count++;
    }

    /*Close the file.*/
    count = fclose(f);
    if (count != 0)
    {
        fatal("Attempt to close file %s failed with return code %d.",
              fname,
              count);
    }

    /*Interpolate and store the values.*/
    for (count=0;count<maxwavenum;count++)
    {
        fp_t w = minw + count*res;
        if (w < wavenums[0])
        {
            AryPtr[count] = buf[0];
        }
        else if (w > wavenums[line_count-1])
        {
            AryPtr[count] = buf[line_count-1];
        }
        else
        {
            /*Binary search.*/
            int left = 0;
            int right = line_count-1;
            int match = 0;
            int mid;
            while (1)
            {
                mid = (right+left)/2;
                if (w == wavenums[mid])
                {
                    match = 1;
                    break;
                }
                else if (w < wavenums[mid])
                {
                    right = mid;
                }
                else
                {
                    left = mid;
                }
                if (right - left == 1)
                {
                    break;
                }
                else if (right-left == 0)
                {
                    fatal("wave number %e not contained in input file %s.",
                          w,
                          fname);
                }
            }
            if (match)
            {
                AryPtr[count] = buf[mid];
            }
            else
            {
                /*Linear interpolation.*/
                fp_t m = (buf[right]-buf[left])/
                         (wavenums[right]-wavenums[left]);
                fp_t b = buf[right] - m*wavenums[right];
                AryPtr[count] = w*m + b;
            }
        }
    }
    free(buf);
    free(wavenums);
    return SUCCESS;
}


int get_h2o_continuum_coefs(ContinuumCoefs_t *h2o,
                            unsigned int const nws,
                            int const w,
                            double const res,
                            int put_on_device)
{
    not_null(h2o);

    /*Set file names.*/
    char **h2o_coef_files = NULL;
    check(malloc_ptr((void **)(&h2o_coef_files),
                     sizeof(*h2o_coef_files)*NUM_COEF));
    int fname_len = 64;
    int i;
    for (i=0;i<NUM_COEF;++i)
    {
        check(malloc_ptr((void **)(&(h2o_coef_files[i])),
                         sizeof(*(h2o_coef_files[i]))*fname_len));
    }
    snprintf(h2o_coef_files[CS],
             fname_len,
             "INPUT/continuum/296MTCKD25_S.ccf");
    snprintf(h2o_coef_files[CF],
             fname_len,
             "INPUT/continuum/296MTCKD25_F.ccf");
    snprintf(h2o_coef_files[T0],
             fname_len,
             "INPUT/continuum/CKDS.ppp");
    snprintf(h2o_coef_files[T0F],
             fname_len,
             "INPUT/continuum/CKDF.ppp");

    /*Read in the coefficients.*/
    check(malloc_ptr((void **)(&(h2o->coefs)),
                     sizeof(*(h2o->coefs))*NUM_COEF));
    for (i=0;i<NUM_COEF;++i)
    {
        fp_t *buf = NULL;
        check(malloc_ptr((void **)(&buf),
                         sizeof(*buf)* nws));
        check(parse_CKD(h2o_coef_files[i],
                        buf,
                        nws,
                        w,
                        res));
        if (put_on_device)
        {
            using_gpu();
#ifdef __NVCC__
            HANDLE_ERROR(cudaMalloc(&(h2o->coefs[i]),
                                    sizeof(*(h2o->coefs[i]))*nws));
            HANDLE_ERROR(cudaMemcpy(h2o->coefs[i],
                                    buf,
                                    nws*sizeof(*buf),
                                    cudaMemcpyHostToDevice));
#endif
            free(buf);
        }
        else
        {
            h2o->coefs[i] = buf;
        }
    }

    /*Clean up.*/
    for (i=0;i<NUM_COEF;++i)
    {
        free(h2o_coef_files[i]);
    }
    free(h2o_coef_files);
    return SUCCESS;
}


int free_continuum_coeffs(ContinuumCoefs_t *c,
                          int const on_device)
{
    not_null(c);
    int i;
    for (i=0;i<NUM_COEF;++i)
    {
        if (on_device)
        {
            using_gpu();
#ifdef __NVCC__
            HANDLE_ERROR(cudaFree(c->coefs[i]));
#endif
        }
        else
        {
            free(c->coefs[i]);
        }
    }
    free(c->coefs);
    c->coefs = NULL;
    return SUCCESS;
}


#ifdef __NVCC__
__global__
void calc_ctm_optdepth(unsigned int const nF,
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


void calc_ctm_optdepth_h(unsigned int const nF,
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
                         private(lyr,tid) \
                         shared(numLayers,nF,optdepth,CS,tref, \
                                T,Ps, \
                                Z,T0,kB,AtmToPa,CmToM,CF,P,T0F)
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
