#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "continuum.h"
#include "myreal.h"
#include "omp.h"

static const int MAXCHARSPERLINE=128;

/*---------------------------------------------------------------------------*/
/*Parse the coefficients out of the ascii file and interpolate.*/
void parseCKD(const char fname[],
              REAL_t *AryPtr,
              const int maxwavenum,
              const int minw,
              REAL_t const res)
{
    if (AryPtr == NULL)
    {
        fprintf(stderr,
                "Please first malloc AryPtr.  Aborting\n");
            exit(EXIT_FAILURE);
    }

    /*Open the file.*/
    FILE *F = fopen(fname,
                    "r");
    if (F == NULL)
    {
        fprintf(stderr,
                "Open %s Failed.  Aborting\n",
                fname);
        exit(EXIT_FAILURE);
    }

    /*Count the number of lines in the file.*/
    unsigned int line_count = 0;
    char line[MAXCHARSPERLINE];
    while (fgets(line,MAXCHARSPERLINE,F) != NULL)
    {
        line_count++;
    }

    /*Read in the data.*/
    int *wavenums = (int *)malloc(line_count*sizeof(int));
    REAL_t *buf = (REAL_t *)malloc(line_count*sizeof(REAL_t));
    if (wavenums == NULL || buf == NULL)
    {
        fprintf(stderr,
                "malloc of %zu bytes for wavenums or buf failed.",
                line_count*sizeof(REAL_t));
    }

    int count = 0;
    double v0;
    double v1;
    rewind(F);
    while (fgets(line,MAXCHARSPERLINE,F) != NULL)
    {
        sscanf(line,
               "%lf %lf",
               &v0,
               &v1);
        wavenums[count] = (int)v0;
        buf[count] = (REAL_t)v1;
        count++;
    }

    /*Close the file.*/
    count = fclose(F);
    if (count != 0)
    {
        fprintf(stderr,
                "Attempt to close file %s previously opened with handle"
                    " at %p failed with %d\n\tAborting.\n",
                fname,
                F,
                count);
        exit(EXIT_FAILURE);
    }

    /*Interpolate and store the values.*/
    REAL_t w;
    int left;
    int right;
    int mid;
    int match;
    REAL_t m;
    REAL_t b;
    for (count=0;count<maxwavenum;count++)
    {
        w = minw + count*res;
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
            left = 0;
            right = line_count-1;
            match = 0;
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
                    fprintf(stderr,
                            "wave number %e not contained in input file.",
                            w);
                    exit(EXIT_FAILURE);
                }
            }

            if (match)
            {
                AryPtr[count] = buf[mid];
            }
            else
            {
                m = (buf[right]-buf[left])/(wavenums[right]-wavenums[left]);
                b = buf[right] - m*wavenums[right];
                AryPtr[count] = w*m + b;
            }
        }
    }

    /*Clean up.*/
    free(buf);
    free(wavenums);
}

/*---------------------------------------------------------------------------*/
#ifdef __NVCC__
__global__
void calc_ctm_optdepth(unsigned int const nF,
                       unsigned int const numLayers,
                       REAL_t * const optdepth,
                       REAL_t const * const CS,
                       REAL_t const * const T,
                       REAL_t const * const PS_H2O,
                       REAL_t const * const Z,
                       REAL_t const * const T0,
                       REAL_t const * const CF,
                       REAL_t const * const P,
                       REAL_t const * const T0F)
{
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    unsigned int lyr;
    REAL_t const tref = 296.0;
    REAL_t const kB = 1.3806E-19;
    REAL_t const pconst = 1013.25;

    if (tid < nF)
    {
#pragma unroll
        for (lyr=0;lyr<numLayers;++lyr)
        {
            optdepth[lyr*nF+tid] += (CS[tid]*(tref/T[lyr])*PS_H2O[lyr]*
                                        PS_H2O[lyr]*Z[lyr]*
                                        exp(T0[tid]*(tref-T[lyr])))/
                                        (T[lyr]*kB*pconst) +
                                        (CF[tid]*(tref/T[lyr])*PS_H2O[lyr]*
                                        (P[lyr]-PS_H2O[lyr])*Z[lyr]*
                                        exp(T0F[tid]*(tref-T[lyr])))/
                                        (T[lyr]*kB*pconst);
        }
    }

    return;
}
#endif

/*---------------------------------------------------------------------------*/
void calc_ctm_optdepth_h(unsigned int const nF,
                         unsigned int const numLayers,
                         REAL_t * const optdepth,
                         REAL_t const * const CS,
                         REAL_t const * const T,
                         REAL_t const * const PS_H2O,
                         REAL_t const * const Z,
                         REAL_t const * const T0,
                         REAL_t const * const CF,
                         REAL_t const * const P,
                         REAL_t const * const T0F)
{
    unsigned int tid;
    unsigned int lyr;
    REAL_t const tref = 296.0;
    REAL_t const kB = 1.3806E-19;
    REAL_t const pconst = 1013.25;

#pragma omp parallel for collapse(2) default(none) \
                                     private(lyr,tid) \
                                     shared(numLayers,nF,optdepth,CS,tref, \
                                            T,PS_H2O, \
                                            Z,T0,kB,pconst,CF,P,T0F)
    for (lyr=0;lyr<numLayers;++lyr)
    {
        for (tid=0;tid<nF;++tid)
        {
            optdepth[lyr*nF+tid] += (CS[tid]*(tref/T[lyr])*PS_H2O[lyr]*
                                        PS_H2O[lyr]*Z[lyr]*
                                        exp(T0[tid]*(tref-T[lyr])))/
                                        (T[lyr]*kB*pconst) +
                                        (CF[tid]*(tref/T[lyr])*PS_H2O[lyr]*
                                        (P[lyr]-PS_H2O[lyr])*Z[lyr]*
                                        exp(T0F[tid]*(tref-T[lyr])))/
                                        (T[lyr]*kB*pconst);
        }
    }

    return;
}
