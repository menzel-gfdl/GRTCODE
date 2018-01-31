#include "GasProps.h"
#include "GaussianFuncs.h"
#include "IdaVoigtFuncs.h"
#include "line_shape.h"
#include "LorentzFuncs.h"
#include "myreal.h"
#include "omp.h"
#include "RfmVoigtFuncs.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#else
#include <math.h>
#endif

/*---------------------------------------------------------------------------*/
/*Compute the dimensionless optical depth values at each desired frequency.

  Arguments:
      molId        [in]      Molecule id.
      nL           [in]      Size of the line dimension for the inputted
                                 arrays.
      nF           [in]      Size of the frequency dimension for the outputted
                                 array.
      loWn         [in]      Lowest frequency where the optical depth is
                                 calculated (cm^-1).
      resolution   [in]      Frequency resolution for the optical depth
                                 values (cm^-1).
      numLayers    [in]      Size of the height dimension for the inputted
                                 arrays.
      breadth      [in]      Integer number of frequencies that each molecular
                                 line spans.
      T            [in]      Array of temperatures (K).  This array is stored
                                 as [height].
      Vnn          [in]      Array of spectral line transition frequencies
                                 (cm^-1).  This array is stored as [line].
      Gam          [in]      Array of pressure broadened line halfwidths
                                 (cm^-1).  This array is stored as
                                 [height][line].
      PShift       [in]      Array of pressure-shift corrections of the line
                                 positions (cm^-1).  This array is stored as
                                 [height][line].
      S            [in]      Array of corrected spectral line intensities (cm).
                                 This array is stored as [height][line].
      tauU_d       [in]      Array of integrated number densities (cm^-2).
                                 This array is stored as [height].
      out          [in,out]  Array of dimensionless optical depths.  This
                                 array is stored as [height][frequency].
*/

#ifdef __NVCC__
__global__ void eval_profile(unsigned int const molId,
                             unsigned int const nL,
                             int const nF,
                             REAL_t const loWn,
                             REAL_t const resolution,
                             unsigned int const numLayers,
                             unsigned int const breadth,
                             REAL_t const * const T,
                             REAL_t const * const Gam,
                             REAL_t const * const PShift,
                             REAL_t const * const S,
                             REAL_t const * const tauU_d,
                             REAL_t * const out)
{
    unsigned int ltid = blockIdx.x * blockDim.x + threadIdx.x;

    if (ltid < nL)
    {
        int ftid;
        int fcenterid;
        unsigned int lyr;
        unsigned int loffset;
        const int fsteps = ceil((REAL_t)breadth/resolution);
        REAL_t snn;
        REAL_t tauu;
        REAL_t temp;
        REAL_t molarMass;
        LineShapeInputs_t in;
        REAL_t line_shape;

#pragma unroll
        for (lyr=0;lyr<numLayers;++lyr)
        {
            loffset = lyr*nL + ltid;
            in.lineCenter = PShift[loffset];

            /*Find index of nearest frequency bin to line.*/
            fcenterid = (2*((in.lineCenter-loWn)/resolution)+1)/2;
            if (fcenterid >= 0 && fcenterid < nF)
            {
                snn = S[loffset];
                tauu = tauU_d[lyr];
                temp = T[lyr];
                molarMass = getMolarMass(molId);

                /*Set the necessary values for the line shape input
                  structure.*/
                in.lorHWHM = Gam[loffset];
                in.gauHWHM = gaussian_hwhm(temp,
                                           molarMass,
                                           in.lineCenter);
#ifdef IDA_VOIGT
                in.eta = eta(2.*in.lorHWHM,
                             2.*in.gauHWHM);
#else
                in.eta = -1;
#endif

                /*Calculate the optical depth values from the left edge of
                  the line to the line center.*/
#pragma unroll
                for (ftid=fcenterid-((int)fsteps);ftid<=fcenterid;++ftid)
                {
                    if (ftid >= 0)
                    {
                        in.freq = ((REAL_t)ftid)*resolution + loWn;

                        /*Calculate the value of the line shape function.*/
#if defined(GAUSSIAN)
                        line_shape = gaussian_line_shape(in);
#elif defined(LORENTZ)
                        line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        line_shape = ida_voigt_line_shape(in);
#else
                        line_shape = rfm_voigt_line_shape(in);
#endif

                        /*Atomics must be used for now because of a race on
                          load-alter-write out[ftid].*/
                        atomicAdd(&(out[lyr*nF+ftid]),
                                  snn*tauu*line_shape);
                    }
                }

                /*Calculate the optical depth values from the right edge of
                  the line to the line center.*/
#pragma unroll
                for (ftid=fcenterid+((int)fsteps);ftid>fcenterid;--ftid)
                {
                    if (ftid < nF)
                    {
                        in.freq = ((REAL_t)ftid)*resolution + loWn;

                        /*Calculate the value of the line shape function.*/
#if defined(GAUSSIAN)
                        line_shape = gaussian_line_shape(in);
#elif defined(LORENTZ)
                        line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        line_shape = ida_voigt_line_shape(in);
#else
                        line_shape = rfm_voigt_line_shape(in);
#endif

                        /*Atomics must be used for now because of a race on
                          load-alter-write out[ftid].*/
                        atomicAdd(&(out[lyr*nF+ftid]),
                                  snn*tauu*line_shape);
                    }
                }
            }
        }
    }

    return;
}

#endif

void eval_profile_h(unsigned int const molId,
                    unsigned int const nL,
                    int const nF,
                    REAL_t const loWn,
                    REAL_t const resolution,
                    unsigned int const numLayers,
                    unsigned int const breadth,
                    REAL_t const * const T,
                    REAL_t const * const Gam,
                    REAL_t const * const PShift,
                    REAL_t const * const S,
                    REAL_t const * const tauU_d,
                    REAL_t const * const pathlength_d,
                    REAL_t * const out)
{
    /*CPU specific code.*/
    unsigned int ltid;
    int ftid;
    int fcenterid;
    unsigned int lyr;
    unsigned int loffset;
    const int fsteps = ceil((REAL_t)breadth/resolution);
    REAL_t snn;
    REAL_t tauu;
    REAL_t len;
    REAL_t temp;
    REAL_t molarMass;
    LineShapeInputs_t in;
    REAL_t line_shape;

#pragma omp parallel for schedule(static) \
                         collapse(2) \
                         default(none) \
                         private(ltid,loffset,in,fcenterid,snn, \
                                 molarMass,ftid,line_shape, \
                                 lyr,tauu,len,temp)
/*
                         shared(nL,numLayers,PShift,loWn, \
                                resolution,nF,S, \
                                tauU_d,pathlength_d,T,molId, \
                                Gam,fsteps,out)
*/
    for (lyr=0;lyr<numLayers;++lyr)
    {
        for (ltid=0;ltid<nL;++ltid)
        {
            tauu = tauU_d[lyr];
            len = pathlength_d[lyr];
            temp = T[lyr];
            loffset = lyr*nL + ltid;
            in.lineCenter = PShift[loffset];

            /*Find index of nearest frequency bin to line.*/
            fcenterid = (2*((in.lineCenter-loWn)/resolution)+1)/2;
            if (fcenterid >= 0 && fcenterid < nF)
            {
                snn = S[loffset];
                molarMass = getMolarMass(molId);

                /*Set the necessary values for the line shape input
                  structure.*/
                in.lorHWHM = Gam[loffset];
                in.gauHWHM = gaussian_hwhm(temp,
                                           molarMass,
                                           in.lineCenter);
#ifdef IDA_VOIGT
                in.eta = eta(2.*in.lorHWHM,
                             2.*in.gauHWHM);
#else
                in.eta = -1;
#endif

                /*Calculate the optical depth values from the left edge of
                  the line to the line center.*/
                for (ftid=fcenterid-((int)fsteps);ftid<=fcenterid;++ftid)
                {
                    if (ftid >= 0)
                    {
                        in.freq = ((REAL_t)ftid)*resolution + loWn;

                        /*Calculate the value of the line shape function.*/
#if defined(GAUSSIAN)
                        line_shape = gaussian_line_shape(in);
#elif defined(LORENTZ)
                        line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        line_shape = ida_voigt_line_shape(in);
#else
                        line_shape = rfm_voigt_line_shape(in);
#endif

#pragma omp atomic update
                        out[lyr*nF+ftid] += snn*tauu*len*line_shape;
                    }
                }

                /*Calculate the optical depth values from the right edge of
                  the line to the line center.*/
                for (ftid=fcenterid+((int)fsteps);ftid>fcenterid;--ftid)
                {
                    if (ftid < nF)
                    {
                        in.freq = ((REAL_t)ftid)*resolution + loWn;

                        /*Calculate the value of the line shape function.*/
#if defined(GAUSSIAN)
                        line_shape = gaussian_line_shape(in);
#elif defined(LORENTZ)
                        line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        line_shape = ida_voigt_line_shape(in);
#else
                        line_shape = rfm_voigt_line_shape(in);
#endif

#pragma omp atomic update
                        out[lyr*nF+ftid] += snn*tauu*len*line_shape;
                    }
                }
            }
        }
    }

    return;
}
