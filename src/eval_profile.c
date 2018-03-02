#include "floating_point_type.h"
#include "GasProps.h"
#include "GaussianFuncs.h"
#include "IdaVoigtFuncs.h"
#include "line_shape.h"
#include "LorentzFuncs.h"
#include "omp.h"
#include "RfmVoigtFuncs.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#else
#include <math.h>
#endif


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
      N            [in]      Array of integrated number densities (cm^-2).
                                 This array is stored as [height].
      tau          [in,out]  Array of dimensionless optical depths.  This
                                 array is stored as [height][frequency].
*/
#ifdef __NVCC__
__global__ void eval_profile(int const molId,
                             unsigned int const nL,
                             unsigned int const nF,
                             fp_t const loWn,
                             fp_t const resolution,
                             int const numLayers,
                             int const breadth,
                             fp_t const * const T,
                             fp_t const * const Gam,
                             fp_t const * const PShift,
                             fp_t const * const S,
                             fp_t const * const N,
                             fp_t * const tau)
{
    unsigned int ltid = blockIdx.x * blockDim.x + threadIdx.x;
    if (ltid < nL)
    {
        int const fsteps = ceil((fp_t)breadth/resolution);

#pragma unroll
        int lyr;
        for (lyr=0;lyr<numLayers;++lyr)
        {
            unsigned int loffset = lyr*nL + ltid;
            LineShapeInputs_t in;
            in.lineCenter = PShift[loffset];
            int fcenterid = (2*((in.lineCenter-loWn)/resolution)+1)/2;

            /*Find index of nearest frequency bin to line.*/
            if (fcenterid >= 0 && fcenterid < nF)
            {
                fp_t snn = S[loffset];
                fp_t n = N[lyr];
                fp_t temp = T[lyr];
                fp_t molarMass = getMolarMass(molId);

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

                int ftid;
#pragma unroll
                for (ftid=fcenterid-((int)fsteps);ftid<=fcenterid;++ftid)
                {
                    /*Calculate the optical depth values from the left edge
                      of the line to the line center.*/
                    if (ftid >= 0)
                    {
                        in.freq = ((fp_t)ftid)*resolution + loWn;

                        /*Calculate the value of the line shape function.*/
#if defined(GAUSSIAN)
                        fp_t line_shape = gaussian_line_shape(in);
#elif defined(LORENTZ)
                        fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        fp_t line_shape = ida_voigt_line_shape(in);
#else
                        fp_t line_shape = rfm_voigt_line_shape(in);
#endif

                        /*Atomics must be used for now because of a race on
                          load-alter-write out[ftid].*/
                        atomicAdd(&(tau[lyr*nF+ftid]),
                                  snn*n*line_shape);
                    }
                }

#pragma unroll
                for (ftid=fcenterid+((int)fsteps);ftid>fcenterid;--ftid)
                {
                    /*Calculate the optical depth values from the right edge
                      of the line to the line center.*/
                    if (ftid < nF)
                    {
                        in.freq = ((fp_t)ftid)*resolution + loWn;

                        /*Calculate the value of the line shape function.*/
#if defined(GAUSSIAN)
                        fp_t line_shape = gaussian_line_shape(in);
#elif defined(LORENTZ)
                        fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        fp_t line_shape = ida_voigt_line_shape(in);
#else
                        fp_t line_shape = rfm_voigt_line_shape(in);
#endif

                        /*Atomics must be used for now because of a race on
                          load-alter-write out[ftid].*/
                        atomicAdd(&(tau[lyr*nF+ftid]),
                                  snn*n*line_shape);
                    }
                }
            }
        }
    }
    return;
}
#endif


void eval_profile_h(int const molId,
                    unsigned int const nL,
                    unsigned int const nF,
                    fp_t const loWn,
                    fp_t const resolution,
                    int const numLayers,
                    int const breadth,
                    fp_t const * const T,
                    fp_t const * const Gam,
                    fp_t const * const PShift,
                    fp_t const * const S,
                    fp_t const * const N,
                    fp_t * const tau)
{
    int const fsteps = ceil((fp_t)breadth/resolution);
    int lyr;
    unsigned int ltid;

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
            fp_t n = N[lyr];
            fp_t temp = T[lyr];
            unsigned int loffset = lyr*nL + ltid;
            LineShapeInputs_t in;
            in.lineCenter = PShift[loffset];

            /*Find index of nearest frequency bin to line.*/
            int fcenterid = (2*((in.lineCenter-loWn)/resolution)+1)/2;
            if (fcenterid >= 0 && fcenterid < nF)
            {
                fp_t snn = S[loffset];
                fp_t molarMass = getMolarMass(molId);

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
                int ftid;
                for (ftid=fcenterid-((int)fsteps);ftid<=fcenterid;++ftid)
                {
                    if (ftid >= 0)
                    {
                        in.freq = ((fp_t)ftid)*resolution + loWn;

                        /*Calculate the value of the line shape function.*/
#if defined(GAUSSIAN)
                        fp_t line_shape = gaussian_line_shape(in);
#elif defined(LORENTZ)
                        fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        fp_t line_shape = ida_voigt_line_shape(in);
#else
                        fp_t line_shape = rfm_voigt_line_shape(in);
#endif

#pragma omp atomic update
                        tau[lyr*nF+ftid] += snn*n*line_shape;
                    }
                }

                /*Calculate the optical depth values from the right edge of
                  the line to the line center.*/
                for (ftid=fcenterid+((int)fsteps);ftid>fcenterid;--ftid)
                {
                    if (ftid < nF)
                    {
                        in.freq = ((fp_t)ftid)*resolution + loWn;

                        /*Calculate the value of the line shape function.*/
#if defined(GAUSSIAN)
                        fp_t line_shape = gaussian_line_shape(in);
#elif defined(LORENTZ)
                        fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        fp_t line_shape = ida_voigt_line_shape(in);
#else
                        fp_t line_shape = rfm_voigt_line_shape(in);
#endif

#pragma omp atomic update
                        tau[lyr*nF+ftid] += snn*n*line_shape;
                    }
                }
            }
        }
    }

    return;
}
