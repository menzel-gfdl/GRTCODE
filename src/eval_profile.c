#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "doppler.h"
#include "floating_point_type.h"
#include "gas_properties.h"
#include "ida_voigt.h"
#include "line_shape.h"
#include "lorentz.h"
#include "RFM_voigt.h"

#ifdef __NVCC__
#include "cuda_helpers.cuh"
#else
#include <math.h>
#endif


/*Compute the dimensionless optical depth values at each desired frequency.

  Arguments:
      mol_id        [in]      Molecule id.
      num_lines           [in]      Size of the line dimension for the inputted
                                 arrays.
      num_wpoints           [in]      Size of the frequency dimension for the outputted
                                 array.
      w0         [in]      Lowest frequency where the optical depth is
                                 calculated (cm^-1).
      wres   [in]      Frequency wres for the optical depth
                                 values (cm^-1).
      num_layers    [in]      Size of the height dimension for the inputted
                                 arrays.
      wcutoff      [in]      Integer number of frequencies that each molecular
                                 line spans.
      T            [in]      Array of temperatures (K).  This array is stored
                                 as [height].
      Vnn          [in]      Array of spectral line transition frequencies
                                 (cm^-1).  This array is stored as [line].
      gamma          [in]      Array of pressure broadened line halfwidths
                                 (cm^-1).  This array is stored as
                                 [height][line].
      Pshift       [in]      Array of pressure-shift corrections of the line
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
__global__ void eval_profile(int const mol_id,
                             unsigned int const num_lines,
                             uint64_t const num_wpoints,
                             double const w0,
                             double const wres,
                             int const num_layers,
                             double const wcutoff,
                             fp_t const * const T,
                             fp_t const * const gamma,
                             fp_t const * const Pshift,
                             fp_t const * const s,
                             fp_t const * const N,
                             fp_t * const tau)
{
    unsigned int const ltid = blockIdx.x * blockDim.x + threadIdx.x;
    if (ltid < num_lines)
    {
        int const fsteps = ceil(wcutoff/wres);
        int lyr;
#pragma unroll
        for (lyr=0;lyr<num_layers;++lyr)
        {
            unsigned int loffset = lyr*num_lines + ltid;
            LineShapeInputs_t in;
            in.line_center = Pshift[loffset];
            int fcenterid = (2*((in.line_center-w0)/wres)+1)/2;

            /*Find index of nearest frequency bin to line.*/
            if (fcenterid >= 0 && fcenterid < num_wpoints)
            {
                fp_t snn = s[loffset];
                fp_t n = N[lyr];
                fp_t temp = T[lyr];
                fp_t molar_mass = get_molar_mass(mol_id);

                /*Set the necessary values for the line shape input
                  structure.*/
                in.lorentz_hwhm = gamma[loffset];
                in.doppler_hwhm = doppler_hwhm(temp,
                                               molar_mass,
                                               in.line_center);
#ifdef IDA_VOIGT
                in.eta = eta(2.*in.lorentz_hwhm,
                             2.*in.doppler_hwhm);
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
                        in.w = ((fp_t)ftid)*wres + w0;

                        /*Calculate the value of the line shape function.*/
#if defined(DOPPLER)
                        fp_t line_shape = doppler_line_shape(in);
#elif defined(LORENTZ)
                        fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        fp_t line_shape = ida_voigt_line_shape(in);
#else
                        fp_t line_shape = rfm_voigt_line_shape(in);
#endif

                        /*Atomics must be used for now because of a race on
                          load-alter-write out[ftid].*/
                        atomicAdd(&(tau[lyr*num_wpoints+ftid]),
                                  snn*n*line_shape);
                    }
                }

#pragma unroll
                for (ftid=fcenterid+((int)fsteps);ftid>fcenterid;--ftid)
                {
                    /*Calculate the optical depth values from the right edge
                      of the line to the line center.*/
                    if (ftid < num_wpoints)
                    {
                        in.w = ((fp_t)ftid)*wres + w0;

                        /*Calculate the value of the line shape function.*/
#if defined(DOPPLER)
                        fp_t line_shape = doppler_line_shape(in);
#elif defined(LORENTZ)
                        fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        fp_t line_shape = ida_voigt_line_shape(in);
#else
                        fp_t line_shape = rfm_voigt_line_shape(in);
#endif

                        /*Atomics must be used for now because of a race on
                          load-alter-write out[ftid].*/
                        atomicAdd(&(tau[lyr*num_wpoints+ftid]),
                                  snn*n*line_shape);
                    }
                }
            }
        }
    }
    return;
}
#endif


void eval_profile_h(int const mol_id,
                    unsigned int const num_lines,
                    uint64_t const num_wpoints,
                    double const w0,
                    double const wres,
                    int const num_layers,
                    double const wcutoff,
                    fp_t const * const T,
                    fp_t const * const gamma,
                    fp_t const * const Pshift,
                    fp_t const * const s,
                    fp_t const * const N,
                    fp_t * const tau)
{
    int const fsteps = ceil(wcutoff/wres);
    int lyr;
    unsigned int ltid;

#pragma omp parallel for schedule(static) \
                         collapse(2) \
                         default(none) \
                         private(lyr,ltid)
    for (lyr=0;lyr<num_layers;++lyr)
    {
        for (ltid=0;ltid<num_lines;++ltid)
        {
            fp_t n = N[lyr];
            fp_t temp = T[lyr];
            unsigned int loffset = lyr*num_lines + ltid;
            LineShapeInputs_t in;
            in.line_center = Pshift[loffset];

            /*Find index of nearest frequency bin to line.*/
            unsigned int fcenterid = (2*((in.line_center-w0)/wres)+1)/2;
            if (fcenterid < num_wpoints)
            {
                fp_t snn = s[loffset];
                fp_t molar_mass = get_molar_mass(mol_id);

                /*Set the necessary values for the line shape input
                  structure.*/
                in.lorentz_hwhm = gamma[loffset];
                in.doppler_hwhm = doppler_hwhm(temp,
                                               molar_mass,
                                               in.line_center);
#ifdef IDA_VOIGT
                in.eta = eta(2.*in.lorentz_hwhm,
                             2.*in.doppler_hwhm);
#else
                in.eta = -1;
#endif

                /*Calculate the optical depth values from the left edge of
                  the line to the line center.*/
                int ftid;
                for (ftid=fcenterid-((int)fsteps);ftid<=(int)fcenterid;++ftid)
                {
                    if (ftid >= 0)
                    {
                        in.w = ((fp_t)ftid)*wres + w0;

                        /*Calculate the value of the line shape function.*/
#if defined(DOPPLER)
                        fp_t line_shape = doppler_line_shape(in);
#elif defined(LORENTZ)
                        fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        fp_t line_shape = ida_voigt_line_shape(in);
#else
                        fp_t line_shape = rfm_voigt_line_shape(in);
#endif

#pragma omp atomic update
                        tau[lyr*num_wpoints+ftid] += snn*n*line_shape;
                    }
                }

                /*Calculate the optical depth values from the right edge of
                  the line to the line center.*/
                for (ftid=fcenterid+((int)fsteps);ftid>(int)fcenterid;--ftid)
                {
                    if (ftid < (int)num_wpoints)
                    {
                        in.w = ((fp_t)ftid)*wres + w0;

                        /*Calculate the value of the line shape function.*/
#if defined(DOPPLER)
                        fp_t line_shape = doppler_line_shape(in);
#elif defined(LORENTZ)
                        fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                        fp_t line_shape = ida_voigt_line_shape(in);
#else
                        fp_t line_shape = rfm_voigt_line_shape(in);
#endif

#pragma omp atomic update
                        tau[lyr*num_wpoints+ftid] += snn*n*line_shape;
                    }
                }
            }
        }
    }

    return;
}
