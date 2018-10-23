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


#ifdef __NVCC__
__global__
void eval_profile(int const mol_id, /**< Molecule id.*/
                  unsigned int const num_lines, /**< Number of molecular
                                                     lines.*/
                  uint64_t const num_wpoints_fine, /**< Size of the fine
                                                        mesh.*/
                  uint64_t const num_wpoints_coarse, /**< Size of the coarse
                                                          mesh.*/
                  double const w0, /**< Lowest allowed wavenumber [1/cm].*/
                  double const wres_fine, /**< Resolution of the fine
                                               mesh [1/cm].*/
                  double const wres_coarse, /**< Resoluion of the coarse
                                                 mesh [1/cm].*/
                  int const num_layers, /**< Number of atmospheric layers.*/
                  double const wcutoff, /**< Cutoff from line center [1/cm].*/
                  fp_t const * const T, /**< Layer temperatures [K].*/
                  fp_t const * const gamma, /**< Pressure broadened line
                                                 halfwidths [1/cm].*/
                  fp_t const * const Pshift, /**< Pressure shifted line
                                                  centers [1/cm].*/
                  fp_t const * const s, /**< Line strengths [cm^2].*/
                  fp_t const * const N, /**< Integrated Layer number
                                             densities [cm^-2].*/
                  fp_t * const tau_fine, /**< Optical depths on the fine
                                              mesh.*/
                  fp_t * const tau_coarse, /**< Optical depths on the
                                                coarse mesh.*/
                  fp_t const fine_factor /**< Cutoff from the line center
                                              for the fine mesh in units
                                              of Pressure broadened line
                                              halwidths.*/
                 )
{
    fp_t const molar_mass = get_molar_mass(mol_id);
    unsigned int const ltid = blockIdx.x*blockDim.x + threadIdx.x;
    if (ltid < num_lines)
    {
        int lyr;
        for (lyr=0;lyr<num_layers;++lyr)
        {
            unsigned int loffset = lyr*num_lines + ltid;
            LineShapeInputs_t in;
            in.line_center = Pshift[loffset];
            fp_t n = N[lyr];
            fp_t temp = T[lyr];
            fp_t snn = s[loffset];
            in.lorentz_hwhm = gamma[loffset];
            in.doppler_hwhm = doppler_hwhm(temp,
                                           molar_mass,
                                           in.line_center);
#ifdef IDA_VOIGT
            in.eta = eta(2.*in.lorentz_hwhm,
                         2.*in.doppler_hwhm);
#endif
/*
            double fcutoff = fine_factor*in.lorentz_hwhm;
*/
            double fcutoff = fine_factor;
            if (fcutoff > wcutoff)
            {
                fcutoff = wcutoff;
            }
            double wleft = in.line_center - fcutoff;
            if (wleft < w0)
            {
                wleft = w0;
            }
            double wright = in.line_center + fcutoff;
            if (wright > w0 + wres_fine*num_wpoints_fine)
            {
                wright = w0 + wres_fine*num_wpoints_fine;
            }

            /*Fine grid calculation.*/
            uint64_t leftid = (2*((wleft-w0)/wres_fine)+1)/2;
            if (leftid >= num_wpoints_fine)
            {
                leftid = num_wpoints_fine - 1;
            }
            uint64_t rightid = (2*((wright-w0)/wres_fine)+1)/2;
            if (rightid >= num_wpoints_fine)
            {
                rightid = num_wpoints_fine - 1;
            }
            uint64_t i;
            for (i=leftid;i<=rightid;++i)
            {
                in.w = w0 + i*wres_fine;
#if defined(DOPPLER)
                fp_t line_shape = doppler_line_shape(in);
#elif defined(LORENTZ)
                fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                fp_t line_shape = ida_voigt_line_shape(in);
#else
                fp_t line_shape = rfm_voigt_line_shape(in);
#endif
                atomicAdd(&(tau_fine[lyr*num_wpoints_fine+i]),
                          snn*n*line_shape);
            }

            /*Coarse grid calculation.*/
            /*Left side.*/
            double cleft = in.line_center - wcutoff;
            if (cleft < w0)
            {
                cleft = w0;
            }
            leftid = (2*((cleft-w0)/wres_coarse)+1)/2;
            if (leftid >= num_wpoints_coarse)
            {
                leftid = num_wpoints_coarse - 1;
            }
            rightid = (2*((wleft-w0)/wres_coarse)+1)/2;
            if (rightid >= num_wpoints_coarse)
            {
                rightid = num_wpoints_coarse - 1;
            }
            for (i=leftid;i<=rightid;++i)
            {
                in.w = w0 + i*wres_coarse;
#if defined(DOPPLER)
                fp_t line_shape = doppler_line_shape(in);
#elif defined(LORENTZ)
                fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                fp_t line_shape = ida_voigt_line_shape(in);
#else
                fp_t line_shape = rfm_voigt_line_shape(in);
#endif
                atomicAdd(&(tau_coarse[lyr*num_wpoints_coarse+i]),
                          snn*n*line_shape);
            }

            /*Right side.*/
            double cright = in.line_center + wcutoff;
            if (cright > w0 + wres_coarse*num_wpoints_coarse)
            {
                cright = w0 + wres_coarse*num_wpoints_coarse;
            }
            leftid = (2*((wright-w0)/wres_coarse)+1)/2;
            if (leftid >= num_wpoints_coarse)
            {
                leftid = num_wpoints_coarse - 1;
            }
            rightid = (2*((cright-w0)/wres_coarse)+1)/2;
            if (rightid >= num_wpoints_coarse)
            {
                rightid = num_wpoints_coarse - 1;
            }
            for (i=leftid;i<=rightid;++i)
            {
                in.w = w0 + i*wres_coarse;
#if defined(DOPPLER)
                fp_t line_shape = doppler_line_shape(in);
#elif defined(LORENTZ)
                fp_t line_shape = lorentz_line_shape(in);
#elif defined(IDA_VOIGT)
                fp_t line_shape = ida_voigt_line_shape(in);
#else
                fp_t line_shape = rfm_voigt_line_shape(in);
#endif
                atomicAdd(&(tau_coarse[lyr*num_wpoints_coarse+i]),
                          snn*n*line_shape);
            }
        }
    }
    return;
}
#endif


/**Host version.*/
void eval_profile_h(int const mol_id,
                    unsigned int const num_lines,
                    uint64_t const num_wpoints_fine,
                    uint64_t const num_wpoints_coarse,
                    double const w0,
                    double const wres_fine,
                    double const wres_coarse,
                    int const num_layers,
                    double const wcutoff,
                    fp_t const * const T,
                    fp_t const * const gamma,
                    fp_t const * const Pshift,
                    fp_t const * const s,
                    fp_t const * const N,
                    fp_t * const tau_fine,
                    fp_t * const tau_coarse,
                    fp_t const fine_factor)
{
    fp_t const molar_mass = get_molar_mass(mol_id);
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
            unsigned int loffset = lyr*num_lines + ltid;
            LineShapeInputs_t in;
            in.line_center = Pshift[loffset];
            fp_t n = N[lyr];
            fp_t temp = T[lyr];
            fp_t snn = s[loffset];
            in.lorentz_hwhm = gamma[loffset];
            in.doppler_hwhm = doppler_hwhm(temp,
                                           molar_mass,
                                           in.line_center);
#ifdef IDA_VOIGT
            in.eta = eta(2.*in.lorentz_hwhm,
                         2.*in.doppler_hwhm);
#endif
/*
            double fcutoff = fine_factor*in.lorentz_hwhm;
*/
            double fcutoff = fine_factor;
            if (fcutoff > wcutoff)
            {
                fcutoff = wcutoff;
            }
            double wleft = in.line_center - fcutoff;
            if (wleft < w0)
            {
                wleft = w0;
            }
            double wright = in.line_center + fcutoff;
            if (wright > w0 + wres_fine*num_wpoints_fine)
            {
                wright = w0 + wres_fine*num_wpoints_fine;
            }

            /*Fine grid calculation.*/
            uint64_t leftid = (2*((wleft-w0)/wres_fine)+1)/2;
            if (leftid >= num_wpoints_fine)
            {
                leftid = num_wpoints_fine - 1;
            }
            uint64_t rightid = (2*((wright-w0)/wres_fine)+1)/2;
            if (rightid >= num_wpoints_fine)
            {
                rightid = num_wpoints_fine - 1;
            }
            uint64_t i;
            for (i=leftid;i<=rightid;++i)
            {
                in.w = w0 + i*wres_fine;
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
                tau_fine[lyr*num_wpoints_fine+i] += snn*n*line_shape;
            }

            /*Coarse grid calculation.*/
            /*Left side.*/
            double cleft = in.line_center - wcutoff;
            if (cleft < w0)
            {
                cleft = w0;
            }
            leftid = (2*((cleft-w0)/wres_coarse)+1)/2;
            if (leftid >= num_wpoints_coarse)
            {
                leftid = num_wpoints_coarse - 1;
            }
            rightid = (2*((wleft-w0)/wres_coarse)+1)/2;
            if (rightid >= num_wpoints_coarse)
            {
                rightid = num_wpoints_coarse - 1;
            }
            for (i=leftid;i<=rightid;++i)
            {
                in.w = w0 + i*wres_coarse;
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
                tau_coarse[lyr*num_wpoints_coarse+i] += snn*n*line_shape;
            }

            /*Right side.*/
            double cright = in.line_center + wcutoff;
            if (cright > w0 + wres_coarse*num_wpoints_coarse)
            {
                cright = w0 + wres_coarse*num_wpoints_coarse;
            }
            leftid = (2*((wright-w0)/wres_coarse)+1)/2;
            if (leftid >= num_wpoints_coarse)
            {
                leftid = num_wpoints_coarse - 1;
            }
            rightid = (2*((cright-w0)/wres_coarse)+1)/2;
            if (rightid >= num_wpoints_coarse)
            {
                rightid = num_wpoints_coarse - 1;
            }
            for (i=leftid;i<=rightid;++i)
            {
                in.w = w0 + i*wres_coarse;
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
                tau_coarse[lyr*num_wpoints_coarse+i] += snn*n*line_shape;
            }
        }
    }
    return;
}
