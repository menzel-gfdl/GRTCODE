#ifndef LW_FLUX_H_
#define LW_FLUX_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__global__ void calc_lw_flux(unsigned int const nF,
                             int const numLayers,
                             fp_t * const fluxDown,
                             fp_t * const fluxUp,
                             fp_t const * const T,
                             fp_t const Tsurf,
                             fp_t const * const tau,
                             fp_t const w,
                             fp_t const res,
                             fp_t const emissivity,
                             fp_t const * const TLEV);
#endif


void calc_lw_flux_h(unsigned int const nF,
                    int const numLayers,
                    fp_t * const fluxDown,
                    fp_t * const fluxUp,
                    fp_t const * const T,
                    fp_t const Tsurf,
                    fp_t const * const tau,
                    fp_t const w,
                    fp_t const res,
                    fp_t const emissivity,
                    fp_t const * const TLEV);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t planck_func(fp_t const T,
                 fp_t const w);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t effective_planck(fp_t const Tcenter,
                      fp_t const Tedge,
                      fp_t const w,
                      fp_t const tau);


void integrate_fluxes(unsigned int const nF,
                      int const numLevels,
                      fp_t const * const fluxes,
                      fp_t * const fluxes_accumulated,
                      fp_t const res);


#endif
