#ifndef SET_FLUX_H_
#define SET_FLUX_H_

#include "myreal.h"

#ifdef __NVCC__
__global__ void calcFlux(unsigned int const nF,
                         unsigned int const numLayers,
                         REAL_t * const fluxDown,
                         REAL_t * const fluxUp,
                         REAL_t const * const T,
                         REAL_t const Tsurf,
                         REAL_t const * const tau,
                         REAL_t const w,
                         REAL_t const res,
                         REAL_t const emissivity,
                         REAL_t const * const TLEV);
#endif

void calcFlux_h(unsigned int const nF,
                unsigned int const numLayers,
                REAL_t * const fluxDown,
                REAL_t * const fluxUp,
                REAL_t const * const T,
                REAL_t const Tsurf,
                REAL_t const * const tau,
                REAL_t const w,
                REAL_t const res,
                REAL_t const emissivity,
                REAL_t const * const TLEV);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t planckFunc(REAL_t const T,
                  REAL_t const w);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t effectivePlanck(REAL_t const Tcenter,
                       REAL_t const Tedge,
                       REAL_t const w,
                       REAL_t const tau);

void sum_fluxes(unsigned int const nF,
                unsigned int const numLevels,
                REAL_t const * const fluxes,
                REAL_t * const fluxes_accumulated,
                REAL_t const res);

/*
Old versions.

#ifdef __NVCC__
__global__ void calcFlux(unsigned int const nF,
                         unsigned int const numLayers,
                         REAL_t * const fluxDown,
                         REAL_t * const fluxUp,
                         REAL_t const * const T,
                         REAL_t const Tsurf,
                         REAL_t const * const tau,
                         REAL_t const w,
                         REAL_t const res);
#endif

void calcFlux_h(unsigned int const nF,
                unsigned int const numLayers,
                REAL_t * const fluxDown,
                REAL_t * const fluxUp,
                REAL_t const * const T,
                REAL_t const Tsurf,
                REAL_t const * const tau,
                REAL_t const w,
                REAL_t const res);
*/

#endif
