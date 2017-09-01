#ifndef SET_FLUX_H_
#define SET_FLUX_H_

#include "myreal.h"

#ifdef __NVCC__
__global__ void calcFlux(unsigned int const nF,
                         unsigned int const numLayers,
                         REAL_t * const fluxOut,
                         REAL_t const * const T,
                         REAL_t const Tsurf,
                         REAL_t const * const tau,
                         REAL_t const w,
                         REAL_t const res);
#endif

void calcFlux_h(unsigned int const nF,
                unsigned int const numLayers,
                REAL_t * const fluxOut,
                REAL_t const * const T,
                REAL_t const Tsurf,
                REAL_t const * const tau,
                REAL_t const w,
                REAL_t const res);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t planckFunc(REAL_t const T,
                  REAL_t const w);

#endif
