#ifndef EVAL_GAMMA_H_
#define EVAL_GAMMA_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__global__ void eval_gamma(int const numLayers,
                           unsigned int const nL,
                           fp_t const * const P,
                           fp_t const * const T,
                           fp_t const * const Ps,
                           float const * const Yself,
                           float const * const Yair,
                           float const * const n,
                           fp_t * const Gam);
#endif


void eval_gamma_h(int const numLayers,
                  unsigned int const nL,
                  fp_t const * const P,
                  fp_t const * const T,
                  fp_t const * const Ps,
                  float const * const Yself,
                  float const * const Yair,
                  float const * const n,
                  fp_t * const Gam);


#endif
