#ifndef EVAL_GAMMA_H_
#define EVAL_GAMMA_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__global__ void eval_gamma(int const num_layers,
                           unsigned int const num_lines,
                           fp_t const * const P,
                           fp_t const * const T,
                           fp_t const * const Ps,
                           float const * const yself,
                           float const * const yair,
                           float const * const n,
                           fp_t * const gamma);
#endif


void eval_gamma_h(int const num_layers,
                  unsigned int const num_lines,
                  fp_t const * const P,
                  fp_t const * const T,
                  fp_t const * const Ps,
                  float const * const yself,
                  float const * const yair,
                  float const * const n,
                  fp_t * const gamma);


#endif
