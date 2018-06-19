#ifndef EVAL_PSHIFT_H_
#define EVAL_PSHIFT_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__global__ void eval_pShift(int const numLayers,
                            unsigned int const nL,
                            fp_t const * const P,
                            fp_t const * const Vnn,
                            float const * const d,
                            fp_t * const PShift);
#endif


void eval_pShift_h(int const numLayers,
                   unsigned int const nL,
                   fp_t const * const P,
                   fp_t const * const Vnn,
                   float const * const d,
                   fp_t * const PShift);


#endif
