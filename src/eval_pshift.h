#ifndef EVAL_PSHIFT_H_
#define EVAL_PSHIFT_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__global__ void eval_pshift(int const num_layers,
                            unsigned int const num_lines,
                            fp_t const * const P,
                            fp_t const * const vnn,
                            float const * const d,
                            fp_t * const pshift);
#endif


void eval_pshift_h(int const num_layers,
                   unsigned int const num_lines,
                   fp_t const * const P,
                   fp_t const * const vnn,
                   float const * const d,
                   fp_t * const pshift);


#endif
