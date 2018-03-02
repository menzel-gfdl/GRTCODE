#ifndef EVAL_PROFILE_H_
#define EVAL_PROFILE_H_

#include "floating_point_type.h"


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
                             fp_t * const tau);
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
                    fp_t * const tau);


#endif
