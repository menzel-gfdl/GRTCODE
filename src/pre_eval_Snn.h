#ifndef PRE_EVAL_SNN_H_
#define PRE_EVAL_SNN_H_

#include "floating_point_type.h"

#ifdef __NVCC__
__global__ void pre_eval_Snn(unsigned int const nL,
                             int const molId,
                             int const * const iso,
                             fp_t const * const Vnn,
                             float const * const En,
                             fp_t * const Snn_ref);

#endif

void pre_eval_Snn_h(unsigned int const nL,
                    int const molId,
                    int const * const iso,
                    fp_t const * const Vnn,
                    float const * const En,
                    fp_t * const Snn_ref);

#endif
