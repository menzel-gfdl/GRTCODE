#ifndef PRE_EVAL_SNN_H_
#define PRE_EVAL_SNN_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__global__ void pre_eval_snn(unsigned int const num_lines,
                             int const mol_id,
                             int const * const iso,
                             fp_t const * const vnn,
                             float const * const en,
                             fp_t * const snn_ref);

#endif


void pre_eval_snn_h(unsigned int const num_lines,
                    int const mol_id,
                    int const * const iso,
                    fp_t const * const vnn,
                    float const * const en,
                    fp_t * const snn_ref);


#endif
