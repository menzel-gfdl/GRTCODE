#ifndef EVAL_SNN_CORRECTION_H_
#define EVAL_SNN_CORRECTION_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__global__ void eval_snn_correction(int const num_layers,
                                    unsigned int const num_lines,
                                    int const mol_id,
                                    fp_t const * const T,
                                    int const * const iso,
                                    fp_t const * const vnn,
                                    float const * const en,
                                    fp_t const * const snn_partial,
                                    fp_t * const s);
#endif


void eval_snn_correction_h(int const num_layers,
                           unsigned int const num_lines,
                           int const mol_id,
                           fp_t const * const T,
                           int const * const iso,
                           fp_t const * const vnn,
                           float const * const en,
                           fp_t const * const snn_partial,
                           fp_t * const s);


#endif
