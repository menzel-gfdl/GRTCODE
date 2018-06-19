#ifndef EVAL_SNN_CORRECTION_H_
#define EVAL_SNN_CORRECTION_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__global__ void eval_Snn_correction(int const numLayers,
                                    unsigned int const nL,
                                    int const molId,
                                    fp_t const * const T,
                                    int const * const iso,
                                    fp_t const * const Vnn,
                                    float const * const En,
                                    fp_t const * const Snn_partial,
                                    fp_t * const S);
#endif


void eval_Snn_correction_h(int const numLayers,
                           unsigned int const nL,
                           int const molId,
                           fp_t const * const T,
                           int const * const iso,
                           fp_t const * const Vnn,
                           float const * const En,
                           fp_t const * const Snn_partial,
                           fp_t * const S);


#endif
