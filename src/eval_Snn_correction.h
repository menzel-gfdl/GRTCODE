#ifndef SET_EVAL_SNN_CORRECTION_H_
#define SET_EVAL_SNN_CORRECTION_H_

#include <stdint.h>
#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__global__ void eval_Snn_correction(unsigned int const numLayers,
                                    unsigned int const nL,
                                    uint8_t const molId,
                                    REAL_t const * const T,
                                    uint8_t const * const iso,
                                    REAL_t const * const Vnn,
                                    float const * const En,
                                    REAL_t const * const Snn_partial,
                                    REAL_t * const S);
#endif

void eval_Snn_correction_h(unsigned int const numLayers,
                           unsigned int const nL,
                           uint8_t const molId,
                           REAL_t const * const T,
                           uint8_t const * const iso,
                           REAL_t const * const Vnn,
                           float const * const En,
                           REAL_t const * const Snn_partial,
                           REAL_t * const S);
#endif
