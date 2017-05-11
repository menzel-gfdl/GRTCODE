#ifndef SET_PRE_EVAL_SNN_H_
#define SET_PRE_EVAL_SNN_H_

#include <stdint.h>
#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__global__
void pre_eval_Snn(unsigned int const nL,
                  uint8_t const molId,
                  uint8_t const * const iso,
                  REAL_t const * const Vnn,
                  float const * const En,
                  REAL_t * const Snn_ref);

#endif
void pre_eval_Snn_h(unsigned int const nL,
                    uint8_t const molId,
                    uint8_t const * const iso,
                    REAL_t const * const Vnn,
                    float const * const En,
                    REAL_t * const Snn_ref);

#endif

