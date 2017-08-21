#ifndef SET_EVAL_PSHIFT_H_
#define SET_EVAL_PSHIFT_H_

#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__global__ void eval_pShift(unsigned int const numLayers,
                            unsigned int const nL,
                            REAL_t const * const P,
                            REAL_t const * const Vnn,
                            float const * const d,
                            REAL_t * const PShift);
#endif

void eval_pShift_h(unsigned int const numLayers,
                   unsigned int const nL,
                   REAL_t const * const P,
                   REAL_t const * const Vnn,
                   float const * const d,
                   REAL_t * const PShift);
#endif
