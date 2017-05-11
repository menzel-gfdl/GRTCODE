#ifndef SET_EVAL_GAMMA_H_
#define SET_EVAL_GAMMA_H_

#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__global__
void eval_gamma(unsigned int const numLayers,
                unsigned int const nL,
                REAL_t const * const P,
                REAL_t const * const T,
                REAL_t const * const Ps,
                float const * const Yself,
                float const * const Yair,
                float const * const n,
                REAL_t * const Gam);

#endif
void eval_gamma_h(unsigned int const numLayers,
                  unsigned int const nL,
                  REAL_t const * const P,
                  REAL_t const * const T,
                  REAL_t const * const Ps,
                  float const * const Yself,
                  float const * const Yair,
                  float const * const n,
                  REAL_t * const Gam);

#endif

