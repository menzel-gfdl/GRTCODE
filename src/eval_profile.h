#ifndef SET_EVAL_PROFILE_H_
#define SET_EVAL_PROFILE_H_

#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__global__ void eval_profile(unsigned int const molId,
                             unsigned int const nL,
                             int const nF,
                             REAL_t const loWn,
                             REAL_t const resolution,
                             unsigned int const numLayers,
                             unsigned int const breadth,
                             REAL_t const * const T,
                             REAL_t const * const Gam,
                             REAL_t const * const PShift,
                             REAL_t const * const S,
                             REAL_t const * const tauU_d,
                             REAL_t * const out);
#endif

void eval_profile_h(unsigned int const molId,
                    unsigned int const nL,
                    int const nF,
                    REAL_t const loWn,
                    REAL_t const resolution,
                    unsigned int const numLayers,
                    unsigned int const breadth,
                    REAL_t const * const T,
                    REAL_t const * const Gam,
                    REAL_t const * const PShift,
                    REAL_t const * const S,
                    REAL_t const * const tauU_d,
                    REAL_t const * const pathlength_d,
                    REAL_t * const out);
#endif
