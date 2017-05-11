#ifndef SET_LINESHAPEUTILS_H_
#define SET_LINESHAPEUTILS_H_

#include <stdint.h>
#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t pressureShiftCorrection(REAL_t const Vnn,
                               float const d,
                               REAL_t const P);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t Snn_partialCorrection(uint8_t const molId,
                             uint8_t const iso,
                             REAL_t const Vnn,
                             float const En,
                             REAL_t const Snn_ref);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t Snn_Tcorrection(uint8_t const molId,
                       REAL_t const T,
                       uint8_t const iso,
                       REAL_t const Vnn,
                       float const En,
                       REAL_t const Snn_partial);

#endif

