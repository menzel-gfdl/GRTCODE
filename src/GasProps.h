#ifndef SET_GASPROPS_H_
#define SET_GASPROPS_H_

#include <stdint.h>
#include <stdlib.h>
#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t getMolarMass(int const hitranMolId);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t Q(uint8_t const molId,
         REAL_t const T,
         uint8_t const iso);

void setGlobalPartialPres(double const val,
                          REAL_t * const PS,
                          REAL_t const * const P,
                          unsigned int const molId,
                          size_t const ntime,
                          size_t const nlat,
                          size_t const nlon,
                          size_t const nmol,
                          size_t const nlvl);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t idealGasNumberDensity(REAL_t const P_atm,
                             REAL_t const T_k);

void setGlobalNumberDensity(REAL_t * const N,
                            REAL_t const * const PartialPres,
                            REAL_t const * const T,
                            unsigned int const hitranMolId,
                            size_t const ntime,
                            size_t const nlat,
                            size_t const nlon,
                            size_t const nmol,
                            size_t const nlvl);

#endif

