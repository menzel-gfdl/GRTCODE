#ifndef GASPROPS_H_
#define GASPROPS_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t getMolarMass(int const molId);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t Q(int const molId,
       fp_t const T,
       int const iso);


#ifdef foo
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


#endif
