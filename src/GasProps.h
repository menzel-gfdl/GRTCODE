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


#endif
