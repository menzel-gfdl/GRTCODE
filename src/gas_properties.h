#ifndef GAS_PROPERTIES_H_
#define GAS_PROPERTIES_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t get_molar_mass(int const mol_id);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t Q(int const mol_id,
       fp_t const T,
       int const iso);


#endif
