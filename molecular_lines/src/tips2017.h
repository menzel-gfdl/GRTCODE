#ifndef TIPS2017_H_
#define TIPS2017_H_

#include "floating_point_type.h"


int inittips_d(void);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t Q(int const mol_id,
       fp_t const T,
       int const iso
      );


#endif
