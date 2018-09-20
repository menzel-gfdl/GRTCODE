#ifndef TIPS_2011_H_
#define TIPS_2011_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__host__ __device__
#endif
void QT(int const molNum, /* HITRAN molecule ID number */
        fp_t const T,       /* temperature in K */
        int const iso,    /* isotope code (HITRAN INDEX) */
        float * const gsi,    /* state independent nuclear degeneracyfactor */
        fp_t * const Qt);    /* Total Internal Partition Function */

#ifdef __NVCC__
__host__
#endif
int initTIPS_d();

#endif
