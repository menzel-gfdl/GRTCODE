#ifndef LINESHAPEUTILS_H_
#define LINESHAPEUTILS_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t pressureShiftCorrection(fp_t const Vnn,
                             float const d,
                             fp_t const P);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t Snn_partialCorrection(int const molId,
                           int const iso,
                           fp_t const Vnn,
                           float const En,
                           fp_t const Snn_ref);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t Snn_Tcorrection(int const molId,
                     fp_t const T,
                     int const iso,
                     fp_t const Vnn,
                     float const En,
                     fp_t const Snn_partial);


#endif
