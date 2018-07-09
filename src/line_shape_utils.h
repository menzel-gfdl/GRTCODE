#ifndef LINE_SHAPE_UTILS_H_
#define LINE_SHAPE_UTILS_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t pressure_shift_correction(fp_t const vnn,
                               float const d,
                               fp_t const P);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t snn_partial_correction(int const mol_id,
                            int const iso,
                            fp_t const vnn,
                            float const en,
                            fp_t const snn_ref);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t snn_T_correction(int const mol_id,
                      fp_t const T,
                      int const iso,
                      fp_t const vnn,
                      float const en,
                      fp_t const snn_partial);


#endif
