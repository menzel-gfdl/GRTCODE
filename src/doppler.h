#ifndef DOPPLER_H_
#define DOPPLER_H_

#include "floating_point_type.h"
#include "line_shape.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t doppler_line_shape(LineShapeInputs_t const vals);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t doppler_alphad(fp_t const T,
                    fp_t const M,
                    fp_t const v0);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t doppler_hwhm(fp_t const T,
                  fp_t const M,
                  fp_t const v0);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t doppler_fwhm(fp_t const T,
                  fp_t const M,
                  fp_t const v0);


#endif
