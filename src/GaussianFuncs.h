#ifndef GAUSSIANFUNCS_H_
#define GAUSSIANFUNCS_H_

#include "floating_point_type.h"
#include "line_shape.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t gaussian_line_shape(LineShapeInputs_t const vals);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t gaussian_alphad(fp_t const T,
                     fp_t const M,
                     fp_t const v0);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t gaussian_hwhm(fp_t const T,
                   fp_t const M,
                   fp_t const v0);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t gaussian_fwhm(fp_t const T,
                   fp_t const M,
                   fp_t const v0);


#endif
