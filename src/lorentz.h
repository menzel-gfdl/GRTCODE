#ifndef LORENTZ_H_
#define LORENTZ_H_

#include "floating_point_type.h"
#include "line_shape.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t lorentz_line_shape(LineShapeInputs_t const vals);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t lorentz_hwhm(fp_t const P,
                  fp_t const T,
                  float const yself,
                  float const yair,
                  float const n,
                  fp_t const Ps);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t lorentz_fwhm(fp_t const P,
                  fp_t const T,
                  float const yself,
                  float const yair,
                  float const n,
                  fp_t const Ps);


#endif
