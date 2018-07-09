#ifndef IDA_VOIGT_H_
#define IDA_VOIGT_H_

#include "floating_point_type.h"
#include "line_shape.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t ida_voigt_line_shape(LineShapeInputs_t const vals);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t eta(fp_t const lorentz_fwhm,
         fp_t const doppler_fwhm);


#endif
