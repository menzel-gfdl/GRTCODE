#ifndef RFM_VOIGT_H_
#define RFM_VOIGT_H_

#include "floating_point_type.h"
#include "line_shape.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t rfm_voigt_line_shape(LineShapeInputs_t const vals);


#endif
