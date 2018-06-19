#ifndef IDAVOIGTFUNCS_H_
#define IDAVOIGTFUNCS_H_

#include "floating_point_type.h"
#include "line_shape.h"


#ifdef __NVCC__
__host__ __device__
#endif
fp_t ida_voigt_line_shape(LineShapeInputs_t const vals);


#ifdef __NVCC__
__host__ __device__
#endif
fp_t eta(fp_t const lorFWHM,
         fp_t const gauFWHM);


#endif
