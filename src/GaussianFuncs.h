#ifndef SET_GAUSSIANFUNCS_H_
#define SET_GAUSSIANFUNCS_H_

#include "line_shape.h"
#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gaussian_line_shape(LineShapeInputs_t const vals);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gaussian_alphad(REAL_t const T,
                       REAL_t const M,
                       REAL_t const v0);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gaussian_hwhm(REAL_t const T,
                     REAL_t const M,
                     REAL_t const v0);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gaussian_fwhm(REAL_t const T,
                     REAL_t const M,
                     REAL_t const v0);

#endif
