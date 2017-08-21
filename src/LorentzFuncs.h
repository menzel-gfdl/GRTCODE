#ifndef SET_LORENTZFUNCS_H_
#define SET_LORENTZFUNCS_H_

#include "line_shape.h"
#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t lorentz_line_shape(LineShapeInputs_t const vals);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t lorentz_hwhm(REAL_t const P,
                    REAL_t const T,
                    float const Yself,
                    float const Yair,
                    float const n,
                    REAL_t const Ps);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t lorentz_fwhm(REAL_t const P,
                    REAL_t const T,
                    float const Yself,
                    float const Yair,
                    float const n,
                    REAL_t const Ps);

#endif
