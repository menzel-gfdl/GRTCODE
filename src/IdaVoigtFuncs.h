#ifndef SET_IDAVOIGTFUNCS_H_
#define SET_IDAVOIGTFUNCS_H_

#include "line_shape.h"
#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t ida_voigt_line_shape(LineShapeInputs_t const vals);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t eta(REAL_t const lorFWHM,
           REAL_t const gauFWHM);

#endif
