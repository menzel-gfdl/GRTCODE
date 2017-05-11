#ifndef SET_RFMVOIGTFUNCS_H_
#define SET_RFMVOIGTFUNCS_H_

#include "line_shape.h"
#include "myreal.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Function prototypes.*/

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t rfm_voigt_line_shape(LineShapeInputs_t const vals);

#endif

