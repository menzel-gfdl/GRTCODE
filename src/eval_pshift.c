#ifdef _OPENMP
#include "omp.h"
#endif
#include "eval_pshift.h"
#include "floating_point_type.h"
#include "line_shape_utils.h"


/*Compute the pressure-shift correction of the line position for each
  transition.

  Arguments:
      num_layers [in]      Size of the height dimension for the inputted
                              arrays.
      num_lines        [in]      Size of the line dimension for the inputted arrays.
      P         [in]      Array of pressures (atm).  This array is stored as
                              [height].
      vnn       [in]      Array of spectral line transition frequencies
                              (cm^-1).  This array is stored as [line].
      d         [in]      Array of air-broadened pressure shifts at
                              (T=296K,p=1atm) of the line transition
                              frequencies (cm^-1*atm^-1).  This array is
                              stored as [line].
      pshift    [in,out]  Array of pressure-shift corrections of the line
                              positions (cm^-1).  This array is stored as
                              [height][line].
*/
#ifdef __NVCC__
__global__ void eval_pshift(int const num_layers,
                            unsigned int const num_lines,
                            fp_t const * const P,
                            fp_t const * const vnn,
                            float const * const d,
                            fp_t * const pshift)
{
    int lyr;
    unsigned int const ltid = blockIdx.x*blockDim.x + threadIdx.x;
    if (ltid < num_lines)
    {
#pragma unroll
        for (lyr=0;lyr<num_layers;++lyr)
        {
            pshift[lyr*num_lines+ltid] = pressure_shift_correction(vnn[ltid],
                                                                   d[ltid],
                                                                   P[lyr]);
        }
    }
    return;
}
#endif


void eval_pshift_h(int const num_layers,
                   unsigned int const num_lines,
                   fp_t const * const P,
                   fp_t const * const vnn,
                   float const * const d,
                   fp_t * const pshift)
{
    int lyr;
    unsigned int ltid;

#pragma omp parallel for schedule(static) \
                         collapse(2) \
                         default(none) \
                         private(ltid,lyr)
    for (lyr=0;lyr<num_layers;++lyr)
    {
        for (ltid=0;ltid<num_lines;++ltid)
        {
            pshift[lyr*num_lines+ltid] = pressure_shift_correction(vnn[ltid],
                                                                   d[ltid],
                                                                   P[lyr]);
        }
    }
    return;
}
