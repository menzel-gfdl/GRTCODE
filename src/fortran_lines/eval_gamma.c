#ifdef _OPENMP
#include <omp.h>
#endif
#include "eval_gamma.h"
#include "floating_point_type.h"
#include "lorentz.h"


/*Compute the pressure broadened line halfwidth for each transition.

  Arguments:
      numLayers [in]      Size of the height dimension for the inputted
                              arrays.
      nL        [in]      Size of the line dimension for the inputted arrays.
      P         [in]      Array of pressures (atm).  This array is stored as
                              [height].
      T         [in]      Array of temperatures (K).  This array is stored as
                              [height].
      Ps        [in]      Array of partial pressures (atm).  This array is
                              stored as [height].
      Yself     [in]      Array of self-broadened halfwidth at half
                              maximum (cm^-1*atm^-1).  This array is stored as
                              [line].
      Yair      [in]      Array of air-broadened halfwidth at half
                              maximum (cm^-1*atm^-1).  This array is stored as
                              [line].
      n         [in]      Array of coefficients of temperature dependence of
                              the air-broadened halfwidth at half maximum.
                              This array is stored as [line].
      Gam       [in,out]  Array of pressure broadened line halfwidths (cm^-1).
                              This array is stored as [height][line].
*/
#ifdef __NVCC__
__global__
void eval_gamma(int const num_layers,
                unsigned int const num_lines,
                fp_t const * const P,
                fp_t const * const T,
                fp_t const * const Ps,
                float const * const yself,
                float const * const yair,
                float const * const n,
                fp_t * const gamma)
{
    int lyr;
    unsigned int const ltid = blockIdx.x*blockDim.x + threadIdx.x;
    if (ltid < num_lines)
    {
#pragma unroll
        for (lyr=0;lyr<num_layers;++lyr)
        {
            gamma[lyr*num_lines+ltid] = lorentz_hwhm(P[lyr],
                                                     T[lyr],
                                                     yself[ltid],
                                                     yair[ltid],
                                                     n[ltid],
                                                     Ps[lyr]);
        }
    }
    return;
}
#endif


void eval_gamma_h(int const num_layers,
                  unsigned int const num_lines,
                  fp_t const * const P,
                  fp_t const * const T,
                  fp_t const * const Ps,
                  float const * const yself,
                  float const * const yair,
                  float const * const n,
                  fp_t * const gamma)
{
    int lyr;
    unsigned int ltid;

#pragma omp parallel for schedule(static) \
                         collapse(2) \
                         default(none) \
                         private(ltid,lyr)
    for (lyr=0;lyr<num_layers;++lyr)
    {
        for(ltid=0;ltid<num_lines;++ltid)
        {
            gamma[lyr*num_lines+ltid] = lorentz_hwhm(P[lyr],
                                                     T[lyr],
                                                     yself[ltid],
                                                     yair[ltid],
                                                     n[ltid],
                                                     Ps[lyr]);
        }
    }
    return;
}
