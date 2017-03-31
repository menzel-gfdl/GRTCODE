#include "eval_pShift.h"
#include "LineShapeUtils.h"
#include "myreal.h"

#ifdef __NVCC__
/*GPU kernel version.*/

/*---------------------------------------------------------------------------*/
/*Compute the pressure-shift correction of the line position for each
  transition.

  Arguments:
      numLayers [in]      Size of the height dimension for the inputted
                              arrays.
      nL        [in]      Size of the line dimension for the inputted arrays.
      P         [in]      Array of pressures (atm).  This array is stored as
                              [height].
      Vnn       [in]      Array of spectral line transition frequencies
                              (cm^-1).  This array is stored as [line].
      d         [in]      Array of air-broadened pressure shifts at
                              (T=296K,p=1atm) of the line transition
                              frequencies (cm^-1*atm^-1).  This array is
                              stored as [line].
      PShift    [in,out]  Array of pressure-shift corrections of the line
                              positions (cm^-1).  This array is stored as
                              [height][line].
*/
__global__
void eval_pShift(unsigned int const numLayers,
                 unsigned int const nL,
                 REAL_t const * const P,
                 REAL_t const * const Vnn,
                 float const * const d,
                 REAL_t * const PShift)
{
    /*Local variables*/
    unsigned int lyr;
    unsigned int ltid = blockIdx.x*blockDim.x + threadIdx.x;

    if (ltid < nL)
    {
#pragma unroll
        for (lyr=0;lyr<numLayers;++lyr)
        {
            PShift[lyr*nL+ltid] = pressureShiftCorrection(Vnn[ltid],
                                                          d[ltid],
                                                          P[lyr]);
        }
    }

    return;
}

/*---------------------------------------------------------------------------*/

#endif
/*Host only version.*/

/*---------------------------------------------------------------------------*/
/*Compute the pressure-shift correction of the line position for each
  transition.

  Arguments:
      numLayers [in]      Size of the height dimension for the inputted
                              arrays.
      nL        [in]      Size of the line dimension for the inputted arrays.
      P         [in]      Array of pressures (atm).  This array is stored as
                              [height].
      Vnn       [in]      Array of spectral line transition frequencies
                              (cm^-1).  This array is stored as [line].
      d         [in]      Array of air-broadened pressure shifts at
                              (T=296K,p=1atm) of the line transition
                              frequencies (cm^-1*atm^-1).  This array is
                              stored as [line].
      PShift    [in,out]  Array of pressure-shift corrections of the line
                              positions (cm^-1).  This array is stored as
                              [height][line].
*/
void eval_pShift_h(unsigned int const numLayers,
                   unsigned int const nL,
                   REAL_t const * const P,
                   REAL_t const * const Vnn,
                   float const * const d,
                   REAL_t * const PShift)
{
    /*Local variables*/
    unsigned int lyr;
    unsigned int ltid;

    for (ltid=0;ltid<nL;++ltid)
    {
        for (lyr=0;lyr<numLayers;++lyr)
        {
            PShift[lyr*nL+ltid] = pressureShiftCorrection(Vnn[ltid],
                                                          d[ltid],
                                                          P[lyr]);
        }
    }

    return;
}
/*---------------------------------------------------------------------------*/

