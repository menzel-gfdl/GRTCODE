#include <stdint.h>
#include "eval_Snn_correction.h"
#include "LineShapeUtils.h"
#include "myreal.h"
#include "omp.h"

/*---------------------------------------------------------------------------*/
/*Compute the temperature correction of the line intensities for each
  transition.

  Arguments:
      numLayers   [in]      Size of the height dimension for the inputted
                                arrays.
      nL          [in]      Size of the line dimension for the inputted
                                arrays.
      molId       [in]      Molecule id.
      T           [in]      Array of temperatures (K).  This array is stored
                                as [height].
      iso         [in]      Array of isotope indexes.  This array is stored
                                as [line].
      Vnn         [in]      Array of spectral line transition frequencies
                                (cm^-1).  This array is stored as [line].
      En          [in]      Array of lower state energies of the transitions
                                (cm^-1).  This array is stored as [line].
      Snn_partial [in]      Array of partially corrected spectral line
                                intensities (cm).  This array is stored
                                as [line].
      S           [in,out]  Array of corrected spectral line intensities (cm).
                                This array is stored as [height][line].
*/

#ifdef __NVCC__
__global__ void eval_Snn_correction(unsigned int const numLayers,
                                    unsigned int const nL,
                                    uint8_t const molId,
                                    REAL_t const * const T,
                                    uint8_t const * const iso,
                                    REAL_t const * const Vnn,
                                    float const * const En,
                                    REAL_t const * const Snn_partial,
                                    REAL_t * const S)
{
    /*Local variables*/
    unsigned int lyr;
    unsigned int ltid = blockIdx.x*blockDim.x + threadIdx.x;

    if (ltid < nL)
    {
#pragma unroll
        for (lyr=0;lyr<numLayers;++lyr)
        {
            S[lyr*nL+ltid] = Snn_Tcorrection(molId,
                                             T[lyr],
                                             iso[ltid],
                                             Vnn[ltid],
                                             En[ltid],
                                             Snn_partial[ltid]);
        }
    }

    return;
}

#endif

void eval_Snn_correction_h(unsigned int const numLayers,
                           unsigned int const nL,
                           uint8_t const molId,
                           REAL_t const * const T,
                           uint8_t const * const iso,
                           REAL_t const * const Vnn,
                           float const * const En,
                           REAL_t const * const Snn_partial,
                           REAL_t * const S)
{
    /*Local variables*/
    unsigned int lyr;
    unsigned int ltid;

#pragma omp parallel for schedule(dynamic) \
                         collapse(2) \
                         default(none) \
                         private(ltid,lyr)
/*
                         shared(nL,numLayers,S,molId,T,iso, \
                                Vnn,En,Snn_partial)
*/
    for (lyr=0;lyr<numLayers;++lyr)
    {
        for (ltid=0;ltid<nL;++ltid)
        {
            S[lyr*nL+ltid] = Snn_Tcorrection(molId,
                                             T[lyr],
                                             iso[ltid],
                                             Vnn[ltid],
                                             En[ltid],
                                             Snn_partial[ltid]);
        }
    }

    return;
}
