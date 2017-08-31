#include <stdint.h>
#include "LineShapeUtils.h"
#include "myreal.h"
#include "pre_eval_Snn.h"
#include "omp.h"

/*---------------------------------------------------------------------------*/
/*Compute the first part of the temperature correction of the line intensities
  for each transition.  These include all terms dependent on the HITRAN
  reference temperature.

  Arguments:
      nL        [in]    Size of the line dimension for the inputted arrays.
      molId     [in]    Molecule id.
      iso       [in]    Array of isotope indexes.  This array is stored as
                            [line].
      Vnn       [in]    Array of spectral line transition frequencies
                            (cm^-1).  This array is stored as [line].
      En        [in]    Array of lower state energies of the transitions
                            (cm^-1).  This array is stored as [line].
      Snn_ref [in,out]  Array of spectral line intensities (cm).  This
                            array is stored as [line].
*/

#ifdef __NVCC__
__global__ void pre_eval_Snn(unsigned int const nL,
                             uint8_t const molId,
                             uint8_t const * const iso,
                             REAL_t const * const Vnn,
                             float const * const En,
                             REAL_t * const Snn_ref)
{
    int ltid = blockIdx.x*blockDim.x + threadIdx.x;

    if (ltid < nL)
    {
        Snn_ref[ltid] = Snn_partialCorrection(molId,
                                              iso[ltid],
                                              Vnn[ltid],
                                              En[ltid],
                                              Snn_ref[ltid]);
    }

    return;
}

#endif

void pre_eval_Snn_h(unsigned int const nL,
                    uint8_t const molId,
                    uint8_t const * const iso,
                    REAL_t const * const Vnn,
                    float const * const En,
                    REAL_t * const Snn_ref)
{
    unsigned int ltid;

#pragma omp parallel for schedule(static) \
                         default(none) \
                         private(ltid)
/*
                         shared(nL,Snn_ref,molId,iso,Vnn,En)
*/
    for (ltid=0;ltid<nL;++ltid)
    {
        Snn_ref[ltid] = Snn_partialCorrection(molId,
                                              iso[ltid],
                                              Vnn[ltid],
                                              En[ltid],
                                              Snn_ref[ltid]);
    }

    return;
}
