#include "omp.h"
#include "floating_point_type.h"
#include "line_shape_utils.h"
#include "pre_eval_snn.h"


/*Compute the first part of the temperature correction of the line intensities
  for each transition.  These include all terms dependent on the HITRAN
  reference temperature.

  Arguments:
      num_lines        [in]    Size of the line dimension for the inputted arrays.
      mol_id     [in]    Molecule id.
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
__global__ void pre_eval_snn(unsigned int const num_lines,
                             int const mol_id,
                             int const * const iso,
                             fp_t const * const vnn,
                             float const * const en,
                             fp_t * const snn_ref)
{
    int const ltid = blockIdx.x*blockDim.x + threadIdx.x;
    if (ltid < num_lines)
    {
        snn_ref[ltid] = snn_partial_correction(mol_id,
                                               iso[ltid],
                                               vnn[ltid],
                                               en[ltid],
                                               snn_ref[ltid]);
    }
    return;
}
#endif


void pre_eval_snn_h(unsigned int const num_lines,
                    int const mol_id,
                    int const * const iso,
                    fp_t const * const vnn,
                    float const * const en,
                    fp_t * const snn_ref)
{
    unsigned int ltid;

#pragma omp parallel for schedule(static) \
                         default(none) \
                         private(ltid)
    for (ltid=0;ltid<num_lines;++ltid)
    {
        snn_ref[ltid] = snn_partial_correction(mol_id,
                                               iso[ltid],
                                               vnn[ltid],
                                               en[ltid],
                                               snn_ref[ltid]);
    }
    return;
}
