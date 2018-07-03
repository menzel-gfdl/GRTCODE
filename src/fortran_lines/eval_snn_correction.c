#ifdef _OPENMP
#include "omp.h"
#endif
#include "eval_snn_correction.h"
#include "floating_point_type.h"
#include "line_shape_utils.h"


/*Compute the temperature correction of the line intensities for each
  transition.

  Arguments:
      num_layers   [in]      Size of the height dimension for the inputted
                                arrays.
      num_lines          [in]      Size of the line dimension for the inputted
                                arrays.
      mol_id       [in]      Molecule id.
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
__global__ void eval_snn_correction(int const num_layers,
                                    unsigned int const num_lines,
                                    int const mol_id,
                                    fp_t const * const T,
                                    int const * const iso,
                                    fp_t const * const vnn,
                                    float const * const en,
                                    fp_t const * const snn_partial,
                                    fp_t * const s)
{
    int lyr;
    unsigned int const ltid = blockIdx.x*blockDim.x + threadIdx.x;
    if (ltid < num_lines)
    {
#pragma unroll
        for (lyr=0;lyr<num_layers;++lyr)
        {
            s[lyr*num_lines+ltid] = snn_T_correction(mol_id,
                                                     T[lyr],
                                                     iso[ltid],
                                                     vnn[ltid],
                                                     en[ltid],
                                                     snn_partial[ltid]);
        }
    }
    return;
}
#endif


void eval_snn_correction_h(int const num_layers,
                           unsigned int const num_lines,
                           int const mol_id,
                           fp_t const * const T,
                           int const * const iso,
                           fp_t const * const vnn,
                           float const * const en,
                           fp_t const * const snn_partial,
                           fp_t * const s)
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
            s[lyr*num_lines+ltid] = snn_T_correction(mol_id,
                                                     T[lyr],
                                                     iso[ltid],
                                                     vnn[ltid],
                                                     en[ltid],
                                                     snn_partial[ltid]);
        }
    }
    return;
}
