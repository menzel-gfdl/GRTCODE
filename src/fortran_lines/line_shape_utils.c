#include <math.h>
#include "floating_point_type.h"
#include "gas_properties.h"
#include "HITRAN_constants.h"
#include "line_shape_utils.h"


/*Calculate the pressure-shift correction of the line position.  See
  equation A13 from:

  Rothman, L. S. et. al (1998). J. Quant. Spectrosc. Radiat. Transfer. 60,
      665-710.

  Arguments:
      Vnn [in]  Spectral line transition frequency (cm^-1).
      d   [in]  Air-broadened pressure shift at (T=296K,p=1atm) of the line
                    transition frequency (cm^-1*atm^-1).
      P   [in]  Pressure (atm).

  Return:
      Pressure-shift correction of the line position (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t pressure_shift_correction(fp_t const vnn,
                               float const d,
                               fp_t const P)
{
    return (vnn + ((fp_t)d)*P);
}


/*Part of the temperature correction of the line intensity.  Includes all
  terms dependent on the reference temperature.  See equation A11 from:

  Rothman, L. S. et. al (1998). J. Quant. Spectrosc. Radiat. Transfer. 60,
      665-710.

  Arguments:
      molId   [in]  Molecule id.
      iso     [in]  Isotope index.
      Vnn     [in]  Spectral line transition frequency (cm^-1).
      En      [in]  Lower state energy of the transition (cm^-1).
      Snn_ref [in]  Spectral line intensity at reference
                        temperature 296K (cm).

  Return:
      A partial correction to the line intensity (cm).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t snn_partial_correction(int const mol_id,
                            int const iso,
                            fp_t const vnn,
                            float const en,
                            fp_t const snn_ref)
{
    return (snn_ref*Q(mol_id,TREF,iso))/
               (exp(-c2*en/TREF)*(1-exp(-c2*(vnn/TREF))));
}


/*Part of the temperature correction of the line intensity due to the current
  temperature.  See equation A11 from:

  Rothman, L. S. et. al (1998). J. Quant. Spectrosc. Radiat. Transfer. 60,
      665-710.

  Arguments:
      molId       [in]  Molecule id.
      T           [in]  Temperature (K).
      iso         [in]  Isotope index.
      Vnn         [in]  Spectral line transition frequency (cm^-1).
      En          [in]  Lower state energy of the transition (cm^-1).
      Snn_partial [in]  Partial correction to the line intensity (cm).
  Return:
      Temperature correction of the line intensity (cm).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t snn_T_correction(int const mol_id,
                      fp_t const T,
                      int const iso,
                      fp_t const vnn,
                      float const en,
                      fp_t const snn_partial)
{
    return (snn_partial/Q(mol_id,T,iso))*(((fp_t)1)-exp(-c2*vnn/T))*
               exp(-c2*en/T);
}
