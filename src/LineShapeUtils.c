#include <math.h>
#include "GasProps.h"
#include "HitranConstants.h"
#include "LineShapeUtils.h"
#include "myreal.h"

/*---------------------------------------------------------------------------*/
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
REAL_t pressureShiftCorrection(REAL_t const Vnn,
                               float const d,
                               REAL_t const P)
{
    return (Vnn + ((REAL_t)d)*P);
}

/*---------------------------------------------------------------------------*/
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
REAL_t Snn_partialCorrection(uint8_t const molId,
                             uint8_t const iso,
                             REAL_t const Vnn,
                             float const En,
                             REAL_t const Snn_ref)
{
    return (Snn_ref*Q(molId,TREF,iso))/
               (exp(-c2*En/TREF)*(1-exp(-c2*(Vnn/TREF))));
}

/*---------------------------------------------------------------------------*/
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
REAL_t Snn_Tcorrection(uint8_t const molId,
                       REAL_t const T,
                       uint8_t const iso,
                       REAL_t const Vnn,
                       float const En,
                       REAL_t const Snn_partial)
{
    return (Snn_partial/Q(molId,T,iso))*(((REAL_t)1)-exp(-c2*Vnn/T))*
               exp(-c2*En/T);
}

/*---------------------------------------------------------------------------*/

