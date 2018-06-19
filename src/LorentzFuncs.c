#include <math.h>
#include "floating_point_type.h"
#include "HitranConstants.h"
#include "line_shape.h"
#include "LorentzFuncs.h"

#ifndef M_1_PI
#define M_1_PI 0.31830988618379067154
#endif


/*Calculate the normalized line shape function assuming a Lorentz profile.
  See equation A14 from:

  Rothman, L. S. et. al (1998). J. Quant. Spectrosc. Radiat. Transfer. 60,
      665-710.

  Arguments:
      vals [in]  Structure containing line shape input values.

  Return:
      Value of the normalized line shape function assuming a Lorentz
          profile (cm).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t lorentz_line_shape(LineShapeInputs_t const vals)
{
    fp_t const gam2 = vals.lorHWHM*vals.lorHWHM;
    fp_t const del = vals.freq - vals.lineCenter;
    return ((fp_t)M_1_PI)*(vals.lorHWHM/(gam2+(del*del)));
}


/*Calculate the pressure broadened line half-width.  See equation A12 from:

  Rothman, L. S. et. al (1998). J. Quant. Spectrosc. Radiat. Transfer. 60,
      665-710.

  Arguments:
      P     [in]  Pressure (atm).
      T     [in]  Temperature (K).
      Yself [in]  Self-broadened halfwidth at half maximum (cm^-1*atm^-1).
      Yair  [in]  Air-broadened halfwidth at half maximum (cm^-1*atm^-1).
      n     [in]  Coefficient of temperature dependence of the air-broadened
                      halfwidth at half maximum.
      Ps    [in]  Partial pressure (atm).

  Return:
      Pressure broadened line half-width (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t lorentz_hwhm(fp_t const P,
                  fp_t const T,
                  float const Yself,
                  float const Yair,
                  float const n,
                  fp_t const Ps)
{
    return pow((((fp_t)TREF)/T),(fp_t)n)*
               ((((fp_t)Yair)*(P-Ps)) + (((fp_t)Yself)*Ps));
}


/*Calculate the pressure broadened line full-width.  See equation A12 from:

  Arguments:
      P     [in]  Pressure (atm).
      T     [in]  Temperature (K).
      Yself [in]  Self-broadened halfwidth at half maximum (cm^-1*atm^-1).
      Yair  [in]  Air-broadened halfwidth at half maximum (cm^-1*atm^-1).
      n     [in]  Coefficient of temperature dependence of the air-broadened
                      halfwidth at half maximum.
      Ps    [in]  Partial pressure (atm).

  Return:
      Pressure broadened line full-width (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t lorentz_fwhm(fp_t const P,
                  fp_t const T,
                  float const Yself,
                  float const Yair,
                  float const n,
                  fp_t const Ps)
{
    return 2.0*lorentz_hwhm(P,T,Yself,Yair,n,Ps);
}
