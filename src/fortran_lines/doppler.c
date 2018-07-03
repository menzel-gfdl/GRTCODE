#include <math.h>
#include "doppler.h"
#include "floating_point_type.h"
#include "line_shape.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef SQRT_LN2
#define SQRT_LN2 0.83255461115
#endif


/*Calculate a doppler line shape function.

  Arguments:
      vals [in]  Structure containing line shape input values.

  Return:
      Doppler line shape function (cm).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t doppler_line_shape(LineShapeInputs_t const vals)
{
    fp_t const x = vals.w - vals.line_center;
    fp_t const alphad = vals.doppler_hwhm/SQRT_LN2;
    return expf(-x*x/(alphad*alphad))/(alphad*sqrtf(M_PI));
}


/*Calculate "alphad" for a doppler profile.  alphad is defined so that:

  Doppler half-width at half max = sqrt(ln(2))*alphad

  Arguments:
      T  [in]  Temperature (K).
      M  [in]  Molar mass (g/mol).
      v0 [in]  Spectral line center frequency (cm^-1).

  Return:
      alphad for the doppler profile (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t doppler_alphad(fp_t const T,
                    fp_t const M,
                    fp_t const v0)
{
    fp_t const m = M/6.023E23;    /*Mass (g).*/
    fp_t const c = 2.99792458E10; /*Speed of light (cm/s).*/
    fp_t const kb = 1.380658E-16; /*Boltzmann constant (erg/K)*/
    return v0*sqrtf((2*kb*T)/(m*c*c));
}


/*Calculate the half-width at half-max for a doppler profile.

  Arguments:
      T  [in]  Temperature (K).
      M  [in]  Molar mass (g/mol).
      v0 [in]  Spectral line center frequency (cm^-1).

  Return:
      Half-width at half max for the doppler profile (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t doppler_hwhm(fp_t const T,
                  fp_t const M,
                  fp_t const v0)
{
    return SQRT_LN2*doppler_alphad(T,M,v0);
}


/*Calculate the full-width at half-max for a doppler profile.

  Arguments:
      T  [in]  Temperature (K).
      M  [in]  Molar mass (g/mol).
      v0 [in]  Spectral line center frequency (cm^-1).

  Return:
      Full-width at half max for the doppler profile (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t doppler_fwhm(fp_t const T,
                  fp_t const M,
                  fp_t const v0)
{
    return 2.0*doppler_hwhm(T,M,v0);
}
