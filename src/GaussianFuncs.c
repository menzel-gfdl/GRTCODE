#include <math.h>
#include "GaussianFuncs.h"
#include "line_shape.h"
#include "myreal.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef SQRT_LN2
#define SQRT_LN2 0.83255461115
#endif

/*---------------------------------------------------------------------------*/
/*Calculate a gaussian line shape function.

  Arguments:
      vals [in]  Structure containing line shape input values.

  Return:
      Gaussian line shape function (cm).
*/
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gaussian_line_shape(LineShapeInputs_t const vals)
{
    /*Local variables*/
    REAL_t const x = vals.freq - vals.lineCenter;
    REAL_t const alphad = vals.gauHWHM/SQRT_LN2;

    return expf(-x*x/(alphad*alphad))/(alphad*sqrtf(M_PI));
}

/*---------------------------------------------------------------------------*/
/*Calculate "alphad" for a gaussian.  alphad is defined so that:

  Gaussian half-width at half max = sqrt(ln(2))*alphad

  Arguments:
      T  [in]  Temperature (K).
      M  [in]  Molar mass (g/mol).
      v0 [in]  Spectral line center frequency (cm^-1).

  Return:
      alphad for the gaussian (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gaussian_alphad(REAL_t const T,
                       REAL_t const M,
                       REAL_t const v0)
{
    /*Local variables*/
    REAL_t const m = M/6.023E23;    /*Mass (g).*/
    REAL_t const c = 2.99792458E10; /*Speed of light (cm/s).*/
    REAL_t const kb = 1.380658E-16; /*Boltzmann constant (erg/K)*/

    return v0*sqrtf((2*kb*T)/(m*c*c));
}

/*---------------------------------------------------------------------------*/
/*Calculate the half-width at half-max for a gaussian.

  Arguments:
      T  [in]  Temperature (K).
      M  [in]  Molar mass (g/mol).
      v0 [in]  Spectral line center frequency (cm^-1).

  Return:
      Half-width at half max for the gaussian (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gaussian_hwhm(REAL_t const T,
                     REAL_t const M,
                     REAL_t const v0)
{
    return SQRT_LN2*gaussian_alphad(T,M,v0);
}

/*---------------------------------------------------------------------------*/
/*Calculate the full-width at half-max for a gaussian.

  Arguments:
      T  [in]  Temperature (K).
      M  [in]  Molar mass (g/mol).
      v0 [in]  Spectral line center frequency (cm^-1).

  Return:
      Full-width at half max for the gaussian (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gaussian_fwhm(REAL_t const T,
                     REAL_t const M,
                     REAL_t const v0)
{
    return 2.0*gaussian_hwhm(T,M,v0);
}

/*---------------------------------------------------------------------------*/
