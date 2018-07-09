#include <math.h>
#include "doppler.h"
#include "floating_point_type.h"
#include "ida_voigt.h"
#include "line_shape.h"
#include "lorentz.h"


/*Calculate the parameter described in equation 10 in:

  Ida, T., Ando, M. & Toraya H. (2000). J. Appl. Cryst. 33, 1311-1316.

  Arguments:
      lorFWHM [in]  Full-width at half maximum for a Lorentz profile (cm^-1).
      gauFWHM [in]  Full-width at half maximum for a Gaussian profile (cm^-1).

  Return:
      Parameter needed for the pseudovoigt profile (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
static fp_t f(fp_t const lorentz_fwhm,
              fp_t const doppler_fwhm)
{
    fp_t const lorentz_fwhm2 = lorentz_fwhm*lorentz_fwhm; /*lorentz_fwhm^2.*/
    fp_t const lorentz_fwhm3 = lorentz_fwhm2*lorentz_fwhm; /*lorentz_fwhm^3.*/
    fp_t const lorentz_fwhm4 = lorentz_fwhm2*lorentz_fwhm2; /*lorentz_fwhm^4.*/
    fp_t const lorentz_fwhm5 = lorentz_fwhm3*lorentz_fwhm2; /*lorentz_fwhm^5.*/
    fp_t const doppler_fwhm2 = doppler_fwhm*doppler_fwhm; /*doppler_fwhm^2.*/
    fp_t const doppler_fwhm3 = doppler_fwhm2*doppler_fwhm; /*doppler_fwhm^3.*/
    fp_t const doppler_fwhm4 = doppler_fwhm2*doppler_fwhm2; /*doppler_fwhm^4.*/
    fp_t const doppler_fwhm5 = doppler_fwhm3*doppler_fwhm2; /*doppler_fwhm^5.*/
    fp_t const c1 = 2.69269; /*Constant defined in the paper.*/
    fp_t const c2 = 2.42843; /*Constant defined in the paper.*/
    fp_t const c3 = 4.47163; /*Constant defined in the paper.*/
    fp_t const c4 = 0.07842; /*Constant defined in the paper.*/
    fp_t f = doppler_fwhm5 + c1*doppler_fwhm4*lorentz_fwhm
             + c2*doppler_fwhm3*lorentz_fwhm2 + c3*doppler_fwhm2*lorentz_fwhm3
             + c4*doppler_fwhm*lorentz_fwhm4 + lorentz_fwhm5;
    f = powf(f,0.2f);
    return f;
}


/*Calculate the pseudovoigt line shape function.  See equation 6 from:

  Ida, T., Ando, M. & Toraya H. (2000). J. Appl. Cryst. 33, 1311-1316.

  Arguments:
      vals [in]  Structure containing line shape input values.

  Return:
      Pseudovoigt line shape function (cm).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t ida_voigt_line_shape(LineShapeInputs_t const vals)
{
    return vals.eta*lorentz_line_shape(vals) +
               (1.-vals.eta)*doppler_line_shape(vals);
}


/*Calculate the parameter described in equation 9 in:

  Ida, T., Ando, M. & Toraya H. (2000). J. Appl. Cryst. 33, 1311-1316.

  Arguments:
      lorFWHM [in]  Full-width at half maximum for a Lorentz profile (cm^-1).
      gauFWHM [in]  Full-width at half maximum for a Gaussian profile (cm^-1).

  Return:
      Parameter needed for the pseudovoigt profile.
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t eta(fp_t const lorentz_fwhm,
         fp_t const doppler_fwhm)
{
    fp_t const c1 = 1.36603; /*Constant defined in the paper.*/
    fp_t const c2 = -0.47719; /*Constant defined in the paper.*/
    fp_t const c3 = 0.11116; /*Constant defined in the paper.*/
    fp_t const gam = f(lorentz_fwhm,
                       doppler_fwhm);
    fp_t const flf = lorentz_fwhm/gam;
    return flf*(c1 + c2*flf + c3*flf*flf);
}
