#include <math.h>
#include "floating_point_type.h"
#include "GaussianFuncs.h"
#include "IdaVoigtFuncs.h"
#include "line_shape.h"
#include "LorentzFuncs.h"


#ifdef __NVCC__
__host__ __device__
#endif
static fp_t f(fp_t const lorFWHM,
              fp_t const gauFWHM);


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
               (1.-vals.eta)*gaussian_line_shape(vals);
}


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
static fp_t f(fp_t const lorFWHM,
              fp_t const gauFWHM)
{
    fp_t const lorFWHM2 = lorFWHM*lorFWHM;   /*lorFWHM^2.*/
    fp_t const lorFWHM3 = lorFWHM2*lorFWHM;  /*lorFWHM^3.*/
    fp_t const lorFWHM4 = lorFWHM2*lorFWHM2; /*lorFWHM^4.*/
    fp_t const lorFWHM5 = lorFWHM3*lorFWHM2; /*lorFWHM^5.*/
    fp_t const gauFWHM2 = gauFWHM*gauFWHM;   /*gauFWHM^2.*/
    fp_t const gauFWHM3 = gauFWHM2*gauFWHM;  /*gauFWHM^3.*/
    fp_t const gauFWHM4 = gauFWHM2*gauFWHM2; /*gauFWHM^4.*/
    fp_t const gauFWHM5 = gauFWHM3*gauFWHM2; /*gauFWHM^5.*/
    fp_t const c1 = 2.69269;                 /*Constant defined in the
                                                   paper.*/
    fp_t const c2 = 2.42843;                 /*Constant defined in the
                                                   paper.*/
    fp_t const c3 = 4.47163;                 /*Constant defined in the
                                                   paper.*/
    fp_t const c4 = 0.07842;                 /*Constant defined in the
                                                   paper.*/

    fp_t f = gauFWHM5 + c1*gauFWHM4*lorFWHM + c2*gauFWHM3*lorFWHM2
             + c3*gauFWHM2*lorFWHM3 + c4*gauFWHM*lorFWHM4 + lorFWHM5;

    f = powf(f,0.2f);
    return f;
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
fp_t eta(fp_t const lorFWHM,
         fp_t const gauFWHM)
{
    fp_t const c1 = 1.36603;    /*Constant defined in the paper.*/
    fp_t const c2 = -0.47719;   /*Constant defined in the paper.*/
    fp_t const c3 = 0.11116;    /*Constant defined in the paper.*/
    fp_t const gam = f(lorFWHM,
                         gauFWHM);
    fp_t const flf = lorFWHM/gam;
    return flf*(c1 + c2*flf + c3*flf*flf);
}
