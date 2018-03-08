#include "debug.h"
#include "floating_point_type.h"
#include "molecules.h"
#include "TIPS_2011.h"


/*Return the molar mass of the molecule specified by the input molecule
  id.

  Arguments:
      molId [in]  Molecule id.

  Return:
      Mass of the molecule (g/mol).
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t getMolarMass(int const molId)
{
    fp_t res; /*Molar mass of the molecule.*/
    switch(molId)
    {
        case H2O:
            res = 18.01528;
            break;
        case CO2:
            res = 44.01;
            break;
        case O3:
            res = 48.;
            break;
        case N2O:
            res = 44.013;
            break;
        case CO:
            res = 28.01;
            break;
        case CH4:
            res = 16.04;
            break;
        case O2:
            res = 32.;
            break;
        default:
            kernel_err("the molecular with molId=%d is not implemented.",
                       molId);
            break;
    }
    return res;
}


/*Calculate the total internal partition function for the input molecule
  using the method located in TIPS_2011.c.

  Arguments:
      moldId [in]  A molecule id.
      T      [in]  Temperature (K).
      iso    [in]  Isotope index.

  Return:
      The total internal partition function.
*/
#ifdef __NVCC__
__host__ __device__
#endif
fp_t Q(int const molId,
       fp_t const T,
       int const iso)
{
    float gsi; /*State independent nuclear degeneracy factor.*/
    fp_t Qt; /*Total internal partition function.*/
    QT(molId,
       T,
       iso,
       &gsi,
       &Qt);
    return Qt;
}
