#include "debug.h"
#include "floating_point_type.h"
#include "gas_properties.h"
#include "molecules.h"
#include "tips2017.h"


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
fp_t get_molar_mass(int const mol_id)
{
    fp_t res = 0.; /*Molar mass of the molecule.*/
    switch(mol_id)
    {
        case H2O:
            res = 18.010565;
            break;
        case CO2:
            res = 43.98983;
            break;
        case O3:
            res = 47.984745;
            break;
        case N2O:
            res = 44.001062;
            break;
        case CO:
            res = 27.994915;
            break;
        case CH4:
            res = 16.0313;
            break;
        case O2:
            res = 31.98983;
            break;
        case NO:
            res = 29.997989;
            break;
        case SO2:
            res = 63.961901;
            break;
        case NO2:
            res = 45.992904;
            break;
        case NH3:
            res = 17.026549;
            break;
        case HNO3:
            res = 62.995644;
            break;
        case OH:
            res = 17.00274;
            break;
        case HF:
            res = 20.006229;
            break;
        case HCl:
            res = 35.976678;
            break;
        case HBr:
            res = 79.92616;
            break;
        case HI:
            res = 127.912297;
            break;
        case ClO:
            res = 50.963768;
            break;
        case OCS:
            res = 59.966986;
            break;
        case H2CO:
            res = 30.010565;
            break;
        case HOCl:
            res = 51.971593;
            break;
        case N2:
            res = 28.006148;
            break;
        case HCN:
            res = 27.010899;
            break;
        case CH3Cl:
            res = 49.992328;
            break;
        case H2O2:
            res = 34.00548;
            break;
        case C2H2:
            res = 26.01565;
            break;
        case C2H6:
            res = 30.04695;
            break;
        case PH3:
            res = 33.997238;
            break;
        case COF2:
            res = 65.991722;
            break;
        case SF6:
            res = 145.962492;
            break;
        case H2S:
            res = 33.987721;
            break;
        case HCOOH:
            res = 46.00548;
            break;
        case HO2:
            res = 32.997655;
            break;
        case O:
            res = 15.994915;
            break;
        case ClONO2:
            res = 96.956672;
            break;
        case NOp:
            res = 29.997989;
            break;
        case HOBr:
            res = 95.921076;
            break;
        case C2H4:
            res = 28.0313;
            break;
        case CH3OH:
            res = 32.026215;
            break;
        case CH3Br:
            res = 93.941811;
            break;
        case CH3CN:
            res = 41.026549;
            break;
        case CF4:
            res = 87.993616;
            break;
        case C4H2:
            res = 50.01565;
            break;
        case HC3N:
            res = 51.010899;
            break;
        case H2:
            res = 2.01565;
            break;
        case CS:
            res = 43.971036;
            break;
        case SO3:
            res = 79.95682;
            break;
        case C2N2:
            res = 52.006148;
            break;
        case COCl2:
            res = 97.9326199796;
            break;
        case SO:
            res = 48.0644;
            break;
        case C3H4:
            res = 40.0639;
            break;
        case CH3:
            res = 15.035;
            break;
        case CS2:
            res = 76.139;
            break;
        default:
            kernel_err("the molecular with id=%d is not implemented.",
                       mol_id);
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
fp_t Q(int const mol_id,
       fp_t const T,
       int const iso)
{
    fp_t Qt; /*Total internal partition function.*/
    QT(mol_id,
       T,
       iso,
       &Qt);
    return Qt;
}
