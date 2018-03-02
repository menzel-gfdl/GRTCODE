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



#ifdef foo

/*For a given molecule, set the partial pressure at each time, latitude,
  longitude, and height.

  Arguments:
      val   [in]      Molecular concentration in pressure units (ppmv).
      PS    [in,out]  Array of partial pressures (atm).  This array is stored
                          as [time][lat][lon][molecule][height].
      P     [in]      Array of atmospheric pressures (atm).  This array is
                          stored as [time][lat][lon][height].
      molId [in]      Id of the molecule whose partial pressure is being
                          calculated.
      ntime [in]      Size of the time dimension for the pressure arrays.
      nlat  [in]      Size of the latitude dimension for the pressure
                          arrays.
      nlon  [in]      Size of the longitude dimension for the pressure
                          arrays.
      nmol  [in]      Size of the molecule dimension for the pressure
                          arrays.
      nlvl  [in]      Size of the height dimension for the pressure arrays.
*/
void setGlobalPartialPres(double const val,
                          REAL_t * const PS,
                          REAL_t const * const P,
                          unsigned int const molId,
                          size_t const ntime,
                          size_t const nlat,
                          size_t const nlon,
                          size_t const nmol,
                          size_t const nlvl)
{
    /*Local variables*/
    unsigned int itr;                         /*Loop variable.*/
    unsigned int lat;                         /*Loop variable.*/
    unsigned int lon;                         /*Loop variable.*/
    unsigned int time;                        /*Loop variable.*/
    size_t ps_off;                            /*Array offset for partial
                                                  pressure.*/
    size_t p_off;                             /*Array offset for pressure.*/
    const size_t ps_off_mol = (molId-1)*nlvl; /*Used to calculate the offset
                                                  for the given molecule.*/

    /*Loop through the arrays.*/
    for (time=0;time<ntime;++time)
    {
        for (lat=0;lat<nlat;++lat)
        {
            for (lon=0;lon<nlon;++lon)
            {
                /*Calculate the offsets.*/
                ps_off = time*nlat*nlon*nmol*nlvl + lat*nlon*nmol*nlvl +
                             lon*nmol*nlvl + ps_off_mol;
                p_off = time*(nlat*nlon*nlvl) + lat*(nlon*nlvl) + lon*nlvl;

                /*Calculate the partial pressures.*/
                for (itr=0;itr<nlvl;++itr)
                {
                    PS[ps_off+itr] = (val/1.e6)*P[p_off+itr];
                }
            }
        }
    }

    return;
}

/*---------------------------------------------------------------------------*/
/*Calculate the number density of the air from the ideal gas law.

  Arguments:
      P_atm [in]  Pressure (atm).
      T_k   [in]  Temperature (K).

  Return:
      (P_atm*6.022E23)/(T_k*82.057338)  Number density (cm^-3).
*/
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t idealGasNumberDensity(REAL_t const P_atm,
                             REAL_t const T_k)
{
    /*Local variables*/
    REAL_t const R = 82.057338; /*Gas constant (cm^3*atm*K^-1*mol^-1).*/
    REAL_t const AV = 6.022E23; /*Avagadro's number (1/mol).*/

    return (P_atm*AV)/(R*T_k);
}

/*---------------------------------------------------------------------------*/
/*Calculate the number density for a molecular species.

  Arguments:
      N           [in,out]  Number density (cm^-3).  This array is stored as
                                [time][lat][lon][molecule][height].
      PartialPres [in]      Partial pressure (atm).  This array is stored as
                                [time][lat][lon][molecule][height].
      T           [in]      Temperature (K).  This array is stored as
                                [time][lat][lon][height].
      hitranMolId [in]      Molecule id from the hitran database.
      ntime       [in]      Size of the time dimension for the partial
                                pressure, temperature, and number density
                                arrays.
      nlat        [in]      Size of the latitude dimension for the partial
                                pressure, temperature, and number density
                                arrays.
      nlon        [in]      Size of the longitude dimension for the partial
                                pressure, temperature, and number density
                                arrays.
      nmol        [in]      Size of the molecule dimension for the pressure
                                and number density arrays.
      nlvl        [in]      Size of the height dimension for the partial
                                pressure, temperature, and number density
                                arrays.
*/
void setGlobalNumberDensity(REAL_t * const N,
                            REAL_t const * const PartialPres,
                            REAL_t const * const T,
                            unsigned int const hitranMolId,
                            size_t const ntime,
                            size_t const nlat,
                            size_t const nlon,
                            size_t const nmol,
                            size_t const nlvl)
{
    /*Local variables*/
    unsigned int itr;                              /*Loop variable.*/
    unsigned int lat;                              /*Loop variable.*/
    unsigned int lon;                              /*Loop variable.*/
    unsigned int time;                             /*Loop variable.*/
    size_t off;                                    /*Array offset.*/
    size_t n_off;                                  /*Array offset.*/
    const size_t n_off_mol = (hitranMolId-1)*nlvl; /*Used to calculate the
                                                       offset for the inputted
                                                       molecule.*/

    /*Loop through the arrays.*/
    for(time=0;time<ntime;++time)
    {
        for(lat=0;lat<nlat;++lat)
        {
            for(lon=0;lon<nlon;++lon)
            {
                /*Calculate the array offsets.*/
                off = time*nlat*nlon*nlvl + lat*nlon*nlvl + lon*nlvl;
                n_off = time*nlat*nlon*nmol*nlvl + lat*nlon*nmol*nlvl +
                        lon*nmol*nlvl + n_off_mol;

                /*Calculate the number densities using the ideal gas law.*/
                for (itr=0;itr<nlvl;++itr)
                {
                    N[n_off+itr] = idealGasNumberDensity(PartialPres[n_off+itr],
                                                         T[off+itr]);
                }
            }
        }
    }

    return;
}
#endif
