/*
    GRTCODE is a GPU-able Radiative Transfer Code
    Copyright (C) 2016  Garrett Wright

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; version 2.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
    USA.
*/

#include <math.h>
#include "voigt.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif

/* Implementation of method:
   Ida, T., Ando, M. & Toraya, H. (2000) J. Appl. Cryst. 33, 1311-1316.
   Authors claim accurate to within 1%.
*/

/*---------------------------------------------------------------------------*/
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
static REAL_t f(REAL_t const lorFWHM,
                REAL_t const gauFWHM)
{
    /*Local variables*/
    REAL_t f;                                  /*Parameter needed for the
                                                   pseudovoigt profile.*/
    const REAL_t gauFWHM2 = gauFWHM*gauFWHM;   /*gauFWHM^2.*/
    const REAL_t gauFWHM3 = gauFWHM2*gauFWHM;  /*gauFWHM^3.*/
    const REAL_t gauFWHM4 = gauFWHM2*gauFWHM2; /*gauFWHM^4.*/
    const REAL_t gauFWHM5 = gauFWHM3*gauFWHM2; /*gauFWHM^5.*/
    const REAL_t lorFWHM2 = lorFWHM*lorFWHM;   /*lorFWHM^2.*/
    const REAL_t lorFWHM3 = lorFWHM2*lorFWHM;  /*lorFWHM^3.*/
    const REAL_t lorFWHM4 = lorFWHM2*lorFWHM2; /*lorFWHM^4.*/
    const REAL_t lorFWHM5 = lorFWHM3*lorFWHM2; /*lorFWHM^5.*/
    const REAL_t c1 = 2.69269;                 /*Constant defined in the
                                                   paper.*/
    const REAL_t c2 = 2.42843;                 /*Constant defined in the
                                                   paper.*/
    const REAL_t c3 = 4.47163;                 /*Constant defined in the
                                                   paper.*/
    const REAL_t c4 = 0.07842;                 /*Constant defined in the
                                                   paper.*/

    f = gauFWHM5 + c1*gauFWHM4*lorFWHM + c2*gauFWHM3*lorFWHM2
        + c3*gauFWHM2*lorFWHM3 + c4*gauFWHM*lorFWHM4 + lorFWHM5;

    f = powf(f,0.2f);

    return f;
}

/*---------------------------------------------------------------------------*/
/*Calculate the parameter described in equation 9 in:

  Ida, T., Ando, M. & Toraya H. (2000). J. Appl. Cryst. 33, 1311-1316.

  Arguments:
      lorFWHM [in]  Full-width at half maximum for a Lorentz profile (cm^-1).
      f       [in]  Paremeter described in equation 10 of the paper (cm^-1).

  Return:
      Parameter needed for the pseudovoigt profile.
*/
#ifdef __NVCC__
__host__ __device__
#endif
static REAL_t eta_(REAL_t const lorFWHM,
                   REAL_t const f)
{
    /*Local variables*/
    REAL_t eta;                   /*Parameter needed for the pseudovoigt
                                      profile.*/
    const REAL_t flf = lorFWHM/f; /*Ratio of inputs.*/
    const REAL_t c1 = 1.36603;    /*Constant defined in the paper.*/
    const REAL_t c2 = -0.47719;   /*Constant defined in the paper.*/
    const REAL_t c3 = 0.11116;    /*Constant defined in the paper.*/

    eta = flf*(c1 + c2*flf + c3*flf*flf);

    return eta;
}

/*---------------------------------------------------------------------------*/
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
REAL_t eta(REAL_t lorFWHM,
           REAL_t gauFWHM)
{
    return eta_(lorFWHM,f(lorFWHM,gauFWHM));
}

/*---------------------------------------------------------------------------*/
/*Calculate "alphad" for a gaussian.  alphad is defined so that:

  Gaussian half-width at half max = sqrt(ln(2))*alphad

  Arguments:
      T  [in]  Temperature (K).
      M  [in]  Molar mass (g/mol).
      v0 [in]  Spectral line frequency (cm^-1).

  Return:
      alphad for the gaussian (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gauAlphad(REAL_t const T,
                 REAL_t const M,
                 REAL_t const v0)
{
    /*Local variables*/
    const REAL_t m = M/6.023E23;    /*Mass (g).*/
    const REAL_t c = 2.99792458E10; /*Speed of light (cm/s).*/
    const REAL_t kb = 1.380658E-16; /*Boltzmann constant (erg/K)*/

    const REAL_t alphad = v0*sqrtf((2*kb*T)/(m*c*c));

    return alphad;
}

/*---------------------------------------------------------------------------*/
/*Calculate the full-width at half maximum for a gaussian profile.

  Arguments:
      T  [in]  Temperature(K).
      M  [in]  Molar mass (g/mol).
      v0 [in]  Spectral line frequency (cm^-1).

  Return:
      Full width at half maximum (cm^-1).
*/
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gauFWHM(REAL_t const T,
               REAL_t const M,
               REAL_t const v0)
{
    /*Local variables*/
    const REAL_t sqrtln2 = 0.83255461115;

    const REAL_t HWHM = sqrtln2*gauAlphad(T,M,v0);

    return 2*HWHM;
}

/*---------------------------------------------------------------------------*/
/*Calculate a gaussian line shape function.

  Arguments:
      v      [in]  Frequency (cm^-1).
      v0     [in]  Spectral line frequency (cm^-1).
      alphad [in]  alphad value for the gaussian (cm^-1).

  Return:
      Gaussian line shape function (cm).
*/
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gauKernel(REAL_t const v,
                 REAL_t const v0,
                 REAL_t const alphad)
{
    /*Local variables*/
    const REAL_t x = (v-v0); /*Frequency difference.*/

    return expf(-x*x/(alphad*alphad))/(alphad*sqrtf(M_PI));
}

/*---------------------------------------------------------------------------*/
/*Calculate the pseudovoigt line shape function.  See equation 6 from:

  Ida, T., Ando, M. & Toraya H. (2000). J. Appl. Cryst. 33, 1311-1316.

  Arguments:
      eta  [in]  Mixing parameter for the pseudovoigt profile.
      lory [in]  Lorentz line shape function (cm).
      gauy [in]  Gaussian line shape function (cm).

  Return:
      Pseudovoigt line shape function (cm).
*/
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t pseudoVoigt(REAL_t const eta,
                   REAL_t const lory,
                   REAL_t const gauy)
{
    return eta*lory + (1.-eta)*gauy;
}

/*---------------------------------------------------------------------------*/

