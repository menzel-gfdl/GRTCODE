/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef SET_PARSENETCDFRADIATION_H_
#define SET_PARSENETCDFRADIATION_H_

#include "myreal.h"

/*Notes on units in the input NetCDF files:
    time = Days since 01/01/1982 (days).
    lat = Latitude grid points (degrees eastward from ???).
    lon = longitude grid points (degrees northward from ???).
    pfull = Array of pressure layers (hPa).
    phalf = Array of pressure layer interfaces (hPa). This array is one bigger
            than pfull.

    The input arrays are stored as (time,pressure,lat,lon).
*/
typedef struct radiationInputFields_t
{
    float *RH2O;   /*Layer water vapor mole fractions (mol/mol) or mixing
                         ratios (kg/kg).*/
    float *RCO2;   /*Layer carbon dioxide mole fractions (mol/mol).*/
    float *QO3;    /*Layer ozone mole fractions (mol/mol) or mixing
                         ratios (kg/kg).*/
    float *RN2O;   /*Layer nitrous oxide mole fractions (mol/mol).*/
    float *RCO;    /*Layer carbon monoxide mole fractions (mol/mol).*/
    float *RCH4;   /*Layer methane mole fractions (mol/mol).*/
    float *RO2;    /*Layer oxygen mole fractions (mol/mol).*/
    float *PRESSM; /*Layer pressures (Pa).*/
    float *TEMP;   /*Layer Temperatures (K).*/
    float *DELTAZ; /*Layer thicknesses (m).*/
    float *TSURF; /*Surface temperature (K).*/
    float *TLEV; /*Level temperatures (K).*/
    float *EMIS; /*Surface emissivity.*/
    size_t nlat;   /*Number of latitude grid points.*/
    size_t nlon;   /*Number of longitude grid points.*/
    size_t npfull; /*Number of pressure layers.*/
    size_t nphalf; /*Number of pressure layer interfaces.*/
    size_t ntime;  /*Number of time grid points.*/
} radiationInputFields_t;

/*Notes on units in the radiation output fields.:
    The output arrays are stored as (time,lat,lon,pressure) or
    (time,lat,lon,molecule,pressure).
*/
typedef struct radiationOutputFields_t
{
    REAL_t *N;      /*Layer molecular number densities (1/cm^3).*/
    REAL_t *P;      /*Layer pressures (atm).*/
    REAL_t *T;      /*Layer temperatures (K).*/
    REAL_t *DELTAZ; /*Layer thicknesses (cm).*/
    REAL_t *PS;     /*Layer parital pressures (atm).*/
    REAL_t *TSURF; /*Surface temperature (K).*/
    REAL_t *TLEV; /*Level temperatures (K).*/
    REAL_t *EMIS; /*Surface emissivity.*/
    size_t nlat;    /*Number of latitude grid points.*/
    size_t nlon;    /*Number of longitude grid points.*/
    size_t npfull;  /*Number of pressure layers.*/
    size_t nphalf;  /*NUmber of pressure layer interfaces.*/
    size_t ntime;   /*NUmber of time grid points.*/
} radiationOutputFields_t;

int readRfmipFieldsFromFile(char fname[],
                            radiationInputFields_t* in);

int setOutputFieldsFromRfmip(radiationInputFields_t *in,
                             radiationOutputFields_t *out);

int readGfdlFieldsFromFile(char fname[],
                           radiationInputFields_t* in);

int setOutputFieldsFromGfdl(radiationInputFields_t *in,
                            radiationOutputFields_t *out);

int radiationOutputFieldsMalloc(radiationOutputFields_t* out);

int radiationOutputFieldsFree(radiationOutputFields_t* out);

int getAndSetAtmosFieldsFromFile(char fname[],
                                 char *f_format,
                                 radiationOutputFields_t* out);

#ifndef SKIPMAIN
int test(char fname[],
         char *f_format);

int test2(char fname[],
          char *f_format);
#endif

#endif
