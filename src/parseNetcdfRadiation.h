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
    float *RH2O;   /*Layer water vapor mixing ratios (kg/kg).*/
    float *RCO2;   /*Layer carbon dioxide mixing ratios (kg/kg).*/
    float *QO3;    /*Layer ozone mixing ratios (kg/kg).*/
    float *RN2O;   /*Layer nitrous oxide mixing ratios (kg/kg).*/
    float *RCO;    /*Layer carbon monoxide mixing ratios (kg/kg).*/
    float *RCH4;   /*Layer methane mixing ratios (kg/kg).*/
    float *RO2;    /*Layer oxygen mixing ratios (kg/kg).*/
    float *DPFLUX; /*Radiation flus layer thicknesses [(dP/dz)*delta_z] (hPa).*/
    float *PRESSM; /*Layer pressures (Pa).*/
    float *TEMP;   /*Layer Temperatures (K).*/
    float *DELTAZ; /*Layer thicknesses [delta_z] (m).*/
    size_t nlat;   /*Number of latitude grid points.*/
    size_t nlon;   /*Number of longitude grid points.*/
    size_t npfull; /*Number of pressure layers.*/
    size_t nphalf; /*Number of pressure layer interfaces.*/
    size_t ntime;  /*Number of time grid points.*/
} radiationInputFields_t;

/*Notes on units in the radiation output fields.:
    The output arrays are stored as (time,lat,lon,pressure).
*/
typedef struct radiationOutputFields_t
{
    REAL_t* N;      /*Layer molecular number densities (1/cm^3).*/
    REAL_t* P;      /*Layer pressures (atm).*/
    REAL_t* T;      /*Layer temperatures (K).*/
    REAL_t* DELTAZ; /*Layer thicknesses (cm).*/
    REAL_t* PS;     /*Layer parital pressures (atm).*/
    size_t nlat;    /*Number of latitude grid points.*/
    size_t nlon;    /*Number of longitude grid points.*/
    size_t npfull;  /*Number of pressure layers.*/
    size_t nphalf;  /*NUmber of pressure layer interfaces.*/
    size_t ntime;   /*NUmber of time grid points.*/
} radiationOutputFields_t;

int radiationInputFieldsMalloc(radiationInputFields_t* in);

int radiationOutputFieldsMalloc(radiationOutputFields_t* out);

int radiationInputFieldsFree(radiationInputFields_t* in);

int radiationOutputFieldsFree(radiationOutputFields_t* out);

int readInputFieldsFromFile(char fname[],
                            radiationInputFields_t* in);

REAL_t getNumberDensity(const REAL_t rh2o,
                        const REAL_t dpflux,
                        const int hitranMolId);

REAL_t getPartialPres(REAL_t rh2o,
                      REAL_t pressm,
                      const REAL_t molarMass);

int setOutputFields(radiationInputFields_t *in,
                    radiationOutputFields_t *out);

int getAndSetAtmosFieldsFromFile(char fname[],
                                 char *f_format,
                                 radiationOutputFields_t* out);

int test(char fname[]);

#endif
