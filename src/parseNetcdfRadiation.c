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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <netcdf.h>
#include "parseNetcdfRadiation.h"

static const int NUM_MOLS = 7;
static REAL_t const PaToAtm = 9.86923E-6; /*Pa to Atm conversion factor.*/
static REAL_t const MToCm = 100.0; /*Meters to cm conversion factor.*/
static REAL_t const DryAirM = 29.0; /*Dry air molar mass (g/mol)*/
static REAL_t const WaterM = 18.02; /*Water vapor molar mass (g/mol)*/
static REAL_t const OzoneM = 48.0; /*Ozone molar mass (g/mol).*/

/*---------------------------------------------------------------------------*/
/*Handle errors by printing an error message and exiting with a
 *non-zero status. */
#define NCERR(e) {fprintf(stderr, "Error: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}

/*---------------------------------------------------------------------------*/
/*Read in the atmospheric fields from a "rfmip" formatted NetCDF file.*/
int readRfmipFieldsFromFile(char fname[],
                            radiationInputFields_t* in)
{
    /*Local variables*/
    int retval;
    int ncid;
    int varid;
    int dimid;
    size_t i;
    size_t j;
    size_t k;
    size_t zoffset;
    size_t loffset;
    float dp;
    float g = 9.80665; /*(m/s^2)*/
    float R = 8.3144598; /*(m^3*Pa)/(K*mol)*/
    float M = 0.02897; /*kg/mol*/
    float Mkg = M/6.0221409E23; /*kg*/

    /*Open the inputted file. */
    if ((retval = nc_open(fname,
                          NC_NOWRITE,
                          &ncid)))
    {
        NCERR(retval);
    }

    /*Get the number of "profiles" (lat-lon points) from the NetCDF
      dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "profile",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &(in->nlat))))
    {
        NCERR(retval);
    }
    in->nlon = 1;

    /*Get the pressure layer midpoint bounds from the NetCDF
      dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "layer",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &(in->npfull))))
    {
        NCERR(retval);
    }

    /*Get the pressure layer interface bounds from the NetCDF
      dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "level",
                                &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &(in->nphalf))))
    {
        NCERR(retval);
    }

    /*Get the number of experiments (time points) from the NetCDF
      dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "expt",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &(in->ntime))))
    {
        NCERR(retval);
    }

    /*Make sure that the number of pressure layer interfaces is one greater
      than the number of pressure layers.*/
    assert(in->nphalf == in->npfull+1);

    /*Null initialize the radiation input fields.*/
    in->RH2O = NULL;
    in->RCO2 = NULL;
    in->QO3 = NULL;
    in->RN2O = NULL;
    in->RCO = NULL;
    in->RCH4 = NULL;
    in->RO2 = NULL;
    in->PRESSM = NULL;
    in->TEMP = NULL;
    in->DELTAZ = NULL;
    in->TSURF = NULL;
    in->TLEV = NULL;
    in->EMIS = NULL;
    in->N = NULL;

    /*Malloc space for the radiation input fields.*/
    in->RH2O = (float *)malloc((in->ntime)*(in->nlat)*(in->nlon)*(in->npfull)*sizeof(float));
    in->RCO2 = (float *)malloc((in->ntime)*sizeof(float));
    in->QO3 = (float *)malloc((in->ntime)*(in->nlat)*(in->nlon)*(in->npfull)*sizeof(float));
    in->RN2O = (float *)malloc((in->ntime)*sizeof(float));
    in->RCO = (float *)malloc((in->ntime)*sizeof(float));
    in->RCH4 = (float *)malloc((in->ntime)*sizeof(float));
    in->RO2 = (float *)malloc((in->ntime)*sizeof(float));
    in->PRESSM = (float *)malloc((in->nlat)*(in->nlon)*(in->npfull)*sizeof(float));
    in->TEMP = (float *)malloc((in->ntime)*(in->nlat)*(in->nlon)*(in->npfull)*sizeof(float));
    in->DELTAZ = (float *)malloc((in->ntime)*(in->nlat)*(in->nlon)*(in->npfull)*sizeof(float));
    in->TSURF = (float *)malloc((in->ntime)*(in->nlat)*(in->nlon)*sizeof(float));
    in->TLEV = (float *)malloc((in->ntime)*(in->nlat)*(in->nlon)*(in->nphalf)*sizeof(float));
    in->EMIS = (float *)malloc((in->nlat)*(in->nlon)*sizeof(float));
    in->N = (float *)malloc((in->nlat)*(in->nlon)*(in->npfull)*sizeof(float));

    /*Make sure that the mallocs succeeded.*/
    if (in->RH2O == NULL || in->RCO2 == NULL || in->QO3 == NULL ||
            in->RN2O == NULL || in->RCO == NULL || in->RCH4 == NULL ||
            in->RO2 == NULL || in->PRESSM == NULL || in->TEMP == NULL ||
            in->DELTAZ == NULL || in->TSURF == NULL || in->TLEV == NULL ||
            in->EMIS == NULL || in->N == NULL)
    {
        fprintf(stderr,
                "Error(radiationInputFieldsMalloc): malloc failed for the"
                " radiation input fields.\n");
        exit(EXIT_FAILURE);
    }

    /*Read in the water vapor mole fraction in air (mol/mol) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "water_vapor",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->RH2O)))
    {
        NCERR(retval);
    }

    /*Read in the carbon dioxide mole fraction in air (mol/mol) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "carbon_dioxide_GM",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->RCO2)))
    {
        NCERR(retval);
    }

    /*Read in the ozone mole fraction in air (mol/mol) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "ozone",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->QO3)))
    {
        NCERR(retval);
    }

    /*Read in the nitrous oxide mole fraction in air (mol/mol) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "nitrous_oxide_GM",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->RN2O)))
    {
        NCERR(retval);
    }

    /*Read in the carbon monoxide mole fraction in air (mol/mol) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "carbon_monoxide_GM",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->RCO)))
    {
        NCERR(retval);
    }

    /*Read in the methane mole fraction in air (mol/mol) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "methane_GM",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->RCH4)))
    {
        NCERR(retval);
    }

    /*Read in the oxygen mole fraction in air (mol/mol) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "oxygen_GM",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->RO2)))
    {
        NCERR(retval);
    }

    /*Read in the layer temperature (K) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "temp_layer",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->TEMP)))
    {
        NCERR(retval);
    }

    /*Read in all the pressure level (Pa) values.*/
    float *plevs = (float *)malloc((in->nlat)*(in->nphalf)*sizeof(float));
    if (plevs == NULL)
    {
        fprintf(stderr,
                "Error(radiationInputFieldsMalloc): malloc failed for the"
                " plevs array.\n");
        exit(EXIT_FAILURE);
    }
    if ((retval = nc_inq_varid(ncid,
                               "pres_level",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   plevs)))
    {
        NCERR(retval);
    }

    /*Calculate the average layer pressures (pa) and the number density
      [1/(m*m)] integrated across the layer, assuming that the layer is
      in hydrostatic equilibrium.*/
    for (i=0;i<in->nlat;i++)
    {
        for (j=0;j<in->npfull;j++)
        {
            zoffset = i*(in->npfull) + j;
            loffset = i*(in->nphalf) + j;
            in->PRESSM[zoffset] = 0.5*(plevs[loffset] + plevs[loffset+1]);
            dp = plevs[loffset+1] - plevs[loffset];
            if (dp < 0)
            {
                dp *= -1.;
            }
            in->N[zoffset] = dp/(Mkg*g);
        }
    }

    /*Next calculate the layer thicknesses (m).  Assume a hydrostatic
      atmosphere (dp/dz = m*g/V) and the ideal gas law (pV = n*R*T).  Assuming
      that the temperature (T) and mass (m) are constant in a layer, the
      layer thickness (dz) is approximately:
      dz = (dp*R*T)/(p*M*g),
      where R is the gas constant, g is the acceleration due to gravity, and
      M is the molar mass.  Taking dp to be the difference between the
      pressures at the top and bottom of the level and p to pressure at the
      middle of the layer, we can calculate dz.*/
    for (i=0;i<in->ntime;i++)
    {
        for (j=0;j<in->nlat;j++)
        {
            for (k=0;k<in->npfull;k++)
            {
                zoffset = i*(in->nlat)*(in->npfull) + j*(in->npfull) + k;
                loffset = j*(in->nphalf) + k;
                dp = log(plevs[loffset+1]) - log(plevs[loffset]);
                if (dp < 0)
                {
                    dp = -1.0*dp;
                }
                (in->DELTAZ)[zoffset] = (dp*R*((in->TEMP)[zoffset]))/(M*g);
            }
        }
    }
    free(plevs);
    plevs = NULL;

    /*Read in the surface temperature (K) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "surface_temperature",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->TSURF)))
    {
        NCERR(retval);
    }

    /*Read in the level temperatures.*/
    if ((retval = nc_inq_varid(ncid,
                               "temp_level",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->TLEV)))
    {
        NCERR(retval);
    }

    /*Read in the surface emissivity.*/
    if ((retval = nc_inq_varid(ncid,
                               "surface_emissivity",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->EMIS)))
    {
        NCERR(retval);
    }

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Convert read in the atmospheric fields from a "rfmip" formatted NetCDF file
  to the correct units required by grtcode.*/
int setOutputFieldsFromRfmip(radiationInputFields_t *in,
                             radiationOutputFields_t *out)
{
    /*Local variables.*/
    size_t offset;
    size_t poffset;
    size_t h2o_offset;
    size_t co2_offset;
    size_t o3_offset;
    size_t n2o_offset;
    size_t co_offset;
    size_t ch4_offset;
    size_t o2_offset;
    size_t t;
    size_t i;
    size_t k;
    REAL_t xh2o;

    /*Copy metadata into the output fields from the input fields.*/ 
    out->nlat = in->nlat;
    out->nlon = in->nlon;
    out->npfull = in->npfull;
    out->nphalf = in->nphalf;
    out->ntime = in->ntime;

    /*Malloc space for the radiation output fields.*/
    out->N = NULL;
    out->P = NULL;
    out->T = NULL;
    out->DELTAZ = NULL;
    out->PS = NULL;
    out->TSURF = NULL;
    out->TLEV = NULL;
    out->EMIS = NULL;
    radiationOutputFieldsMalloc(out);

    /*Set the radiation output fields.  The input fields are generally
      stored as (time,lat,pressure) (with lon=1), while the output fields
      are stored as (time,lat,lon,pressure) or
      (time,lat,lon,NUM_MOLS,pressure).*/
    for (t=0;t<in->ntime;++t)
    {
        for (i=0;i<in->nlat;++i)
        {
            out->TSURF[t*(in->nlat)+i] = in->TSURF[t*(in->nlat)+i];
            for (k=0;k<in->npfull;++k)
            {
                /*Calculate the offsets.*/
                offset = t*(in->nlat)*(in->npfull) + i*(in->npfull) + k;

                poffset = i*(in->npfull) + k;

                h2o_offset = t*(in->nlat)*NUM_MOLS*(in->npfull) +
                                 i*NUM_MOLS*(in->npfull) +
                                 0*(in->npfull) + k;

                co2_offset = t*(in->nlat)*NUM_MOLS*(in->npfull) +
                                 i*NUM_MOLS*(in->npfull) +
                                 1*(in->npfull) + k;

                o3_offset = t*(in->nlat)*NUM_MOLS*(in->npfull) +
                                i*NUM_MOLS*(in->npfull) +
                                2*(in->npfull) + k;

                n2o_offset = t*(in->nlat)*NUM_MOLS*(in->npfull) +
                                 i*NUM_MOLS*(in->npfull) +
                                 3*(in->npfull) + k;

                co_offset = t*(in->nlat)*NUM_MOLS*(in->npfull) +
                                i*NUM_MOLS*(in->npfull) +
                                4*(in->npfull) + k;

                ch4_offset = t*(in->nlat)*NUM_MOLS*(in->npfull) +
                                 i*NUM_MOLS*(in->npfull) +
                                 5*(in->npfull) + k;

                o2_offset = t*(in->nlat)*NUM_MOLS*(in->npfull) +
                                i*NUM_MOLS*(in->npfull) +
                                6*(in->npfull) + k;

                /*Pack the data into the radiation output fields.  The
                  units for each of the fields are:

                  input fields:                  output fields:
                  -----------------              ---------------
                  pressure (PRESSM) = Pa         pressure (P)  = atm
                  temperature (TEMP) = K         temperature (T) = K
                  deltaz (DELTAZ) = m            deltaz (DELTAZ) = cm
                  molecules/Area (N) = m^-2      molecules/area (N) = cm^-2
                                                 partial pressure (PS) = atm
                */
                out->P[offset] = ((REAL_t)(in->PRESSM[poffset]))*PaToAtm;

                out->T[offset] = ((REAL_t)(in->TEMP[offset]));

                out->DELTAZ[offset] = ((REAL_t)(in->DELTAZ[offset]))*MToCm;

                xh2o = ((REAL_t)(in->RH2O[offset]));

                out->N[h2o_offset] = ((xh2o*((REAL_t)(in->N[poffset])))/
                                          (1.0 + xh2o))*1.e-4;

                out->N[co2_offset] = ((((REAL_t)(in->RCO2[t]))*
                                          ((REAL_t)(in->N[poffset])))/
                                          (1.0 + xh2o))*1.e-4*1.e-6;

                out->N[o3_offset] = ((((REAL_t)(in->QO3[offset]))*
                                         ((REAL_t)(in->N[poffset])))/
                                         (1.0 + xh2o))*1.e-4;

                out->N[n2o_offset] = ((((REAL_t)(in->RN2O[t]))*
                                          ((REAL_t)(in->N[poffset])))/
                                          (1.0 + xh2o))*1.e-4*1.e-9;

                out->N[co_offset] = ((((REAL_t)(in->RCO[t]))*
                                         ((REAL_t)(in->N[poffset])))/
                                         (1.0 + xh2o))*1.e-4;

                out->N[ch4_offset] = ((((REAL_t)(in->RCH4[t]))*
                                          ((REAL_t)(in->N[poffset])))/
                                          (1.0 + xh2o))*1.e-4*1.e-9;

                out->N[o2_offset] = ((((REAL_t)(in->RO2[t]))*
                                         ((REAL_t)(in->N[poffset])))/
                                         (1.0 + xh2o))*1.e-4;

                out->PS[h2o_offset] = ((xh2o*((REAL_t)(in->PRESSM[poffset])))/
                                          (1.0 + xh2o))*PaToAtm;

                out->PS[co2_offset] = ((((REAL_t)(in->RCO2[t]))*
                                          ((REAL_t)(in->PRESSM[poffset])))/
                                          (1.0 + xh2o))*PaToAtm*1.e-6;

                out->PS[o3_offset] = ((((REAL_t)(in->QO3[offset]))*
                                         ((REAL_t)(in->PRESSM[poffset])))/
                                         (1.0 + xh2o))*PaToAtm;

                out->PS[n2o_offset] = ((((REAL_t)(in->RN2O[t]))*
                                          ((REAL_t)(in->PRESSM[poffset])))/
                                          (1.0 + xh2o))*PaToAtm*1.e-9;

                out->PS[co_offset] = ((((REAL_t)(in->RCO[t]))*
                                         ((REAL_t)(in->PRESSM[poffset])))/
                                         (1.0 + xh2o))*PaToAtm;

                out->PS[ch4_offset] = ((((REAL_t)(in->RCH4[t]))*
                                          ((REAL_t)(in->PRESSM[poffset])))/
                                          (1.0 + xh2o))*PaToAtm*1.e-9;

                out->PS[o2_offset] = ((((REAL_t)(in->RO2[t]))*
                                         ((REAL_t)(in->PRESSM[poffset])))/
                                         (1.0 + xh2o))*PaToAtm;
            }
        }
    }

    for (t=0;t<in->ntime;++t)
    {
        for (i=0;i<in->nlat;++i)
        {
            for (k=0;k<in->nphalf;++k)
            {
                offset = t*(in->nlat)*(in->nphalf) + i*(in->nphalf) + k;
                out->TLEV[offset] = (REAL_t)(in->TLEV[offset]);
            }
        }
    }

    for (i=0;i<in->nlat;++i)
    {
        out->EMIS[i] = in->EMIS[i];
    }

    /*Free the radiation input fields.*/
    free(in->RH2O);
    free(in->RCO2);
    free(in->QO3);
    free(in->RN2O);
    free(in->RCO);
    free(in->RCH4);
    free(in->RO2);
    free(in->PRESSM);
    free(in->TEMP);
    free(in->DELTAZ);
    free(in->TSURF);
    free(in->TLEV);
    free(in->EMIS);
    in->RH2O = NULL;
    in->RCO2 = NULL;
    in->QO3 = NULL;
    in->RN2O = NULL;
    in->RCO = NULL;
    in->RCH4 = NULL;
    in->RO2 = NULL;
    in->PRESSM = NULL;
    in->TEMP = NULL;
    in->DELTAZ = NULL;
    in->TSURF = NULL;
    in->TLEV = NULL;
    in->EMIS = NULL;

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Read in the atmospheric fields from a "gfdl" formatted NetCDF file.*/
int readGfdlFieldsFromFile(char fname[],
                           radiationInputFields_t* in)
{
    /*Local variables*/
    int retval;
    int ncid;
    int varid;
    int dimid;

    /*Open the inputted file. */
    if ((retval = nc_open(fname,
                          NC_NOWRITE,
                          &ncid)))
    {
        NCERR(retval);
    }

    /*Get the number of latitude points from the NetCDF dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "lat",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &in->nlat)))
    {
        NCERR(retval);
    }

    /*Get the number of longitude points from the NetCDF dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "lon",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &in->nlon)))
    {
        NCERR(retval);
    }

    /*Get the pressure layer midpoint (pfull) bounds from the NetCDF
      dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "pfull",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &in->npfull)))
    {
        NCERR(retval);
    }

    /*Get the pressure layer interface (phalf) bounds from the NetCDF
      dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "phalf",
                                &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &in->nphalf)))
    {
        NCERR(retval);
    }

    /*Get the time bounds from the NetCDF dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "time",
                               &dimid )))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &in->ntime)))
    {
        NCERR(retval);
    }

    /*Make sure that the number of pressure layer interfaces is one greater
      than the number of pressure layers.*/
    assert(in->nphalf == in->npfull+1);

    /*Null initialize the radiation input fields.*/
    in->RH2O = NULL;
    in->QO3 = NULL;
    in->PRESSM = NULL;
    in->TEMP = NULL;
    in->DELTAZ = NULL;
    in->TSURF = NULL;

    in->RCO2 = NULL;
    in->RN2O = NULL;
    in->RCO = NULL;
    in->RCH4 = NULL;
    in->RO2 = NULL;
    in->TLEV = NULL;
    in->EMIS = NULL;

    /*Malloc space for the radiation input fields.*/
    size_t size1 = (in->ntime)*(in->npfull)*(in->nlat)*(in->nlon);
    size_t size2 = (in->ntime)*(in->nphalf)*(in->nlat)*(in->nlon);
    in->RH2O = (float *)malloc(size1*sizeof(float));
    in->QO3 = (float *)malloc(size1*sizeof(float));
    in->PRESSM = (float *)malloc(size1*sizeof(float));
    in->TEMP = (float *)malloc(size1*sizeof(float));
    in->DELTAZ = (float *)malloc(size2*sizeof(float));
    in->TSURF = (float *)malloc((in->ntime)*(in->nlat)*(in->nlon)*sizeof(float));

    /*Make sure that the mallocs succeeded.*/
    if (in->RH2O == NULL || in->QO3 == NULL || in->PRESSM == NULL ||
            in->TEMP == NULL || in->DELTAZ == NULL || in->TSURF == NULL)
    {
        fprintf(stderr,
                "Error(radiationInputFieldsMalloc): malloc failed for the"
                " radiation input fields.\n");
        exit(EXIT_FAILURE);
    }

    /*Read in the water vapor mixing ratio (kg/kg) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "rh2o",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->RH2O)))
    {
        NCERR(retval);
    }

    /*Read in the ozone mixing ratio (kg/kg) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "qo3",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->QO3)))
    {
        NCERR(retval);
    }

    /*Read in the model layer pressure (Pa) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "pressm",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->PRESSM)))
    {
        NCERR(retval);
    }

    /*Read in the model layer temperature (K) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "temp",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->TEMP)))
    {
        NCERR(retval);
    }

    /*Read in the model layer thickness (m) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "deltaz",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->DELTAZ)))
    {
        NCERR(retval);
    }

    /*There does not seem to a surface temperature in GFDL style input
      files.  Use the lowest atmospheric level temperature?*/
    int i;
    int j;
    int k;
    for (i=0;i<(in->ntime);++i)
    {
        for (j=0;j<(in->nlat);++j)
        {
            for (k=0;k<(in->nlon);++k)
            {
                int toff = i*(in->npfull)*(in->nlat)*(in->nlon) +
                           (in->npfull-1)*(in->nlat)*(in->nlon) +
                           j*(in->nlon) + k;

                int soff = i*(in->nlat)*(in->nlon) +
                           j*(in->nlon) + k;

                in->TSURF[soff] = in->TEMP[toff];
            }
        }
    }

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Convert read in the atmospheric fields from a "rfmip" formatted NetCDF file
  to the correct units required by grtcode.*/
int setOutputFieldsFromGfdl(radiationInputFields_t *in,
                            radiationOutputFields_t *out)
{
    /*Local variables.*/
    size_t in_offset;
    size_t in_loffset;
    size_t out_offset;
    size_t h2o_offset;
    size_t o3_offset;
    size_t t;
    size_t k;
    size_t i;
    size_t j;
    REAL_t xh2o;

    /*Copy metadata into the output fields from the input fields.*/ 
    out->nlat = in->nlat;
    out->nlon = in->nlon;
    out->npfull = in->npfull;
    out->nphalf = in->nphalf;
    out->ntime  = in->ntime;

    /*Malloc space for the radiation output fields.*/
    out->N = NULL;
    out->P = NULL;
    out->T = NULL;
    out->DELTAZ = NULL;
    out->PS = NULL;
    out->TSURF = NULL;
    out->TLEV = NULL;
    out->EMIS = NULL;
    radiationOutputFieldsMalloc(out);

    /*Set the radiation output fields.  The input fields are generally
      stored as (time,pressure,lat,lon), while the output fields
      are stored as (time,lat,lon,pressure) or
      (time,lat,lon,NUM_MOLS,pressure).*/
    for (t=0;t<in->ntime;++t)
    {
        for (k=0;k<in->npfull;++k)
        {
            for (i=0;i<in->nlat;++i)
            {
                for (j=0;j<in->nlon;++j)
                {
                    /*Calculate the offsets.*/
                    in_offset = t*(in->npfull)*(in->nlat)*(in->nlon) +
                                    k*(in->nlat)*(in->nlon) + i*(in->nlon) + j;

                    in_loffset = t*(in->nphalf)*(in->nlat)*(in->nlon) +
                                     k*(in->nlat)*(in->nlon) + i*(in->nlon) + j;

                    out_offset = t*(out->nlat)*(out->nlon)*(out->npfull) +
                                     i*(out->nlon)*(out->npfull) +
                                     j*(out->npfull) + k;

                    h2o_offset = t*(out->nlat)*(out->nlon)*NUM_MOLS*(out->npfull) +
                                     i*(out->nlon)*NUM_MOLS*(out->npfull) +
                                     j*NUM_MOLS*(out->npfull) +
                                     0*(out->npfull) + k;

                    o3_offset = t*(out->nlat)*(out->nlon)*NUM_MOLS*(out->npfull) +
                                    i*(out->nlon)*NUM_MOLS*(out->npfull) +
                                    j*NUM_MOLS*(out->npfull) +
                                    2*(out->npfull) + k;

                    /*Pack the data into the radiation output fields.  The
                      units for each of the fields are:

                      input fields:                  output fields:
                      -----------------              ---------------
                      pressure (PRESSM) = Pa         pressure (P)  = atm
                      temperature (TEMP) = K         temperature (T) = K
                      deltaz (DELTAZ) = m            deltaz (DELTAZ) = cm
                                                     partial pressure (PS) = atm
                    */

                    out->P[out_offset] = ((REAL_t)(in->PRESSM[in_offset]))*PaToAtm;

                    out->T[out_offset] = (REAL_t)(in->TEMP[in_offset]);

                    out->DELTAZ[out_offset] = ((REAL_t)(in->DELTAZ[in_loffset]))*MToCm;

                    xh2o = ((REAL_t)(in->RH2O[in_offset]))*(DryAirM/WaterM);

                    out->PS[h2o_offset] = ((xh2o*((REAL_t)(in->PRESSM[in_offset])))/
                                              (1.0 + xh2o))*PaToAtm;

                    out->PS[o3_offset] = ((((REAL_t)(in->QO3[in_offset]))*
                                             ((REAL_t)(in->PRESSM[in_offset]))*
                                             (DryAirM/OzoneM))/
                                             (1.0 + xh2o))*PaToAtm;
                }
            }
        }
    }

    for (t=0;t<in->ntime;++t)
    {
        for (i=0;i<in->nlat;++i)
        {
            for (j=0;j<in->nlon;++j)
            {
                int off = t*(in->nlat)*(in->nlon) + i*(in->nlon) + j;
                out->TSURF[off] = in->TSURF[off];
            }
        }
    }

    printf("Level temperature values are not given in GFDL input files."
               "  Aborting.\n");
    exit(EXIT_FAILURE);

    /*Free the radiation input fields.*/
    free(in->RH2O);
    free(in->QO3);
    free(in->PRESSM);
    free(in->TEMP);
    free(in->DELTAZ);
    free(in->TSURF);
    in->RH2O = NULL;
    in->QO3 = NULL;
    in->PRESSM = NULL;
    in->TEMP = NULL;
    in->DELTAZ = NULL;
    in->TSURF = NULL;

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Malloc space for the radiation output fields.  This routine assumes that the
  time, latitude, longitude, and pressure layer dimensions have already been
  stored in the inputted structure.*/
int radiationOutputFieldsMalloc(radiationOutputFields_t* out)
{
    /*Set the size for values defined at pressure layer midpoints.*/
    size_t pfulllen = out->ntime*out->npfull*out->nlat*out->nlon;

    /*Make sure that the radiation output fields are null.*/
    if (out->N != NULL || out->P != NULL || out->T != NULL ||
        out->DELTAZ != NULL || out->PS != NULL || out->TSURF != NULL ||
        out->TLEV != NULL || out->EMIS != NULL)
    {
        fprintf(stderr,
                "Error(radiationOutputFieldsMalloc): please first initialize"
                " all radiation output fields to null.\n");
        exit(EXIT_FAILURE);
    }

    /*Malloc the radiation output fields.*/
    out->N = (REAL_t *)malloc(NUM_MOLS*pfulllen*sizeof(REAL_t));
    out->P = (REAL_t *)malloc(pfulllen*sizeof(REAL_t));
    out->T = (REAL_t *)malloc(pfulllen*sizeof(REAL_t));
    out->DELTAZ = (REAL_t *)malloc(pfulllen*sizeof(REAL_t));
    out->PS = (REAL_t *)malloc(NUM_MOLS*pfulllen*sizeof(REAL_t));
    out->TSURF = (REAL_t *)malloc(out->ntime*out->nlat*out->nlon*sizeof(REAL_t));
    out->TLEV = (REAL_t *)malloc(out->ntime*out->nlat*out->nlon*out->nphalf*sizeof(REAL_t));
    out->EMIS = (REAL_t *)malloc(out->nlat*out->nlon*sizeof(REAL_t));

    /*Make sure that the mallocs succeeded.*/
    if (out->N == NULL || out->P == NULL || out->T == NULL ||
        out->DELTAZ == NULL || out->PS == NULL || out->TSURF == NULL ||
        out->TLEV == NULL || out->EMIS == NULL)
    {
        fprintf(stderr,
                "Error(radiationOutputFieldsMalloc): malloc failed for the"
                " radiation output fields.\n");
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Free the radiation output fields.*/
int radiationOutputFieldsFree(radiationOutputFields_t* out)
{
    /*Make sure that the radiation output fields are not null.*/
    if (out->N == NULL || out->P == NULL || out->T == NULL ||
        out->DELTAZ == NULL || out->PS == NULL || out->TSURF == NULL ||
        out->TLEV == NULL || out->EMIS == NULL)
    {
        fprintf(stderr,
                "Error(radiationOutputFieldsFree): the output radiation"
                " fields are null.\n");
        exit(EXIT_FAILURE);
    }

    /*Free the radiation output fields.*/
    free(out->N);
    free(out->P);
    free(out->T);
    free(out->DELTAZ);
    free(out->PS);
    free(out->TSURF);
    free(out->TLEV);
    free(out->EMIS);

    /*Nullify the radiation output field pointers.*/
    out->N = NULL;
    out->P = NULL;
    out->T =  NULL;
    out->DELTAZ = NULL;
    out->PS = NULL;
    out->TSURF = NULL;
    out->TLEV = NULL;
    out->EMIS = NULL;

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Get the radiation input fields from the NetCDF field, and set the values
  for the radiation output fields.*/
int getAndSetAtmosFieldsFromFile(char fname[],
                                 char *f_format,
                                 radiationOutputFields_t* out)
{
    /*Local variables*/
    radiationInputFields_t in;

    /*Make sure that a valid file name was passed in.*/
    if (fname == NULL || strlen(fname) == 0)
    {
        fprintf(stderr,
                "Error(getAndSetAtmosFieldsFromFile): invalid inputted"
                " atmosphere file.\n");
        exit(EXIT_FAILURE);
    }

    /*Make sure that the file format is valid.*/
    if (f_format == NULL || strlen(f_format) == 0)
    {
        fprintf(stderr,
                "Error(%s,%s:%d): null or empty inputted atmosphere file"
                    " format.\n",
                __func__,
                __FILE__,
                __LINE__);
        exit(EXIT_FAILURE);
    }

    /*Read in the data and reshape it as necessary.*/
    if (strcmp("rfmip",f_format) == 0)
    {
        printf("Attempting to open and read input from rfmip formatted file"
                   " %s.\n",
                   fname);
        readRfmipFieldsFromFile(fname,
                                &in);
        printf("Attempting to reshape rfmip data and set output.\n");
        setOutputFieldsFromRfmip(&in,
                                 out);
    }
    else if (strcmp("gfdl",f_format) == 0)
    {
        printf("Attempting to open and read input from gfdl formatted file"
                   " %s.\n",
                   fname);
        readGfdlFieldsFromFile(fname,
                               &in);
        printf("Attempting to reshape gfdl data and set output.\n");
        setOutputFieldsFromGfdl(&in,
                                out);
    }
    else
    {
        fprintf(stderr,
                "Error(%s,%s:%d): invalid atmosphere file format %s.  This"
                    " must be one of either:\n\trfmip\n\tgfdl\n",
                __func__,
                __FILE__,
                __LINE__,
                f_format);
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
#ifndef SKIPMAIN
/*Unit test.*/
int test(char fname[],
         char *f_format)
{
    /*Local variables*/
    radiationOutputFields_t out;
    unsigned int i; /* = 0; */
    unsigned int j; /* = 0; */
    unsigned int k; /* = 47; */
    unsigned int t;
    unsigned int offset;

    getAndSetAtmosFieldsFromFile(fname,
                                 f_format,
                                 &out);

    printf("\nTest Results:\n\n");
    for (t=0;t<out.ntime;++t)
    {
        for (i=0;i<out.nlat;++i)
        {
            for (j=0;j<out.nlon;++j)
            {
                for (k=0;k<out.npfull;++k)
                {
                    offset = t*(out.nlat*out.nlon*out.npfull) +
                                 i*(out.nlon*out.npfull) + j*(out.npfull) +k;
                    printf("N(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n",
                           t,
                           i,
                           j,
                           k,
                           out.N[offset]);
                }
            }
        }
    }

    radiationOutputFieldsFree(&out);

    return EXIT_SUCCESS;
}

int test2(char fname[],
          char *f_format)
{
    /*Local variables*/
    radiationOutputFields_t out;
    unsigned int i = 0; 
    unsigned int j = 0;
    unsigned int k = 47;
    unsigned int t = 0;
    unsigned int offset;

    getAndSetAtmosFieldsFromFile(fname,
                                 f_format,
                                 &out);

    printf("\nTest Results:\n\n");
/*
    for (t=0;t<out.ntime;++t)
    {
        for (i=0;i<out.nlat;++i)
        {
            for (j=0;j<out.nlon;++j)
            {
                for (k=0;k<out.npfull;++k)
                {
*/
                    offset = t*(out.nlat*out.nlon*out.npfull) +
                                 i*(out.nlon*out.npfull) + j*(out.npfull) +k;
                    printf("DELTAZ(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n",
                           t,
                           i,
                           j,
                           k,
                           out.DELTAZ[offset]);
                    printf("P(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n",
                           t,
                           i,
                           j,
                           k,
                           out.P[offset]);
                    printf("T(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n",
                           t,
                           i,
                           j,
                           k,
                           out.T[offset]);
                    printf("N(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n",
                           t,
                           i,
                           j,
                           k,
                           out.N[offset]);
                    printf("PS(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n",
                           t,
                           i,
                           j,
                           k,
                           out.PS[offset]);
/*
                }
            }
        }
    }
*/

    radiationOutputFieldsFree(&out);

    return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        fprintf(stderr,
                "Usage: ./%s fname.nc format\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }

    if ( /* (test(argv[1],argv[2])) || */ (test2(argv[1],argv[2])) )
    {
        fprintf(stderr,
                "Tests FAILED!\n");
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
#endif

/*---------------------------------------------------------------------------*/
