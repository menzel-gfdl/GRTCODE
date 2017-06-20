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
#include <netcdf.h>
#include "parseNetcdfRadiation.h"

static const int NUM_MOLS = 7;

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
    int ndims_in;
    int nvars_in;
    int ngatts_in;
    int unlimdimid_in;

    /*Open the inputted file. */
    if ((retval = nc_open(fname,
                          NC_NOWRITE,
                          &ncid)))
    {
        NCERR(retval);
    }

    /*Inquire about the file.  Get the number of dimensions, variables,
      and global attributes, and id of the unlimited dimension.*/
    if ((retval = nc_inq(ncid,
                         &ndims_in,
                         &nvars_in,
                         &ngatts_in,
                         &unlimdimid_in)))
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

    /*Malloc space for the radiation input fields.*/
    in->RH2O = (float *)malloc((in->ntime)*(in->nlat)*(in->npfull)*sizeof(float));
    in->RCO2 = (float *)malloc((in->ntime)*sizeof(float));
    in->QO3 = (float *)malloc((in->ntime)*(in->nlat)*(in->npfull)*sizeof(float));
    in->RN2O = (float *)malloc((in->ntime)*sizeof(float));
    in->RCO = (float *)malloc((in->ntime)*sizeof(float));
    in->RCH4 = (float *)malloc((in->ntime)*sizeof(float));
    in->RO2 = (float *)malloc((in->ntime)*sizeof(float));
    in->PRESSM = (float *)malloc((in->nlat)*(in->npfull)*sizeof(float));
    in->TEMP = (float *)malloc((in->ntime)*(in->nlat)*(in->npfull)*sizeof(float));
    in->DELTAZ = (float *)malloc((in->ntime)*(in->nlat)*(in->npfull)*sizeof(float));

    /*Make sure that the mallocs succeeded.*/
    if (in->RH2O == NULL || in->RCO2 == NULL || in->QO3 == NULL ||
        in->RN2O == NULL || in->RCO == NULL || in->RCH4 == NULL ||
        in->RO2 == NULL || in->PRESSM == NULL || in->TEMP == NULL ||
        in->DELTAZ == NULL)
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

    /*Read in the layer pressure (pa) values.*/
    if ((retval = nc_inq_varid(ncid,
                               "pres_layer",
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

    /*Calculate the layer thicknesses (m).*/

    /*First read in all the pressure level (Pa) values.*/
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

    /*Next calculate the layer thicknesses.  Assume a hydrostatic atmosphere
      (dp/dz = m*g/V) and the ideal gas law (pV = n*R*T).  Assuming that
      the temperature (T) and mass (m) are constant in a layer, the
      layer thickness (dz) is approximately:
      dz = (dp*R*T)/(p*M*g),
      where R is the gas constant, g is the acceleration due to gravity, and
      M is the molar mass.  Taking dp to be the difference between the
      pressures at the top and bottom of the level and p to pressure at the
      middle of the layer, we can calculate dz.*/
    size_t i;
    size_t j;
    size_t k;
    size_t offset;
    size_t poffset;
    size_t loffset;
    float dp;
    float g = 9.80665; /*(m/s^2)*/
    float R = 8.3144598; /*(m^3*Pa)/(K*mol)*/
    float M = 0.02897; /*kg/mol*/
    for (i=0;i<in->ntime;i++)
    {
        for (j=0;j<in->nlat;j++)
        {
            for (k=0;k<in->npfull;k++)
            {
                zoffset = i*(in->nlat)*(in->npfull) + j*(in->npull) + k;
                poffset = j*(in->npfull) + k;
                loffset = j*(in->nphalf) + k;
                dp = (plevs[loffset+1] - plevs[loffset])
                if (dp < 0)
                {
                    dp = -1.0*dp;
                }
                (in->DELTAZ)[zoffset] = (dp*R*((in->TEMP)[zoffset]))/
                                            (M*g*((in->PRESSM)[poffset]));
            }
        }
    }
    free(plevs);
    plevs = NULL;

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

    /*Copy metadata into the output fields from the input fields.*/ 
    out->nlat = in->nlat;
    out->nlon = in->nlon;
    out->npfull = in->npfull;
    out->ntime = in->ntime;

    /*Malloc space for the radiation output fields.*/
    out->N = NULL;
    out->P = NULL;
    out->T = NULL;
    out->DELTAZ = NULL;
    out->PS = NULL;
    radiationOutputFieldsMalloc(out);

    /*Set the radiation output fields.  The input fields are generally
      stored as (time,lat,pressure) (with lon=1), while the output fields
      are stored as (time,lat,lon,pressure)*/
    for (t=0;t<in->ntime;++t)
    {
        for (i=0;i<in->nlat;++i)
        {
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
                                                 partial pressure (PS) = atm
                */
                out->P[offset] = ((REAL_t)(in->PRESSM[poffset]))*9.86923E-6;
                out->T[offset] = ((REAL_t)(in->TEMP[offset]));
                out->DELTAZ[offset] = (in->DELTAZ[offset])*100.;

                out->PS[h2o_offset] = getPartialPres(((REAL_t)(in->RH2O[offset])),
                                                     ((REAL_t)(in->PRESSM[poffset])),
                                                     18.)*9.86923E-6;

                out->PS[co2_offset] = getPartialPres(((REAL_t)(in->RCO2[t])),
                                                     ((REAL_t)(in->PRESSM[poffset])),
                                                     44.)*9.86923E-6;

                out->PS[o3_offset] = getPartialPres(((REAL_t)(in->QO3[offset])),
                                                    ((REAL_t)(in->PRESSM[poffset])),
                                                    48.)*9.86923E-6;

                out->PS[n2o_offset] = getPartialPres(((REAL_t)(in->RN2O[t])),
                                                     ((REAL_t)(in->PRESSM[poffset])),
                                                     44.)*9.86923E-6;

                out->PS[co_offset] = getPartialPres(((REAL_t)(in->RCO[t])),
                                                    ((REAL_t)(in->PRESSM[poffset])),
                                                    28.)*9.86923E-6;

                out->PS[ch4_offset] = getPartialPres(((REAL_t)(in->RCH4[t])),
                                                     ((REAL_t)(in->PRESSM[poffset])),
                                                     16.)*9.86923E-6;

                out->PS[o2_offset] = getPartialPres(((REAL_t)(in->RO2[t])),
                                                    ((REAL_t)(in->PRESSM[poffset])),
                                                    16.)*9.86923E-6;
            }
        }
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

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/



























/*---------------------------------------------------------------------------*/
/*Read in the atmospheric fields from a "gfdl" formatted NetCDF file.*/
int readInputFieldsFromFile(char fname[],
                            radiationInputFields_t* in)
{
    /*Local variables*/
    int retval;        /*NetCDF error code.*/
    int ncid;          /*NetCDF generated dataset id.*/
    int varid;         /*NetCDF variable id.*/
    int dimid;         /*NetCDF dimension id.*/
    int ndims_in;      /*Number of dimensions defined for the NetCDF dataset.*/
    int nvars_in;      /*Number of variables defined for the NetCDF dataset.*/
    int ngatts_in;     /*Number of global attributes defined for the NetCDF dataset.*/
    int unlimdimid_in; /*Id of the unlimited dimension for the NetCDF dataset.*/

    /*Open the inputted file. */
    if ((retval = nc_open(fname,
                          NC_NOWRITE,
                          &ncid)))
    {
        NCERR(retval);
    }

    /*Inquire about the file.  Get the number of dimensions, variables,
      and global attributes, and id of the unlimited dimension.*/
    if ((retval = nc_inq(ncid,
                         &ndims_in,
                         &nvars_in,
                         &ngatts_in,
                         &unlimdimid_in)))
    {
        NCERR(retval);
    }

    /*If we needed to test things regarding the file,
     they'd go here.*/

    /*Get the latitude coordinate bounds from the NetCDF dimension.*/
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

    /*Get the longitude coordinate bounds from the NetCDF dimension.*/
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
    in->QO3  = NULL;
    in->DPFLUX = NULL;
    in->PRESSM = NULL;
    in->TEMP   = NULL;
    in->DELTAZ = NULL;

    /*Malloc space for the radiation input fields.*/
    radiationInputFieldsMalloc(in);

    /*Get the rh2o variable id and read in the data.  This quantity is the
      water vapor mixing ratio (kg/kg).*/
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

    /*Get the qo3 variable id and read in the data.  This quantity is the
      ozone mixing ratio (kg/kg).*/
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

    /*Get the DPFLUX variable id and read in the data. This variable
      represents (dP/dz)*DELTAZ (hPa).*/
    if ((retval = nc_inq_varid(ncid,
                               "dpflux",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   in->DPFLUX)))
    {
        NCERR(retval);
    }

    /*Get the PRESSM variable id and read in the data. This quantity is the
      layer pressure (Pa).*/
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

    /*Get the TEMP variable id and read in the data. This quantity is the
      temperature (K).*/
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

    /*Get the DELTAZ variable id and read in the data.  This quantity is the
      layer thickness (m)*/
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

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Malloc space for the radiation input fields.  This routine assumes that the
  time latitude, longitude, and pressure layer dimensions have already been
  read in from the NetCDF file and stored in the inputted structure.*/
int radiationInputFieldsMalloc(radiationInputFields_t* in)
{
    /*Set the size for values defined at pressure layer midpoints.*/
    size_t pfulllen = in->ntime*in->npfull*in->nlat*in->nlon;

    /*Set the size for values defined at pressure layer interfaces.*/
    size_t phalflen = in->ntime*in->nphalf*in->nlat*in->nlon;

    /*Make sure that the radiation input fields are null.*/
    if (in->RH2O != NULL || in->QO3 != NULL || in->DPFLUX != NULL ||
        in->PRESSM != NULL || in->TEMP != NULL || in->DELTAZ != NULL)
    {
        fprintf(stderr,
                "Error(radiationInputFieldsMalloc): please first initialize"
                " all radiation input fields to null.\n");
        exit(EXIT_FAILURE);
    }

    /*Malloc space for the radiation input fields.*/
    in->RH2O = (float*)malloc(pfulllen*sizeof(*(in->RH2O)));
    in->QO3  = (float*)malloc(pfulllen*sizeof(*(in->QO3)));
    in->DPFLUX = (float*)malloc(pfulllen*sizeof(*(in->DPFLUX)));
    in->PRESSM = (float*)malloc(pfulllen*sizeof(*(in->PRESSM)));
    in->TEMP   = (float*)malloc(pfulllen*sizeof(*(in->TEMP)));
    in->DELTAZ = (float*)malloc(phalflen*sizeof(*(in->DELTAZ)));

    /*Make sure that the mallocs succeeded.*/
    if (in->RH2O == NULL || in->QO3 == NULL || in->DPFLUX == NULL ||
        in->PRESSM == NULL || in->TEMP == NULL || in->DELTAZ == NULL)
    {
        fprintf(stderr,
                "Error(radiationInputFieldsMalloc): malloc failed for the"
                " radiation input fields.\n");
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Free the radiation input fields.*/
int radiationInputFieldsFree(radiationInputFields_t* in)
{
    /*Make sure that the radiation input fields are not null.*/
    if (in->RH2O == NULL || in->QO3 == NULL || in->DPFLUX == NULL ||
        in->PRESSM == NULL || in->TEMP == NULL || in->DELTAZ == NULL)
    {
        fprintf(stderr,
                "Error(radiationInputFieldsFree): the input radiation"
                " fields are null.\n");
        exit(EXIT_FAILURE);
    }

    /*Free the radiation input fields.*/
    free(in->RH2O);
    free(in->QO3);
    free(in->DPFLUX);
    free(in->PRESSM);
    free(in->TEMP);
    free(in->DELTAZ);

    /*Nullify the radiation input field pointers.*/
    in->RH2O = NULL;
    in->QO3 = NULL;
    in->DPFLUX = NULL;
    in->PRESSM = NULL;
    in->TEMP = NULL;
    in->DELTAZ = NULL;

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Free the radiation output fields.*/
int radiationOutputFieldsFree(radiationOutputFields_t* out)
{
    /*Make sure that the radiation output fields are not null.*/
    if (out->N == NULL || out->P == NULL || out->T == NULL ||
        out->DELTAZ == NULL || out->PS == NULL)
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

    /*Nullify the radiation output field pointers.*/
    out->N = NULL;
    out->P = NULL;
    out->T =  NULL;
    out->DELTAZ = NULL;
    out->PS = NULL;

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Malloc space for the radiation output fields.  This routine assumes that the
  time latitude, longitude, and pressure layer dimensions have already been
  stored in the inputted structure.*/
int radiationOutputFieldsMalloc(radiationOutputFields_t* out)
{
    /*Set the size for values defined at pressure layer midpoints.*/
    size_t pfulllen = out->ntime*out->npfull*out->nlat*out->nlon;

    /*Make sure that the radiation output fields are null.*/
    if (out->N != NULL || out->P != NULL || out->T != NULL ||
        out->DELTAZ != NULL || out->PS != NULL)
    {
        fprintf(stderr,
                "Error(radiationOutputFieldsMalloc): please first initialize"
                " all radiation output fields to null.\n");
        exit(EXIT_FAILURE);
    }

    /*Malloc the radiation output fields.*/
    out->N = (REAL_t*)malloc(NUM_MOLS*pfulllen*sizeof(*(out->N)));
    out->P = (REAL_t*)malloc(pfulllen*sizeof(*(out->P)));
    out->T = (REAL_t*)malloc(pfulllen*sizeof(*(out->T)));
    out->DELTAZ = (REAL_t*)malloc(pfulllen*sizeof(*(out->DELTAZ)));
    out->PS = (REAL_t*)malloc(NUM_MOLS*pfulllen*sizeof(*(out->PS)));

    /*Make sure that the mallocs succeeded.*/
    if (out->N == NULL || out->P == NULL || out->T == NULL ||
        out->DELTAZ == NULL || out->PS == NULL)
    {
        fprintf(stderr,
                "Error(radiationOutputFieldsMalloc): malloc failed for the"
                " radiation output fields.\n");
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Calculate the partial pressure (Pa) for a molecular species.  Ratio is a
  dimensionless molecular mixing ratio (such as kg/kg), pressm is a layer
  pressure (Pa), and molarMass is the molar mass of a molecular species
  (g/mol).*/
REAL_t getPartialPres(const REAL_t ratio,
                      const REAL_t pressm,
                      const REAL_t molarMass)
{
    /*Local variables*/
    const REAL_t molarMassDryAir= 29; /*Molar mass of dry air (g/mol).*/

    return ratio*(molarMassDryAir/molarMass)*pressm;
}

/*---------------------------------------------------------------------------*/
/*Set the values of the radiation output fields.*/
int setOutputFields(radiationInputFields_t *in,
                    radiationOutputFields_t *out)
{
    /*Local variables.*/
    size_t ifoffset;  /*Index offset for the radiation input field.*/
    size_t ofoffset;  /*Index offset for the radiation output field.*/
    size_t ihoffset;  /*Index offseft for the radiation input field.*/
    size_t h2ooffset; /*Index offset for the radiation output field.*/
    size_t o3offset;  /*Index offset for the radiation output field.*/
    size_t i;         /*Loop variable.*/
    size_t j;         /*Loop variable.*/
    size_t k;         /*Loop variable.*/
    size_t t;         /*Loop variable.*/

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
    radiationOutputFieldsMalloc(out);

    /*Set the radiation output fields.  The input fields are stored as
      (time,pressure,lat,lon), while the output fields are stored as
      (time,lat,lon,pressure)*/
    for (t=0;t<in->ntime;++t)
    {
        for (k=0;k<in->npfull;++k)
        {
            for (i=0;i<in->nlat;++i)
            {
                for (j=0;j<in->nlon;++j)
                {
                    /*Calculate the offsets.*/
                    ifoffset = t*(in->npfull*in->nlat*in->nlon) +
                               k*(in->nlat*in->nlon) + i*(in->nlon) + j;
                    ofoffset = t*(in->nlat*in->nlon*in->npfull) +
                               i*(in->nlon*in->npfull) + j*(in->npfull) + k;
                    ihoffset = t*(in->nphalf*in->nlat*in->nlon) +
                               k*(in->nlat*in->nlon) + i*(in->nlon) + j;
                    h2ooffset = t*(in->nlat*in->nlon*in->npfull*NUM_MOLS)
                                + i*(in->nlon*in->npfull*NUM_MOLS)
                                + j*(in->npfull*NUM_MOLS)
                                + 0*in->npfull + k;
                    o3offset = t*(in->nlat*in->nlon*in->npfull*NUM_MOLS)
                               + i*(in->nlon*in->npfull*NUM_MOLS)
                               + j*(in->npfull*NUM_MOLS)
                               + 2*in->npfull + k;

                    /*Pack the data into the radiation output fields.  The
                      units for each of the fields are:

                      input fields:                  output fields:
                      -----------------              ---------------
                      pressure (P)    = Pa           pressure (PRESSM)  = atm
                      temperature (T) = K            temperature (TEMP) = K
                      deltaz (DELTAZ) = m            deltaz (DELTAZ)    = cm
                      partial pressure (PS) = atm
                    */
                    out->P[ofoffset] = ((REAL_t)in->PRESSM[ifoffset])*9.86923E-6;
                    out->T[ofoffset] = ((REAL_t)in->TEMP[ifoffset]);
                    /* if: method is already packing level path integral (ie, dpflux, cough) then: */
                    /* out->DELTAZ[ofoffset] = ((REAL_t)1); */
                    /* else: */
                    out->DELTAZ[ofoffset] = in->DELTAZ[ihoffset]*100.;
                    out->PS[h2ooffset] = getPartialPres(((REAL_t)in->RH2O[ifoffset]),
                                                        ((REAL_t)in->PRESSM[ifoffset]),
                                                        18.)*9.86923E-6;
                    out->PS[o3offset] = getPartialPres(((REAL_t)in->QO3[ifoffset]),
                                                       ((REAL_t)in->PRESSM[ifoffset]),
                                                       48.)*9.86923E-6;
                }
            }
        }
    }

    /*Free the radiation input fields.*/
    radiationInputFieldsFree(in);

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
                __LINE__,);
        exit(EXIT_FAILURE);
    }
    if (f_format != "rfmip" && f_format != "gfdl")
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

    /*Read in the radiation input fields.*/
    printf("Attempting to open and input from %s.\n",
           fname);
    if (f_format == "rfmip")
    {
        readRfmipFieldsFromFile(fname,
                                &in);
    }
    else
    {
        readInputFieldsFromFile(fname,
                                &in);
    }

    /*Set the radiation output fields.*/
    printf("Attempting to reshape, compute, and set output.\n");
    if (f_format == "rfmip")
    {
        setOutputFieldsFromRfmip(&in,
                                 out);
    }
    else
    {
        setOutputFields(&in,
                        out);
    }

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Unit test.*/
int test(char fname[]){
  radiationOutputFields_t out;
  unsigned int i; /* = 0; */
  unsigned int j; /* = 0; */
  unsigned int k; /* = 47; */
  unsigned int t;

  unsigned int offset;

  getAndSetAtmosFieldsFromFile(fname, &out);

  printf("\nTest Results:\n\n");
  for(t=0; t<out.ntime; ++t){
    for(i=0; i<out.nlat; ++i){
      for(j=0; j<out.nlon; ++j){
        for(k=0; k<out.npfull; ++k){
          offset = t*(out.nlat*out.nlon*out.npfull) + i*(out.nlon*out.npfull) + j*(out.npfull) +k;
          printf("N(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.N[offset]);
        }
      }
    }
  }

  radiationOutputFieldsFree(&out);

  return EXIT_SUCCESS;

}

int test2(char fname[]){
  radiationOutputFields_t out;
  unsigned int i = 0; 
  unsigned int j = 0;
  unsigned int k = 47;
  /* unsigned int kn1; */
  unsigned int t = 0;

  unsigned int offset;

  getAndSetAtmosFieldsFromFile(fname, &out);

  printf("\nTest Results:\n\n");
  /* for(t=0; t<out.ntime; ++t){ */
    /* for(i=0; i<out.nlat; ++i){ */
      /* for(j=0; j<out.nlon; ++j){ */
        /* for(k=0; k<out.npfull; ++k){ */
       /* for(k=out.npfull; k>0; --k){ */
       /*   kn1 = k-1; */
         offset = t*(out.nlat*out.nlon*out.npfull) + i*(out.nlon*out.npfull) + j*(out.npfull) +k;
         printf("DELTAZ(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.DELTAZ[offset]);
         printf("P(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.P[offset]);
         printf("T(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.T[offset]);
         printf("N(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.N[offset]);
         printf("PS(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.PS[offset]);
         
        /* } */
  /*     } */
  /*   } */
  /* } */

  radiationOutputFieldsFree(&out);

  return EXIT_SUCCESS;

}


#ifndef SKIPMAIN
int main(int argc, char* argv[]){
  if(argc != 2){
    fprintf(stderr,"./%s fname.nc\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  if ( /* (test(argv[1])) || */ (test2(argv[1])) ){  
    fprintf(stderr,"Tests FAILED!\n");
    exit(EXIT_FAILURE);
  }
    
  return EXIT_SUCCESS;

}
#endif

/*---------------------------------------------------------------------------*/
