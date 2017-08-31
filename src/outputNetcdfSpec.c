#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <netcdf.h>
#include "outputNetcdfSpec.h"

/*Handle errors by printing an error message and exiting with a non-zero
  status.*/
#define NCERR(e) \
{ \
    fprintf(stderr, \
            "Error: %s\n", \
            nc_strerror(e)); \
    exit(EXIT_FAILURE); \
}

void closeOpticalDepthOutput(int const ncid)
{
    int retval;

    /*Close the file. This frees up any internal netCDF resources
      associated with the file, and flushes any buffers.*/
    if ((retval = nc_close(ncid)))
    {
        NCERR(retval);
    }
}

void openOpticalDepthOutput(int * const ncid,
                            int * const varid,
                            char const FNAME[],
                            size_t const nlat,
                            size_t const nlon,
                            size_t const nlayers,
                            size_t const nF)
{
    int retval;
    int t_dimid;
    int lat_dimid;
    int lon_dimid;
    int lay_dimid;
    int lev_dimid;
    int f_dimid;
    const int ndims=5;
    int did[1];
    int dimids[ndims];

    /*Create the file.  The NC_CLOBBER parameter tells netCDF to overwrite
      this file, if it already exists. */
    if ((retval = nc_create(FNAME,NC_CLOBBER,ncid)))
    {
        NCERR(retval);
    }

    /*NOFILL does not prefill file, so avoids extraneous writing.*/
    ncsetfill(*ncid,NC_NOFILL); 

    /*Define the dimensions. NetCDF will hand back an ID for each.*/
    if ((retval = nc_def_dim(*ncid,"time",NC_UNLIMITED,&t_dimid)))
    {
        NCERR(retval);
    }

    if ((retval = nc_def_dim(*ncid,"lat",nlat,&lat_dimid)))
    {
        NCERR(retval);
    }

    if ((retval = nc_def_dim(*ncid,"lon",nlon,&lon_dimid)))
    {
        NCERR(retval);
    }

    if ((retval = nc_def_dim(*ncid,"pfull",nlayers,&lay_dimid)))
    {
        NCERR(retval);
    }

    if ((retval = nc_def_dim(*ncid,"phalf",nlayers+1,&lev_dimid)))
    {
        NCERR(retval);
    }

    if ((retval = nc_def_dim(*ncid,"wavenumber",nF,&f_dimid)))
    {
        NCERR(retval);
    }

    /*Define the variables.*/
    did[0] = t_dimid;
    if ((retval = nc_def_var(*ncid,"time",NC_FLOAT,1,did,&(varid[0]))))
    {
        NCERR(retval);
    }

    did[0] = lat_dimid;
    if ((retval = nc_def_var(*ncid,"lat",NC_FLOAT,1,did,&(varid[1]))))
    {
        NCERR(retval);
    }

    did[0] = lon_dimid;
    if ((retval = nc_def_var(*ncid,"lon",NC_FLOAT,1,did,&(varid[2]))))
    {
        NCERR(retval);
    }

    did[0] = lay_dimid;
    if ((retval = nc_def_var(*ncid,"pfull",NC_FLOAT,1,did,&(varid[3]))))
    {
        NCERR(retval);
    }

    did[0] = lev_dimid;
    if ((retval = nc_def_var(*ncid,"phalf",NC_FLOAT,1,did,&(varid[4]))))
    {
        NCERR(retval);
    }

    did[0] = f_dimid;
    if ((retval = nc_def_var(*ncid,"wavenumber",NC_FLOAT,1,did,&(varid[5]))))
    {
        NCERR(retval);
    }

    dimids[0] = t_dimid;
    dimids[1] = lat_dimid;
    dimids[2] = lon_dimid;
    dimids[3] = lay_dimid;
    dimids[4] = f_dimid;
    if ((retval = nc_def_var(*ncid,"OpticalDepth",NC_FLOAT,ndims,dimids,&(varid[6]))))
    {
        NCERR(retval);
    }

    dimids[3] = lev_dimid;
    if ((retval = nc_def_var(*ncid,"Fluxes",NC_FLOAT,ndims,dimids,&(varid[7]))))
    {
        NCERR(retval);
    }

    /*End define mode. This tells netCDF we are done defining metadata.*/
    if ((retval = nc_enddef(*ncid)))
    {
        NCERR(retval)
    }

    return;
}

void writeDimensionData(int const ncid,
                        int const varid,
                        size_t const dim_size,
                        float const * const dim_data)
{
    int retval;
    const int ndims = 1;
    size_t start[ndims];
    size_t count[ndims];

    count[0] = dim_size;
    start[0] = 0;

    if ((retval = nc_put_vara_float(ncid,varid,start,count,dim_data)))
    {
        NCERR(retval);
    }
}

void writeOpticalDepthOutputByColumn(int const ncid,
                                     int const * const varid,
                                     int const t,
                                     int const lat,
                                     int const lon,
                                     int const nlayers,
                                     int const nF,
                                     float const * const spectra,
                                     float const * const fluxes)
{
    int retval;
    const int ndims = 5;
    size_t start[ndims];
    size_t count[ndims];

    /*Write out optical depths.*/
    count[0] = 1;
    count[1] = 1;
    count[2] = 1;
    count[3] = nlayers;
    count[4] = nF;
    start[0] = t;
    start[1] = lat;
    start[2] = lon;
    start[3] = 0;
    start[4] = 0;
    if ((retval = nc_put_vara_float(ncid,varid[6],start,count,spectra)))
    {
        NCERR(retval);
    }

    /*Write out fluxes.*/
    count[3] += 1;
    if ((retval = nc_put_vara_float(ncid,varid[7],start,count,fluxes)))
    {
        NCERR(retval);
    }

    return;
}
