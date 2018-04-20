#include "debug.h"
#include "floating_point_type.h"
#include "netcdf.h"
#include "write_output.h"

static char * time_name = "time";
static char * lon_name = "longitude";
static char * lat_name = "latitude";
static char * layer_name = "layer";
static char * level_name = "level";
static char * w_name = "wavenumber";
static char * lw_flux_down_name = "lw_flux_down";
static char * lw_flux_up_name = "lw_flux_up";
static char * sw_flux_down_name = "sw_flux_down";
static char * sw_flux_up_name = "sw_flux_up";
static char * tau_name = "optical_depth";
static int t_dimid;
static int lon_dimid;
static int lat_dimid;
static int layer_dimid;
static int level_dimid;
static int lw_flux_down_varid;
static int lw_flux_up_varid;
static int sw_flux_down_varid;
static int sw_flux_up_varid;
static int w_dimid;
static int tau_varid;
static nc_type type;


int init_output_file(char const * const filename,
                     int * const ncid,
                     int const nlons,
                     int const nlats,
                     int const nlevels,
                     int const nws,
                     int const output_spectra)
{
    not_null(filename);
    not_null(ncid);

    /*Create the output file.*/
    netcdf_check(nc_create(filename,
                           NC_NETCDF4 | NC_CLOBBER,
                           ncid));
    netcdf_check(ncsetfill(*ncid,
                           NC_NOFILL));
    log_mesg("Opening output file %s (ncid=%d).",
             filename,
             *ncid);

    /*Define dimensions.*/
    netcdf_check(nc_def_dim(*ncid,
                            time_name,
                            NC_UNLIMITED,
                            &t_dimid));

    netcdf_check(nc_def_dim(*ncid,
                            lon_name,
                            nlons,
                            &lon_dimid));

    netcdf_check(nc_def_dim(*ncid,
                            lat_name,
                            nlats,
                            &lat_dimid));

    netcdf_check(nc_def_dim(*ncid,
                            layer_name,
                            nlevels-1,
                            &layer_dimid));

    netcdf_check(nc_def_dim(*ncid,
                            level_name,
                            nlevels,
                            &level_dimid));

    /*Define variables.*/
    if (sizeof(fp_t) == sizeof(float))
    {
        type = NC_FLOAT;
    }
    else if (sizeof(fp_t) == sizeof(double))
    {
        type = NC_DOUBLE;
    }
    else
    {
        fatal("size of floating point type fp_t (%zu) must be equal"
                 " to %zu (float) or %zu (double).",
              sizeof(fp_t),
              sizeof(float),
              sizeof(double));
    }
    int dimids[6] = {t_dimid,lon_dimid,lat_dimid,level_dimid};

    netcdf_check(nc_def_var(*ncid,
                            lw_flux_down_name,
                            type,
                            4,
                            dimids,
                            &lw_flux_down_varid));

    netcdf_check(nc_def_var(*ncid,
                            lw_flux_up_name,
                            type,
                            4,
                            dimids,
                            &lw_flux_up_varid));

    netcdf_check(nc_def_var(*ncid,
                            sw_flux_down_name,
                            type,
                            4,
                            dimids,
                            &sw_flux_down_varid));

    netcdf_check(nc_def_var(*ncid,
                            sw_flux_up_name,
                            type,
                            4,
                            dimids,
                            &sw_flux_up_varid));

    if (output_spectra)
    {
        netcdf_check(nc_def_dim(*ncid,
                                w_name,
                                nws,
                                &w_dimid));

        dimids[3] = layer_dimid;
        dimids[4] = w_dimid;

        netcdf_check(nc_def_var(*ncid,
                                tau_name,
                                type,
                                5,
                                dimids,
                                &tau_varid));
    }
    return SUCCESS;
}


int close_output_file(int const ncid)
{
    log_mesg("Closing output file (ncid=%d).",
             ncid);
    netcdf_check(nc_close(ncid));
    return SUCCESS;
}


int write_data_column(int const ncid,
                      fp_t *lw_flux_down,
                      fp_t *lw_flux_up,
                      fp_t *sw_flux_down,
                      fp_t *sw_flux_up,
                      fp_t *tau,
                      int const time,
                      int const lon,
                      int const lat,
                      int const nlevels,
                      int const nws,
                      int const output_spectra)
{
    not_null(lw_flux_down);
    not_null(lw_flux_up);
    not_null(sw_flux_down);
    not_null(sw_flux_up);
    size_t start[5] = {time,lon,lat,0,0};
    size_t count[5] = {1,1,1,nlevels,nws};
    if (type == NC_FLOAT)
    {
        netcdf_check(nc_put_vara_float(ncid,
                                       lw_flux_down_varid,
                                       start,
                                       count,
                                       (float *)lw_flux_down));
        netcdf_check(nc_put_vara_float(ncid,
                                       lw_flux_up_varid,
                                       start,
                                       count,
                                       (float *)lw_flux_up));
        netcdf_check(nc_put_vara_float(ncid,
                                       sw_flux_down_varid,
                                       start,
                                       count,
                                       (float *)sw_flux_down));
        netcdf_check(nc_put_vara_float(ncid,
                                       sw_flux_up_varid,
                                       start,
                                       count,
                                       (float *)sw_flux_up));
        if (output_spectra)
        {
            not_null(tau);
            count[3] = nlevels-1;
            netcdf_check(nc_put_vara_float(ncid,
                                           tau_varid,
                                           start,
                                           count,
                                           (float *)tau));
        }
    }
    else if (type == NC_DOUBLE)
    {
        netcdf_check(nc_put_vara_double(ncid,
                                        lw_flux_down_varid,
                                        start,
                                        count,
                                        (double *)lw_flux_down));
        netcdf_check(nc_put_vara_double(ncid,
                                        lw_flux_up_varid,
                                        start,
                                        count,
                                        (double *)lw_flux_up));
        netcdf_check(nc_put_vara_double(ncid,
                                        sw_flux_down_varid,
                                        start,
                                        count,
                                        (double *)sw_flux_down));
        netcdf_check(nc_put_vara_double(ncid,
                                        sw_flux_up_varid,
                                        start,
                                        count,
                                        (double *)sw_flux_up));
        if (output_spectra)
        {
            not_null(tau);
            count[3] = nlevels-1;
            netcdf_check(nc_put_vara_double(ncid,
                                            tau_varid,
                                            start,
                                            count,
                                            (double *)tau));
        }
    }
    return SUCCESS;
}
