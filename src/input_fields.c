#include <stdlib.h>
#include "constants.h"
#include "debug.h"
#include "floating_point_type.h"
#include "input_fields.h"
#include "model_fields.h"
#include "molecules.h"
#include "netcdf.h"
#include "utils.h"

static char const *time_name = "time";
static char const *lon_name = "longitude";
static char const *lat_name = "latitude";
static char const *layer_name = "layer";
static char const *level_name = "level";
static char const *pressure_name = "pressure";
static char const *temp_name = "temperature";
static char const *surf_temp_name = "surface temperature";
static char const *surf_emis_name = "surface emissivity";
static char *in_mol_names[NUM_MOL];

typedef struct dimension
{
    int dimid;
    size_t length;
} dim_t;

static int get_dim(dim_t *d,
                   char const *name,
                   int ncid)
{
    not_null(d);
    not_null(name);
    netcdf_check(nc_inq_dimid(ncid,
                              name,
                              &(d->dimid)));
    netcdf_check(nc_inq_dimlen(ncid,
                               d->dimid,
                               &(d->length)));
    if (d->length < 1)
    {
        fatal("dimension %s must have length (%zu) > 0.",
              name,
              d->length);
    }
    return SUCCESS;
}

static int get_var(double **buf,
                   char const *name,
                   int ncid,
                   dim_t *dims,
                   int num_dims)
{
    not_null(name);
    not_null(dims);
    int varid;
    netcdf_check(nc_inq_varid(ncid,
                              name,
                              &varid));

    /*Make sure that the variable has the correct number of dimensions.*/
    int v_num_dims;
    netcdf_check(nc_inq_varndims(ncid,
                                 varid,
                                 &v_num_dims));
    if (v_num_dims != num_dims)
    {
        fatal("variable %s is expected to depend on %d dimensions, but"
                  " actually depends on %d dimensions.",
              name,
              num_dims,
              v_num_dims);
    }

    /*Make sure that the variable has its dimensions in the correct order.*/
    int *v_dimids = NULL;
    check(malloc_ptr((void **)(&v_dimids),
                     sizeof(*v_dimids)*num_dims));
    netcdf_check(nc_inq_vardimid(ncid,
                                 varid,
                                 v_dimids));
    int i;
    int n = 1;
    for (i=0;i<num_dims;++i)
    {
        if (v_dimids[i] != dims[i].dimid)
        {
            fatal("dimension %d of variable %s is expected to be dimension"
                      " dimid=%d, but is actually dimension dimid=%d.",
                  i,
                  name,
                  dims[i].dimid,
                  v_dimids[i]);
        }
        n *= dims[i].length;
    }
    free(v_dimids);

    /*Make sure that the variable has the correct type.*/
    nc_type v_type;
    netcdf_check(nc_inq_vartype(ncid,
                                varid,
                                &v_type));
    if (v_type != NC_DOUBLE)
    {
        fatal("variable %s must contain doubles.",
              name);
    }

    /*Read in the data.*/
    double *b = NULL;
    check(malloc_ptr((void **)(&b),
                     sizeof(*b)*n));
    not_null(b);
    netcdf_check(nc_get_var_double(ncid,
                                   varid,
                                   b));
    *buf = b;
    return SUCCESS;
}

int get_input_data(req_model_fields_t * const out,
                   char const * const input_file,
                   int const * const mol_ids,
                   int const nMols,
                   double const * const molConc)
{
    not_null(out);
    not_null(input_file);
    not_null(mol_ids);
    not_null(molConc);

    /*Open the file.*/
    int ncid;
    netcdf_check(nc_open(input_file,
                         NC_NOWRITE,
                         &ncid));

    /*Get dimensions.*/
    dim_t time;
    check(get_dim(&time,
                  time_name,
                  ncid));
    dim_t lon;
    check(get_dim(&lon,
                  lon_name,
                  ncid));
    dim_t lat;
    check(get_dim(&lat,
                  lat_name,
                  ncid));
    dim_t layer;
    check(get_dim(&layer,
                  layer_name,
                  ncid));
    dim_t level;
    check(get_dim(&level,
                  level_name,
                  ncid));
    if (level.length != layer.length + 1)
    {
        fatal("The number of levels (%zu) must be equal to the number"
                  " of layers (%zu) + 1.",
              level.length,
              layer.length);
    }

    /*Set molecule names expected in the input file.*/
    in_mol_names[H2O] = "water vapor abundance";
    in_mol_names[CO2] = "carbon dioxide abundance";
    in_mol_names[O3] = "ozone abundance";
    in_mol_names[N2O] = "nitrous oxide abundance";
    in_mol_names[CO] = "carbon monoxide abundance";
    in_mol_names[CH4] = "methane abundance";
    in_mol_names[O2] = "oxygen abundance";

    /*Get variables.*/
    input_fields_t in;
    dim_t dims[4];
    dims[0] = time;
    dims[1] = lon;
    dims[2] = lat;
    check(get_var(&(in.TSURF),
                  surf_temp_name,
                  ncid,
                  dims,
                  3));
    check(get_var(&(in.EMIS),
                  surf_emis_name,
                  ncid,
                  dims,
                  3));
    dims[3] = level;
    check(get_var(&(in.P),
                  pressure_name,
                  ncid,
                  dims,
                  4));
    check(get_var(&(in.T),
                  temp_name,
                  ncid,
                  dims,
                  4));
    dims[3] = layer;
    in.x = NULL;
    check(malloc_ptr((void **)(&(in.x)),
                     sizeof(*(in.x))*nMols));
    not_null(in.x);
    int i;
    for (i=0;i<nMols;++i)
    {
        int id = mol_ids[i];
        check(get_var(&(in.x[i]),
                      in_mol_names[id],
                      ncid,
                      dims,
                      4));
    }
    netcdf_check(nc_close(ncid));

    /*Allocate arrays to hold the model input data.*/
    check(init_req_model_fields(out,
                                time.length,
                                lon.length,
                                lat.length,
                                level.length,
                                nMols));

    /*Convert the data from the file to that needed by the model.*/
    int n = time.length*lon.length*lat.length;
    int j;
    for (i=0;i<n;++i)
    {
        out->TSURF[i] = (fp_t)(in.TSURF[i]);
        out->EMIS[i] = (fp_t)(in.EMIS[i]);
        for (j=0;j<(int)level.length;++j)
        {
            int off = (int)(i*level.length + j);
            out->P[off] = (fp_t)(in.P[off]*PA_TO_ATM);
            out->T[off] = (fp_t)(in.T[off]);
        }
    }

    /*Handle molecular concentrations.*/
    for (i=0;i<nMols;++i)
    {
        if (molConc[i] == MISSING_CONC || molConc[i] == CONC_FROM_FILE)
        {
            for (j=0;j<n;++j)
            {
                check(linear_interp(&((in.x[i])[j*layer.length]),
                                    layer.length,
                                    &((out->x[i])[j*level.length]),
                                    level.length));
            }
        }
        else
        {
            for (j=0;j<n*((int)level.length);++j)
            {
                (out->x[i])[j] = molConc[i];
            }
        }

        /*Convert units.*/
        for (j=0;j<n*((int)level.length);++j)
        {
            (out->x[i])[j] *= FROM_PPMV;
        }
    }

    /*Free memory that is no longer needed.*/
    free(in.P);
    free(in.T);
    free(in.TSURF);
    free(in.EMIS);
    for (i=0;i<nMols;i++)
    {
        free(in.x[i]);
    }
    free(in.x);
    return SUCCESS;
}
