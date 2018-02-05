#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "netcdf.h"

#define FSL 64

enum req_dims
{
    T_ = 0,
    LAT_,
#ifdef REQUIRE_LON
    LON_,
#endif
    LAY_,
    LEV_,
    NUM_REQ_DIMS
};

struct dimension
{
    int id; /*Value in req_dims enum.*/
    int dimid;
    size_t length;
};
typedef struct dimension dimension_t;

int get_dimension(dimension_t **dimension,
                  char const * const name,
                  int const id,
                  int const ncid)
{
    not_null(dimension);
    is_null(*dimension);
    not_null(name);
    if (id < T_ || id > LEV_)
    {
        fatal("input id (%d) must be >= %d or <= %d.",
              id,
              T_,
              LEV_);
    }

    dimension_t *d = (dimension_t *)malloc(sizeof(*d));
    not_null(d);
    d->id = id;
    netcdf_check(nc_inq_dimid(ncid,
                              name,
                              &(d->dimid)));
    netcdf_check(nc_inq_dimlen(ncid,
                               d->dimid,
                               &(d->length)));
    *dimension = d;

    return SUCCESS;
}










struct field
{
    char name[FSL];
    char units[FSL];
    int is_unitless;
    int num_dims;
    int *dimids;
    nc_type type;
    void *data;
};
typedef struct field field_t;

int get_field(field_t **field,
              char const * const name,
              int const ncid,
              int const is_unitless,
              dimension_t **required_dims)
{
    not_null(field);
    is_null(*field);
    not_null(name);
    not_null(required_dims);
    int i;
    for (i=0;i<NUM_REQ_DIMS;++i)
    {
        not_null(required_dims[i]);
    }

    int varid;
    netcdf_check(nc_inq_varid(ncid,
                              name,
                              &varid));

    field_t *f = (field_t *)malloc(sizeof(*f));
    not_null(f);

    snprintf(f->name,
             FSL,
             "%s",
             name);

    snprintf(f->units,
             FSL,
             "");
    int natts;
    netcdf_check(nc_inq_varnatts(ncid,
                                 varid,
                                 &natts));
    for (i=0;i<natts;++i)
    {
        char attname[NC_MAX_NAME+1];
        netcdf_check(nc_inq_attname(ncid,
                                    varid,
                                    i,
                                    attname));
        int match = strcmp("units",
                           attname);
        if (match == 0)
        {
            nc_type atttype;
            netcdf_check(nc_inq_atttype(ncid,
                                        varid,
                                        attname,
                                        &atttype));
            size_t attlen;
            netcdf_check(nc_inq_attlen(ncid,
                                       varid,
                                       attname,
                                       &attlen));
            switch (atttype)
            {
                case NC_CHAR:
                    if (attlen > FSL-1)
                    {
                        fatal("units attribute for variable %s is longer"
                                  " (%zu) than is allowed (%d).",
                              name,
                              attlen,
                              FSL-1);
                    }
                    netcdf_check(nc_get_att_text(ncid,
                                                 varid,
                                                 attname,
                                                 f->units));
                    f->units[attlen] = '\0';
                    break;
                default:
                    fatal("units attribute netcdf type (%d) for variable %s"
                              " is not currently supported.",
                          f->type,
                          name);
            }
            break;
        }
    }
    f->is_unitless = is_unitless;

    netcdf_check(nc_inq_varndims(ncid,
                                 varid,
                                 &(f->num_dims)));
    f->dimids = (int *)malloc(sizeof(*(f->dimids))*(f->num_dims));
    not_null(f->dimids);
    netcdf_check(nc_inq_vardimid(ncid,
                                 varid,
                                 f->dimids));
    size_t total_size = 1;
    for (i=0;i<(f->num_dims);++i)
    {
        int found = 0;
        int j;
        for (j=0;j<NUM_REQ_DIMS;++j)
        {
            if (f->dimids[i] == required_dims[j]->dimid)
            {
                found = 1;
                total_size *= required_dims[j]->length;
                break;
            }
        }
        if (!found)
        {
            fatal("dimension %d of variable %s is unsupported.",
                  i,
                  name);
        }
    }

    netcdf_check(nc_inq_vartype(ncid,
                                varid,
                                &(f->type)));
    switch (f->type)
    {
        case NC_FLOAT:
        {
            float *data = (float *)malloc(sizeof(*data)*total_size);
            netcdf_check(nc_get_var_float(ncid,
                                          varid,
                                          data));
            f->data = data;
            break;
        }
        default:
            fatal("unsupported netcdf type (%d) for variable %s.",
                  f->type,
                  name);
    }

    *field = f;

    return SUCCESS;
}

int set_model_field(field_t const * const in,
                    fp_t *out,
                    int const * const out_dimids)
{


    return SUCCESS;
}









int read_input_data(char const * const filename,
                    char const * const time_name,
                    char const * const lat_name,
#ifdef REQUIRE_LON
                    char const * const lon_name,
#endif
                    char const * const layer_name,
                    char const * const level_name,
                    char const * const h2o_name)
{
    not_null(filename);
    not_null(time_name);
    not_null(lat_name);
#ifdef REQUIRE_LON
    not_null(lon_name);
#endif
    not_null(layer_name);
    not_null(level_name);

    int ncid;
    netcdf_check(nc_open(filename,
                         NC_NOWRITE,
                         &ncid));

    /*Required dimensions are as follows:
        required_dims[T_] = time
        required_dims[LAT_] = lat
        required_dims[LON_] = lon
        required_dims[LAY_] = layer
        required_dims[LEV_] = level
      Fields that depend on other dimensions will throw errors.*/
    dimension_t **required_dims = malloc(sizeof(*required_dims)*NUM_REQ_DIMS);
    int i;
    for (i=0;i<NUM_REQ_DIMS;++i)
    {
        required_dims[i] = NULL;
    }
    check(get_dimension(&(required_dims[T_]),
                        time_name,
                        T_,
                        ncid));
    check(get_dimension(&(required_dims[LAT_]),
                        lat_name,
                        LAT_,
                        ncid));
#ifdef REQUIRE_LON
    check(get_dimension(&(required_dims[LON_]),
                        lon_name,
                        LON_,
                        ncid));
#endif
    check(get_dimension(&(required_dims[LAY_]),
                        layer_name,
                        LAY_,
                        ncid));
    check(get_dimension(&(required_dims[LEV_]),
                        level_name,
                        LEV_,
                        ncid));

    /*Fields required by the model.*/
    














    field_t *xh2o = NULL;
    check(get_field(&xh2o,
                    h2o_name,
                    ncid,
                    1,
                    required_dims));
/*
    check(reshape_data(xh2o,
                       &out));
*/

    for (i=0;i<NUM_REQ_DIMS;++i)
    {
        free(required_dims[i]);
    }

    netcdf_check(nc_close(ncid));

    return SUCCESS;
}











/*

    f->ordering = 0;
    for (i=0;i<(f->num_dims);++i)
    {

    }


















































    *field = f;

    return SUCCESS;
}

int convert_field(field_t const * const in,
                  void * const out)
{
    int in_type;
    int i;
    for (i=0;i<(in->num_dims);++i)
    {
        switch (in->dimids[i])
        {
            case t_dimid:
                break;
            case lat_dimid:
                break;
            case lon_dimid:
                break;
            case layer_dimid:
                break;
            case level_dimid:
                break;
            default:
        }
    }


    switch (type)
    {
        case T_LAT_LON:
        case T_LAT_LON_LAYER:
        case T_LAT_LON_MOL_LAYER:
        case T_LAT_LON_LEVEL:
        default:
    }

    return SUCCESS;
}


























*/


int main()
{
    check(read_input_data("/home/Raymond.Menzel/grtcodev2/run/INPUT/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-0-3.0_none.nc",
                          "expt",
                          "profile",
                          "layer",
                          "level",
                          "water_vapor"));


    return SUCCESS;
}
