#ifndef MODEL_FIELDS_H_
#define MODEL_FIELDS_H_

#include "floating_point_type.h"

typedef struct req_model_fields
{
    int ntime; /*Size of time dimension.*/
    int nlon; /*Size of latitdue dimension.*/
    int nlat; /*Size of longitude dimension.*/
    int nlevel; /*Size of levels dimension (height).*/
    int nmol; /*Number of molecules.*/
    fp_t *P; /*Pressure at layer interfaces (t,lon,lat,lev) [atm].*/
    fp_t *T; /*Temperature at layer interfaces (t,lon,lat,lev) [K].*/
    fp_t *TSURF; /*Surface temperature (t,lon,lat) [K].*/
    fp_t *EMIS; /*Surface emissivity (t,lon,lat).*/
    fp_t **x; /*molecular abundances at layer interfaces (t,lon,lat,lev).*/
} req_model_fields_t;

int init_req_model_fields(req_model_fields_t * const fields,
                          int const nt,
                          int const nlon,
                          int const nlat,
                          int const nlev,
                          int const nmol);

int free_req_model_fields(req_model_fields_t * const fields);

#endif
