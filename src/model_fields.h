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
    fp_t *SFC_DIR_ALB; /*Surface albedo for the direct beam (t,lon,lat).*/
    fp_t *SFC_DIF_ALB; /*Surface albedo for the diffuse beam (t,lon,lat).*/
    fp_t *COS_SOL_ZEN_ANG; /*Cosine of the Solar zenith angle
                             (t,lon,lat,lay).*/
    fp_t COS_DIF_BEAM_ANG; /*Cosine of the angle for the diffuse beam.*/
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
