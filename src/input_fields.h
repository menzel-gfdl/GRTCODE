#ifndef INPUT_FIELDS_H_
#define INPUT_FIELDS_H_

#include "model_fields.h"

typedef struct input_fields
{
    double *P; /*Pressure at layer interfaces (t,lon,lat,lev) [Pa].*/
    double *T; /*Temperature at layer interfaces (t,lon,lat,lev) [K].*/
    double *TSURF; /*Surface temperature (t,lon,lat) [K].*/
    double *EMIS; /*Surface emissivity (t,lon,lat).*/
    double **x; /*Moleculare abundances in layers (t,lon,lat,lay) [ppmv].*/
} input_fields_t;

int get_input_data(req_model_fields_t * const out,
                   char const * const input_file,
                   int const * const mol_ids,
                   int const nMols,
                   double const * const molConc);

#endif
