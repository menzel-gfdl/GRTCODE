#ifndef INPUT_FIELDS_H_
#define INPUT_FIELDS_H_

#include "model_fields.h"


/*Container for all atmospheric input data.*/
typedef struct input_fields
{
    double *P; /*Pressure [Pa] at layer interfaces (t,lon,lat,lev).*/
    double *T; /*Temperature [K] at layer interfaces (t,lon,lat,lev).*/
    double *TSURF; /*Surface temperature [K] (t,lon,lat).*/
    double *EMIS; /*Surface emissivity (t,lon,lat).*/
    double *SFC_DIR_ALB; /*Surface albedo for the direct beam (t,lon,lat).*/
    double *SFC_DIF_ALB; /*Surface albedo for the diffuse beam (t,lon,lat).*/
    double *SOL_ZEN_ANG; /*Solar zenith angle [degrees] (t,lon,lat).*/
    double *TOTAL_SOL_FLUX; /*Total solar flux [W/m^2] (t,lon,lat).*/
    double **x; /*Moleculare abundances [ppmv] in layers (t,lon,lat,lay).*/
} input_fields_t;


/*Read in the atmosphere input data, and convert it into the form required
  by the model.*/
int get_input_data(req_model_fields_t * const out,
                   char const * const input_file,
                   int const * const mol_ids,
                   int const nMols,
                   double const * const molConc);


#endif
