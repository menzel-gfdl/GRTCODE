#ifndef OZONE_CONTINUUM_H_
#define OZONE_CONTINUUM_H_

#include <stdint.h>
#include "floating_point_type.h"


typedef struct OzoneContinuumCoefs
{
    fp_t *cross_section; /*Ozone continuum cross-section [cm^2].*/
    uint64_t num_wpoints;
    int gpu_id;
} OzoneContinuumCoefs_t;


/*Read in the ozone continuum coefficients.*/
int get_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc,
                              char const * const o3_ctm_dir,
                              uint64_t const num_wpoints,
                              double const w0,
                              double const res,
                              int const gpu_id);


int free_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc);


#endif
