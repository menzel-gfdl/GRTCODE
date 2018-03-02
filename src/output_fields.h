#ifndef OUTPUT_FIELDS_H_
#define OUTPUT_FIELDS_H_

#include "floating_point_type.h"


typedef struct OutputFields
{
    fp_t *tau; /**<Optical depth (layer,wavenumber).*/
    fp_t *lw_flux_down; /**<Longwave downward fluxes (level) [J/(m^2 s)].*/
    fp_t *lw_flux_up; /**<Longwave upward fluxes (level) [J/(m^2 s)].*/
} OutputFields_t;


int alloc_output_fields(OutputFields_t * const var,
                        int const nws,
                        int const nlevels);


int free_output_fields(OutputFields_t * const var);


#endif
