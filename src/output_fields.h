#ifndef OUTPUT_FIELDS_H_
#define OUTPUT_FIELDS_H_

#include "floating_point_type.h"


typedef struct OutputFields
{
    fp_t *tau; /**<Optical depth (mechanisms,layer,wavenumber).*/
    fp_t *lw_flux_down; /**<Longwave downward fluxes (level) [J/(m^2 s)].*/
    fp_t *lw_flux_up; /**<Longwave upward fluxes (level) [J/(m^2 s)].*/
    fp_t *lw_flux_down_per_w; /**<Longwave downward fluxes per wavenumber
                                  (level,wavenumber) [(J cm)/(m^2 s)].*/
    fp_t *lw_flux_up_per_w; /**<Longwave upward fluxes per wavenumber
                                (level,wavenumber) [(J cm)/(m^2 s)].*/
    fp_t *sw_flux_down; /**<Shortwave downward fluxes (level) [J/(m^2 s)].*/
    fp_t *sw_flux_up; /**<Shortwave upward fluxes (level) [J/(m^2 s)].*/
    fp_t *sw_flux_down_per_w; /**<Shortwave downward fluxes per wavenumber
                                  (level,wavenumber) [(J cm)/(m^2 s)].*/
    fp_t *sw_flux_up_per_w; /**<Shortwave upward fluxes per wavenumber
                                (level,wavenumber) [(J cm)/(m^2 s)].*/
} OutputFields_t;


int alloc_output_fields(OutputFields_t * const var,
                        int const nws,
                        int const nlevels,
                        int const device_launch);


int free_output_fields(OutputFields_t * const var,
                       int const device_launch);


#endif
