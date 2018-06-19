#ifndef SOLAR_FLUX_H_
#define SOLAR_FLUX_H_

#include "floating_point_type.h"


typedef struct SolarFlux
{
    fp_t *incident_sw_flux; /*Incident solar flux [W/m] per wavenumber.*/
    int nws; /*Size of incident solar flux array.*/
    fp_t total_sw_flux; /*Solar flux [W/m^2] integrated over wavenumber.*/
} SolarFlux_t;


/*Read in the solar flux values.*/
int get_solar_flux(SolarFlux_t *sf,
                   unsigned int const nws,
                   int const w0,
                   double const res);


int free_solar_flux(SolarFlux_t *sf);


int put_solar_flux_on_device(SolarFlux_t const * const in,
                             SolarFlux_t * const out);


int remove_solar_flux_from_device(SolarFlux_t * const in);


#endif
