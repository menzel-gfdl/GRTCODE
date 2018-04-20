#ifndef SOLAR_FLUX_H_
#define SOLAR_FLUX_H_

#include "floating_point_type.h"


typedef struct SolarFlux
{
    fp_t *incident_sw_flux; /*Incident solar flux [W/m] per wavenumber.*/
} SolarFlux_t;


/*Read in the solar flux values.*/
int get_solar_flux(char const * const filepath,
                   SolarFlux_t *sf,
                   unsigned int const nws,
                   int const w0,
                   double const res,
                   int put_on_device);


int free_solar_flux(SolarFlux_t *sf,
                    int const on_device);


#endif
