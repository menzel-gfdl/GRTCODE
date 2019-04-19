#ifndef SOLAR_FLUX_H_
#define SOLAR_FLUX_H_

#include "extern.h"
#include "solar_flux.h"
#include "spectral_grid.h"


/** @brief Solar flux object.*/
typedef struct SolarFlux
{
    SpectralGrid_t grid; /**< Spectral grid.*/
    fp_t *incident_flux; /**< Incident solar flux [cm] (wavenumber).*/
    uint64_t n; /**< Size of spectral grid.*/
} SolarFlux_t;


/** @brief Read in data for the solar flux.
    @return RS_SUCCESS or an error code.*/
EXTERN int create_solar_flux(SolarFlux_t * const solar_flux, /**< Solar flux object.*/
                             SpectralGrid_t const * const grid, /**< Spectral grid.*/
                             char const * const filepath /**< Solar flux csv file.*/
                            );


/** @brief Free memory for the solar flux.
    @return RS_SUCCESS or an error code.*/
EXTERN int destroy_solar_flux(SolarFlux_t * const solar_flux /**< Solar flux object.*/
                             );


#endif
