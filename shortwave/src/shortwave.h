#ifndef SHORTWAVE_H_
#define SHORTWAVE_H_

#include <stdint.h>
#include "device.h"
#include "extern.h"
#include "floating_point_type.h"
#include "optics.h"
#include "spectral_grid.h"


/** @brief Shortwave object.*/
typedef struct Shortwave
{
    int num_levels; /**< Number of atmospheric levels.*/
    SpectralGrid_t grid; /**< Spectral grid.*/
    Device_t device; /**< Device.*/
    fp_t *solar_flux; /**< Incident solar flux [W*cm/m^2].*/
    fp_t *flux_up; /**< Upward radiative flux [W*cm/m^2] (level, wavenumber).*/
    fp_t *flux_down; /**< Downward radiative flux [W*cm/m^2] (level, wavenumber).*/
} Shortwave_t;


/** @brief Reserve memory for the shortwave.
    @return RS_SUCCESS or an error code.*/
EXTERN int create_shortwave(Shortwave_t * const sw, /**< Shortwave object.*/
                            int const num_levels, /**< Number of atmospheric levels.*/
                            SpectralGrid_t const * const grid, /**< Spectral grid.*/
                            Device_t const * const device /**< Device.*/
                           );


/** @brief Free memory for the shortwave.
    @return RS_SUCCESS or an error code.*/
EXTERN int destroy_shortwave(Shortwave_t * const sw /**< Shortwave object.*/
                            );


/** @brief Calculate the shortwave radiative fluxes at each spectral grid point at
           each atmospheric level in the column.
    @return RS_SUCCESS or an error code.*/
EXTERN int calculate_sw_fluxes(Shortwave_t * const sw, /**< Shortwave object.*/
                               Optics_t const * const optics, /**< Optics object.*/
                               fp_t const mu_dir, /**< Cosine of zenith angle for direct beam.*/
                               fp_t const mu_dif, /**< Cosine of zenith angle for diffuse beam.*/
                               fp_t const sfc_alpha_dir, /**< Surface albedo for direct beam.*/
                               fp_t const sfc_alpha_dif, /**< Surface albedo for diffuse beam.*/
                               fp_t * const solar_flux, /**< Solar flux [W*cm/m^2] (wavenumber).*/
                               fp_t * const flux_up, /**< Upward flux [W*cm/m^2] (level, wavenumber).*/
                               fp_t * const flux_down /**< Downward flux [W*cm/m^2] (level, wavenumber).*/
                              );


/** @brief Get the size of the spectral grid.
    @return RS_SUCCESS or an error code.*/
EXTERN int sw_get_spectral_grid_size(Shortwave_t const * const sw, /**< Shortwave object.*/
                                     uint64_t * const n /**< Spectral grid size.*/
                                    );


#endif
