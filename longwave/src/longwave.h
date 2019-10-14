/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef LONGWAVE_H_
#define LONGWAVE_H_


#include <stdint.h>
#include "device.h"
#include "extern.h"
#include "floating_point_type.h"
#include "optics.h"
#include "spectral_grid.h"


/** @brief Longwave object.*/
typedef struct Longwave
{
    int num_levels; /**< Number of atmospheric levels.*/
    SpectralGrid_t grid; /**< Spectral grid.*/
    Device_t device; /**< Device.*/
    fp_t *layer_temperature; /**< Temperature [K] (layer).*/
    fp_t *level_temperature; /**< Temperature [K] (level).*/
    fp_t *emissivity; /**< Emissivity of the surface (wavenumber).*/
    fp_t *flux_up; /**< Upward radiative flux [W*cm/m^2] (level, wavenumber).*/
    fp_t *flux_down; /**< Downward radiative flux [W*cm/m^2] (level, wavenumber).*/
} Longwave_t;


/** @brief Reserve memory for the longwave.
    @return RS_SUCCESS or an error code.*/
EXTERN int create_longwave(Longwave_t * const lw, /**< Longwave object.*/
                           int const num_levels, /**< Number of atmospheric levels.*/
                           SpectralGrid_t const * const grid, /**< Spectral grid.*/
                           Device_t const * const device /**< Device.*/
                          );


/** @brief Free memory for the longwave.
    @return RS_SUCCESS or an error code.*/
EXTERN int destroy_longwave(Longwave_t * const lw /**< Longwave object.*/
                           );


/** @brief Calculate the longwave radiative fluxes at each spectral grid point at
           each atmospheric level in the column.
    @return RS_SUCCESS or an error code.*/
EXTERN int calculate_lw_fluxes(Longwave_t * const lw, /**< Longwave object.*/
                               Optics_t const * const optics, /**< Optics object.*/
                               fp_t const T_surf, /**< Surface temperature.*/
                               fp_t * const T_layers, /**< Temperature [K] (layer).*/
                               fp_t * const T_levels, /**< Temperature [K] (level).*/
                               fp_t * const emis, /**< Surface emissivity (wavenumber).*/
                               fp_t * const flux_up, /**< Upward flux [W*cm/m^2] (level, wavenumber).*/
                               fp_t * const flux_down, /**< Downward flux [W*cm/m^2] (level, wavenumber.*/
                               double const * const x, /**< Spectral lower bound [1/cm].*/
                               double const * const X /**< Spectral upper bound [1/cm].*/
                              );


/** @brief Get the size of the spectral grid.
    @return RS_SUCCESS or an error code.*/
EXTERN int lw_get_spectral_grid_size(Longwave_t const * const lw, /**< Longwave object.*/
                                     uint64_t * const n /**< Spectral grid size.*/
                                    );


#endif
