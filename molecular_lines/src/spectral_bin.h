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

/** @file*/
#ifndef SPECTRAL_BIN_H_
#define SPECTRAL_BIN_H_

#include <stdint.h>
#include "floating_point_type.h"


/*Number of interpolation points per bin.  Here we mimic RFM and use a
  quadratic interpolation scheme, thus 3 points are required.*/
#define NIP 3


/** @brief Spectral bins.*/
typedef struct SpectralBins
{
    int num_layers; /**< Number of layers.*/
    fp_t w0; /**< Lower bound [1/cm] of spectral grid.*/
    fp_t wres; /**< Resolution [1/cm] of spectral grid.*/
    uint64_t num_wpoints; /**< Size of the spectral grid.*/
    uint64_t n; /**< Number of spectral bins.*/
    fp_t width; /**< The width [1/cm] of each bin (except potentially
                     the last one).*/
    uint64_t isize; /**< Size of arrays that will be used to do the
                         interpolation.*/
    int ppb; /**< Spectral points in each bin.*/
    int do_interp; /**< Flag telling if interpolation is necessary (i.e., if
                        there are more than 3 spectral points in each bin).*/
    int last_ppb; /**< Spectral points in the last bin.*/
    int do_last_interp; /**< Flag telling if interpolation is necessary for
                             the last bin (i.e., if it contains more than 3
                             spectral points).*/
    fp_t *w; /**< Wavenumbers [1/cm] in each bin (n,NIP).*/
    fp_t *tau; /**< Optical depths in each bin (layer,n,NIP).*/
    uint64_t *l; /**< Index of left-most spectral point in each bin (n).*/
    uint64_t *r; /**< Index of right-most spectral point in each bin (n).*/
    int gpu_id; /**< Id of gpu where memory is allocated (or HOST_ONLY).*/
} SpectralBins_t;


/** @brief Given spectral grid properties, define the spectral bins.
    @return SUCCESS or an error code.*/
int create_spectral_bins(SpectralBins_t *bins, /**< Spectral bins.*/
                         int const num_layers, /**< Number of atmospheric
                                                    layers.*/
                         double const w0, /**< Lower bound [1/cm] of spectral
                                               grid.*/
                         uint64_t const n, /**< Number of spectral grid points.*/
                         double const wres, /**< Resolution [1/cm] of spectral
                                                 grid.*/
                         double const bin_width, /**< Desired width [1/cm] of
                                                      each spectral bin.*/
                         int const gpu_id /**< GPU id.*/
                        );


/** @brief Free memory stored in the spectral bins.
    @return SUCCESS or an error code.*/
int destroy_spectral_bins(SpectralBins_t *bins /**< Spectral bins.*/
                         );


#endif
