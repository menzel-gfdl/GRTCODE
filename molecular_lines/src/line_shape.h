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

#ifndef LINE_SHAPE_H_
#define LINE_SHAPE_H_

#include "floating_point_type.h"


typedef struct LineShapeInputs
{
    fp_t w; /*Starting wavenumber where line shape is calculated [1/cm].*/
    uint64_t num_wpoints; /*Number of spectral points were the line shape
                            is calculated.*/
    fp_t wres; /*Spectral resolution [1/cm].*/
    fp_t line_center; /*Line center wavenumber [1/cm].*/
    fp_t lorentz_hwhm; /*Lorentz half-width at half-maximum [1/cm].*/
    fp_t doppler_hwhm; /*Doppler half-width at half-maximum [1/cm].*/
    fp_t eta; /*Mixing parameter for Ida voigt algorithm.*/
} LineShapeInputs_t;


#endif
