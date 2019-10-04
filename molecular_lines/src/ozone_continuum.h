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
