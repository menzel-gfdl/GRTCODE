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

#ifndef WATER_VAPOR_CONTINUUM_H_
#define WATER_VAPOR_CONTINUUM_H_

#include <stdint.h>
#include "floating_point_type.h"


enum water_vapor_coefs
{
    MTCKD25_F296 = 0,
    MTCKD25_S296,
    CKDF,
    CKDS,
    NUM_COEFS
};


typedef struct WaterVaporContinuumCoefs
{
    fp_t **coefs;
    uint64_t num_wpoints;
    int gpu_id;
} WaterVaporContinuumCoefs_t;


int get_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc,
                                    char const * const h2o_ctm_dir,
                                    uint64_t const num_wpoints,
                                    double const w0,
                                    double const res,
                                    int const gpu_id);


int free_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc);


#endif
