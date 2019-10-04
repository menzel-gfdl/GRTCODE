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

#ifndef PARSE_HITRAN_FILE_H_
#define PARSE_HITRAN_FILE_H_

#include <stdint.h>
#include "floating_point_type.h"


typedef struct LineParams
{
    int *iso;
    fp_t *vnn;
    fp_t *snn;
    fp_t *yair;
    fp_t *yself;
    fp_t *en;
    fp_t *n;
    fp_t *d;
    uint64_t num_lines;
    int gpu_id;
} LineParams_t;


int free_line_params(LineParams_t * const line_params);


int parse_hitran_file(LineParams_t * const line_params,
                      char const * const filename,
                      int const mol_id,
                      double const w0,
                      double const wn,
                      int const gpu_id);


#endif
