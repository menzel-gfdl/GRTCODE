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

#ifndef LAUNCH_H_
#define LAUNCH_H_

#include "floating_point_type.h"
#include "molecular_lines.h"


/** @brief Driver for optical depth calculation.
    @return RS_SUCCESS or an error code.*/
int launch(MolecularLines_t * const ml, /**< Molecular lines object.*/
           fp_t *p, /**< Pressure [atm] (level).*/
           fp_t *t, /**< Temperature [K] (level).*/
           fp_t * const tau /**< Optical depth (level, wavenumber).*/
          );


#endif
