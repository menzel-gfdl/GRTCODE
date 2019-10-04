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

#ifndef RETURN_CODES_H_
#define RETURN_CODES_H_


enum return_codes
{
    RS_SUCCESS,
    RS_INVALID_ERR,
    RS_DIVBYZERO_ERR,
    RS_OVERFLOW_ERR,
    RS_UNDERFLOW_ERR,
    RS_SENTINEL_ERR,
    RS_NULL_ERR,
    RS_NON_NULL_ERR,
    RS_RANGE_ERR,
    RS_VALUE_ERR,
    RS_COMPILER_ERR,
    RS_IO_ERR,
    RS_GPU_ERR
};


#endif
