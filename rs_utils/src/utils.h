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

#ifndef UTILS_H_
#define UTILS_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include "extern.h"
#include "floating_point_type.h"


/** @brief Malloc memory, with error checks.
    @return RS_SUCCESS or an error code.*/
EXTERN int malloc_ptr(void ** const p,
                      size_t const num_bytes);


/** @brief Free malloced memory, with error checks.
    @return RS_SUCCESS or an error code.*/
EXTERN int free_ptr(void ** const p);


/** @brief Copy a string into a buffer, checking its length.
    @return RS_SUCCESS or an error code.*/
int copy_str(char * const dest,
             char const * const src,
             size_t const len);


/** @brief Open a file, with error checks.
    @return RS_SUCCESS or an error code.*/
int open_file(FILE **file,
              char const * const name,
              char const * const mode);


/** @brief Convert a string to an integer.
    @return RS_SUCCESS or an error code.*/
int to_int(char const * const s,
           int * const i);


/** @brief Convert a string to a double.
    @return RS_SUCCESS or an error code.*/
int to_double(char const * const s,
              double * const d);


/** @brief Convert a double to fp_t.
    @return RS_SUCCESS or an error code.*/
int to_fp_t(double const d,
            fp_t * const f);


/** @brief Find the array indices that bracket the input value.
    @return RS_SUCCESS or an error code.*/
int get_sorted_bounds(fp_t const val, /**< Value to bracket.*/
                      fp_t const * const array, /**< Array.*/
                      int const array_size, /**< Size of array.*/
                      int * const left, /**< Index left of input value.*/
                      int * const right /**< Index right of input value.*/
                     );


/** @brief Perform a linear interpolation.
    @return RS_SUCCESS or an error code.*/
int linear_interpolation(fp_t const * const x, /**< Array of x values.*/
                         fp_t const * const y, /**< Array of y values.*/
                         int const xy_size, /**< Size of arrays.*/
                         fp_t const val, /**< Input x value.*/
                         fp_t * const out /**< Interpolated y value.*/
                        );


/** @brief Perform an integral using a reimann sum.
    @return RS_SUCCESS or an error code.*/
int reimann_sum(fp_t const * const data, /**< Data to be integrated.*/
                int const data_size, /**< Size of input data array.*/
                fp_t const dx, /**< Spacing between data points.*/
                fp_t * const out /**< Integral result.*/
               );


/** @brief Turn on bit in bit field.
    @return RS_SUCCESS or an error code.*/
int activate(uint64_t * const bit_field, /**< Bit field.*/
             int const index /**< Bit index (0 - 63).*/
            );


/** @brief Check if bit is turned on in bit field.
    @return 0 if bit is off, else nonzero.*/
int is_active(uint64_t const bit_field, /**< Bit field.*/
              int const index /**< Bit index (0 - 63).*/
             );


#endif
