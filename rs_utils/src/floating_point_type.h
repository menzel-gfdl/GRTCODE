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

#ifndef FLOATING_POINT_TYPE_H_
#define FLOATING_POINT_TYPE_H_

#ifdef SINGLE_PRECISION
#define TYPE float
/*ln(max(float)) = 88.72284.  Let's use 80 so we have some runway.*/
#define MAX_EXP_ARG 80.f
#define EXP expf
#define POW powf
#define SQRT sqrtf
#define ABS fabsf
#else
#define TYPE double
/*ln(max(double)) = 709.782712893384.  Let's use 700 so we have some runway.*/
#define MAX_EXP_ARG 700.
#define EXP exp
#define POW pow
#define SQRT sqrt
#define ABS fabs
#endif


typedef TYPE fp_t;


#endif
