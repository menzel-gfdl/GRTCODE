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

#include <stdint.h>
#include "debug.h"
#include "extern.h"
#include "rs_config.h"
#include "spectral_grid.h"


/*Initialize a spectral grid.*/
EXTERN int create_spectral_grid(SpectralGrid_t * const grid, double const w0, double const wn,
                                double const dw)
{
    not_null(grid);
    in_range(w0, MIN_WAVENUMBER, MAX_WAVENUMBER);
    grid->w0 = w0;
    in_range(wn, w0, MAX_WAVENUMBER);
    grid->wn = wn;
    in_range(dw, MIN_RESOLUTION, MAX_RESOLUTION);
    grid->dw = dw;
    grid->n = ceil((wn-w0)/dw) + 1.;
    char const *mesg = "Spectral grid properties:\n\tlower bound: %e [1/cm]\n\t"
                       "upper bound: %e [1/cm]\n\tresolution: %e [1/cm]\n\t"
                       "total size: %zu grid points";
    log_info(mesg, w0, wn, dw, grid->n);
    return RS_SUCCESS;
}


/*Determine if two spectral grids are the same.*/
EXTERN int compare_spectral_grids(SpectralGrid_t const * const one,
                                  SpectralGrid_t const * const two, int * const result)
{
    not_null(one);
    not_null(two);
    not_null(result);
    if ((one->w0 == two->w0) && (one->wn == two->wn) && (one->dw == two->dw))
    {
        *result = 1;
    }
    else
    {
        *result = 0;
    }
    return RS_SUCCESS;
}
