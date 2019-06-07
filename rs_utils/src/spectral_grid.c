#include <stdint.h>
#include "debug.h"
#include "rs_config.h"
#include "spectral_grid.h"


/*Initialize a spectral grid.*/
int create_spectral_grid(SpectralGrid_t * const grid, double const w0, double const wn,
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
    char *mesg = "Spectral grid properties:\n\tlower bound: %e [1/cm]\n\t"
                 "upper bound: %e [1/cm]\n\tresolution: %e [1/cm]\n\t"
                 "total size: %zu grid points";
    log_info(mesg, w0, wn, dw, grid->n);
    return RS_SUCCESS;
}


/*Determine if two spectral grids are the same.*/
int compare_spectral_grids(SpectralGrid_t const * const one,
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
