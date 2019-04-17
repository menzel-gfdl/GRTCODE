#include <stdint.h>
#include "config.h"
#include "debug.h"
#include "spectral_grid.h"


/*Initialize a spectral grid.*/
int create_spectral_grid(SpectralGrid_t * const grid, double const w0, double const wn,
                         double const dw, int const * const gpu_id)
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

    int num_devices;
    catch(get_num_gpus(&num_devices, 1));
    if (gpu_id != NULL)
    {
        if (*gpu_id != HOST_ONLY)
        {
            in_range(*gpu_id, 0, num_devices);
        }
        grid->gpu_id = *gpu_id;
    }
    else if (num_devices > 0)
    {
        grid->gpu_id = DEFAULT_GPU;
    }
    else
    {
        grid->gpu_id = HOST_ONLY;
    }

    fp_t *w = NULL;
    gmalloc(w, grid->n, HOST_ONLY);
    uint64_t i;
    for (i=0; i<grid->n; ++i)
    {
        w[i] = w0 + i*dw;
    }
    if (grid->gpu_id == HOST_ONLY)
    {
        grid->w = w;
    }
    else
    {
        gmalloc(grid->w, grid->n, grid->gpu_id);
        gmemcpy(grid->w, w, grid->n, grid->gpu_id, FROM_HOST);
        gfree(w, HOST_ONLY);
    }
    return RS_SUCCESS;
}


/*Free memory stored in a spectral grid object.*/
int destroy_spectral_grid(SpectralGrid_t * const grid)
{
    not_null(grid);
    gfree(grid->w, grid->gpu_id);
    return RS_SUCCESS;
}


/*Determine if two spectral grids are the same.*/
int compare_spectral_grids(SpectralGrid_t const * const one,
                           SpectralGrid_t const * const two, int * const result)
{
    not_null(one);
    not_null(two);
    not_null(result);
    if ((one->w0 == two->w0) && (one->wn == two->wn) && (one->dw == two->dw) &&
        (one->gpu_id == two->gpu_id))
    {
        *result = 1;
    }
    else
    {
        *result = 0;
    }
    return RS_SUCCESS;
}
