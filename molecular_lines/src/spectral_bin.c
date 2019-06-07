#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "floating_point_type.h"
#include "spectral_bin.h"
#include "utils.h"


/*Set parameters and allocate arrays.*/
int create_spectral_bins(SpectralBins_t *bins,
                         int const num_layers,
                         double const w0,
                         uint64_t const n,
                         double const wres,
                         double const bin_width,
                         int const gpu_id)
{
    not_null(bins);

    /*Store the necessary properties that are used to create these spectral
      bins.*/
    bins->num_layers = num_layers;
    bins->w0 = w0;
    bins->wres = wres;
    bins->num_wpoints = n;
    bins->width = bin_width;

    /*Determine the number of spectral points per bin.  Each bin will contain
      at least one spectral point.  Each spectral point may only exist in
      one bin.  Interpolation is only required if there are more than 3
      spectral points per bin.*/
    bins->ppb = floor(bins->width/wres) + 1;
    bins->do_interp = bins->ppb > 3 ? 1 : 0;

    /*The last bin might have a smaller number of spectral points than all
      the rest, if the numbers do not divide evenly.*/
    bins->last_ppb = n % bins->ppb;
    bins->last_ppb = bins->last_ppb == 0 ? bins->ppb : bins->last_ppb;
    bins->do_last_interp = bins->last_ppb > 3 ? 1 : 0;

    /*Allocate arrays to store wavenumber and optical depth values at each
      interpolation point in each bin, as well as the left-most and
      right-most spectral grid indices in each bin.*/
    bins->n = n/bins->ppb;
    if (bins->ppb != bins->last_ppb)
    {
        (bins->n)++;
    }
    bins->isize = NIP*bins->n;
    uint64_t l[bins->n];
    uint64_t r[bins->n];
    fp_t w[bins->isize];

    /*Interpolation wavenumbers defined as follows:
      - First spectral point in the bin
      - Last spectral point in the bin
      - Midpoint between the first and last spectral points (not
        necessarily on the spectral grid.)*/
    uint64_t i;
    for (i=0;i<bins->n;++i)
    {
        l[i] = i*bins->ppb;
        int s = i < (bins->n - 1) ? bins->ppb : bins->last_ppb;
        r[i] = l[i] + s - 1;
        uint64_t o = i*NIP;
        w[o] = w0 + bins->ppb*i*wres;
        w[o+(NIP-1)] = w[o] + (s-1)*wres;
        w[o+1] = 0.5f*(w[o] + w[o+(NIP-1)]);
    }

    gmalloc(bins->l,bins->n,gpu_id);
    gmemcpy(bins->l,l,bins->n,gpu_id,FROM_HOST);
    gmalloc(bins->r,bins->n,gpu_id);
    gmemcpy(bins->r,r,bins->n,gpu_id,FROM_HOST);
    gmalloc(bins->w,bins->isize,gpu_id);
    gmemcpy(bins->w,w,bins->isize,gpu_id,FROM_HOST);
    gmalloc(bins->tau,bins->isize*bins->num_layers,gpu_id);
    bins->gpu_id = gpu_id;
    return RS_SUCCESS;
}


/*Free arrays.*/
int destroy_spectral_bins(SpectralBins_t *bins)
{
    not_null(bins);
    gfree(bins->w,bins->gpu_id);
    gfree(bins->tau,bins->gpu_id);
    gfree(bins->l,bins->gpu_id);
    gfree(bins->r,bins->gpu_id);
    return RS_SUCCESS;
}
