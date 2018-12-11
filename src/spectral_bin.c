#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "floating_point_type.h"
#include "spectral_bin.h"
#include "utils.h"


/** @brief Loop through the bins, interpolate optical depth values, and add
           them to the input optical depth array.
    @return SUCCESS or an error code.*/
static int quad_bin_interp(fp_t const * const wb, /**< Wavenumbers [1/cm] in each bin (n,ppb).*/
                           fp_t const * const taub, /**< Optical depths in each bin (n,ppb).*/
                           uint64_t const * const left, /**< Index of the left-most spectral
                                                             grid point in each bin (n).*/
                           uint64_t const * const right, /**< Index of the right-most spectral
                                                              grid point in each bin (n).*/
                           int const do_interp, /**< Flag telling if interpolation is
                                                     required.*/
                           uint64_t const num_bins, /**< The number of bins.*/
                           fp_t const w0, /**< Lower bound [1/cm] of the spectral grid.*/
                           fp_t const wres, /**< Resolution [1/cm] of the spectral grid.*/
                           fp_t * const tau /**< Optical depths.*/
                          )
{
    if (do_interp)
    {
        uint64_t i;
#pragma omp parallel for default(none) private(i)
        for (i=0;i<num_bins;++i)
        {
            fp_t x[3];
            x[0] = wb[i*NIP];
            x[1] = wb[i*NIP+1];
            x[2] = wb[i*NIP+2];
            fp_t y[3];
            y[0] = taub[i*NIP];
            y[1] = taub[i*NIP+1];
            y[2] = taub[i*NIP+2];
            uint64_t j;
            for (j=left[i];j<=right[i];++j)
            {
                fp_t w = w0 + j*wres;
                fp_t t = (w-x[1])*(w-x[2])*y[0]/((x[0]-x[1])*(x[0]-x[2])) +
                         (w-x[0])*(w-x[2])*y[1]/((x[1]-x[0])*(x[1]-x[2])) +
                         (w-x[0])*(w-x[1])*y[2]/((x[2]-x[0])*(x[2]-x[1]));
                if (t < 0.f)
                {
                    t = 0.f;
                }
                tau[j] += t;
            }
        }
    }
    else
    {
        uint64_t i;
#pragma omp parallel for default(none) private(i)
        for (i=0;i<num_bins;++i)
        {
            uint64_t j;
            for (j=left[i];j<=right[i];++j)
            {
                tau[j] += taub[i*NIP+j-left[i]];
            }
        }
    }
    return SUCCESS;
}


/*Set parameters and allocate arrays.*/
int create_spectral_bins(SpectralBins_t *bins,
                         int const num_layers,
                         double const w0,
                         uint64_t const n,
                         double const wres,
                         double const bin_width)
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
    uint64_t *l;
    gmalloc(l,bins->n,HOST_ONLY);
    uint64_t *r;
    gmalloc(r,bins->n,HOST_ONLY);
    fp_t *w;
    gmalloc(w,bins->isize,HOST_ONLY);

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
    return SUCCESS;
}


/*Free arrays.*/
int destroy_spectral_bins(SpectralBins_t *bins)
{
    not_null(bins);
    not_null(bins->w);
    free(bins->w);
    not_null(bins->tau);
    free(bins->tau);
    not_null(bins->l);
    free(bins->l);
    not_null(bins->r);
    free(bins->r);
    return SUCCESS;
}


/*Do a quadratic interpolation of line wing values in each bin for each layer
  to the spectral grid.*/
int interpolate(SpectralBins_t const * const bins,
                fp_t * const tau)
{
    not_null(bins);
    uint64_t const *l = bins->l;
    uint64_t const *r = bins->r;
    fp_t const *w = bins->w;
    int i;
    for (i=0;i<bins->num_layers;++i)
    {
        fp_t const *taub = &(bins->tau[i*bins->isize]);
        fp_t *t = &(tau[i*bins->num_wpoints]);

        /*Do the interpolation on all but the last bin.*/
        throw(quad_bin_interp(w,
                              taub,
                              l,
                              r,
                              bins->do_interp,
                              bins->n-1,
                              bins->w0,
                              bins->wres,
                              t));

        /*Move pointers so that they reference the last bin.*/
        uint64_t o = bins->n-1;
        uint64_t const *llast = l + o;
        uint64_t const *rlast = r + o;
        fp_t const *wlast = w + o*NIP;
        taub = taub + o*NIP;

        /*Do the interpolation on the last bin.*/
        throw(quad_bin_interp(wlast,
                              taub,
                              llast,
                              rlast,
                              bins->do_last_interp,
                              1,
                              bins->w0,
                              bins->wres,
                              t));
    }
    return SUCCESS;
}
