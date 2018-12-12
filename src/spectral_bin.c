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
HOST DEVICE static int bin_quad_interp(fp_t const * const x, /**< Wavenumbers [1/cm] in each bin (n,ppb).*/
                                       fp_t const * const y, /**< Optical depths in each bin (n,ppb).*/
                                       uint64_t const left, /**< Index of the left-most spectral
                                                                 grid point in each bin (n).*/
                                       uint64_t const right, /**< Index of the right-most spectral
                                                                  grid point in each bin (n).*/
                                       fp_t const w0, /**< Lower bound [1/cm] of the spectral grid.*/
                                       fp_t const wres, /**< Resolution [1/cm] of the spectral grid.*/
                                       fp_t * const tau /**< Optical depths.*/
                                      )
{
    not_null(x);
    not_null(y);
    not_null(tau);
    uint64_t j;
    for (j=left;j<=right;++j)
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
    return SUCCESS;
}


/** @brief Copy optical depths from the coarse mesh to the fine mesh in a
           bin.
    @return SUCCESS or an error code.*/
HOST DEVICE static int bin_no_interp(uint64_t const left,
                                     uint64_t const right,
                                     fp_t const * const taub,
                                     fp_t * const tau)
{
    not_null(taub);
    not_null(tau);
    uint64_t j;
    for (j=left;j<=right;++j)
    {
        tau[j] += taub[j-left];
    }
    return SUCCESS;
}


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
    return SUCCESS;
}


/*Free arrays.*/
int destroy_spectral_bins(SpectralBins_t *bins)
{
    not_null(bins);
    gfree(bins->w,bins->gpu_id);
    gfree(bins->tau,bins->gpu_id);
    gfree(bins->l,bins->gpu_id);
    gfree(bins->r,bins->gpu_id);
    return SUCCESS;
}


/*Do a quadratic interpolation of line wing values in each bin for each layer
  to the spectral grid.*/
interpolate(SpectralBins_t const * const bins,
            fp_t * const tau)
{
    not_null(bins);
    not_null(tau);
    int i;
    uint64_t j;

    /*Handle all but the last bin.*/
    if (bins->do_interp)
    {
#pragma omp parallel for collapse(2) default(none) private(i,j)
        for (i=0;i<bins->num_layers;++i)
        {
            for (j=0;j<bins->n-1;++j)
            {
                fp_t *t = &(tau[i*bins->num_wpoints]);
                fp_t const *x = &(bins->w[j*NIP]);
                fp_t const *y = &(bins->tau[i*bins->isize + j*NIP]);
                bin_quad_interp(x,
                                y,
                                bins->l[j],
                                bins->r[j],
                                bins->w0,
                                bins->wres,
                                t);
            }
        }
    }
    else
    {
#pragma omp parallel for collapse(2) default(none) private(i,j)
        for (i=0;i<bins->num_layers;++i)
        {
            for (j=0;j<bins->n-1;++j)
            {
                fp_t *t = &(tau[i*bins->num_wpoints]);
                fp_t const *y = &(bins->tau[i*bins->isize + j*NIP]);
                bin_no_interp(bins->l[j],
                              bins->r[j],
                              y,
                              t);
            }
        }
    }

    /*Handle the last bin.*/
    j = bins->n - 1;
    if (bins->do_last_interp)
    {
#pragma omp parallel for default(none) shared(j) private(i)
        for (i=0;i<bins->num_layers;++i)
        {
            fp_t *t = &(tau[i*bins->num_wpoints]);
            fp_t const *x = &(bins->w[j*NIP]);
            fp_t const *y = &(bins->tau[i*bins->isize + j*NIP]);
            bin_quad_interp(x,
                            y,
                            bins->l[j],
                            bins->r[j],
                            bins->w0,
                            bins->wres,
                            t);
        }
    }
    else
    {
#pragma omp parallel for default(none) shared(j) private(i)
        for (i=0;i<bins->num_layers;++i)
        {
            fp_t *t = &(tau[i*bins->num_wpoints]);
            fp_t const *y = &(bins->tau[i*bins->isize + j*NIP]);
            bin_no_interp(bins->l[j],
                          bins->r[j],
                          y,
                          t);
        }
    }
    return SUCCESS;
}
