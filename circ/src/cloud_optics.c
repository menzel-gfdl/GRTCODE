#include <math.h>
#include <stdint.h>
#include "debug.h"
#include "extern.h"
#include "floating_point_type.h"
#include "hu_stamnes_1993_tables.h"
#include "optics.h"
#include "spectral_grid.h"
#include "utilities.h"


/** @brief Calculate cloud optics using the parameterization described in
           https://doi.org/10.1175/1520-0442(1993)006<0728:AAPOTR>2.0.CO;2
    @return RS_SUCCESS or an error code.*/
static int hu_stamnes_1993(fp_t const liquid_water_path, /**< Liquid water path [g/m^3].*/
                           fp_t const droplet_equivalent_radius, /**< Droplet equivalent radius [microns].*/
                           SpectralGrid_t const grid, /**< Spectral grid object.*/
                           fp_t * const optical_depth, /**< Optical depth (wavenumber).*/
                           fp_t * const single_scatter_albedo, /**< Single-scatter albedo (wavenumber).*/
                           fp_t * const asymmetry_factor /**< Asymmetry factor (wavenumber).*/
                          )
{
    int const num_radius_bins = 3;
    fp_t const radius_bin_bounds[4] = {2.5, 12.5, 30., 60.};
    fp_t const min_radius = 3.;
    fp_t const max_radius = 60.;
    fp_t r = droplet_equivalent_radius;
    if (r < min_radius)
    {
        r = min_radius;
    }
    else if (r > max_radius)
    {
        r = max_radius;
    }
    int i = num_radius_bins - 1;
    int j;
    for (j=0; j<num_radius_bins; ++j)
    {
        if (r < radius_bin_bounds[j+1])
        {
            i = j;
            break;
        }
    }
    i *= num_spectral_bands;
    fp_t tau[num_spectral_bands+2];
    fp_t omega[num_spectral_bands+2];
    fp_t g[num_spectral_bands+2];
    fp_t const lwp = liquid_water_path/1.e3;
    for (j=0; j<num_spectral_bands; ++j)
    {
        /*Equations 13, 14, and 15.*/
        tau[j+1] = lwp*(a1[i+j]*pow(r, b1[i+j]) + c1[i+j]);
        omega[j+1] = 1. - (a2[i+j]*pow(r, b2[i+j]) + c2[i+j]);
        if (omega[j+1] < 0.)
        {
            /*Some longwave bands can return negative single-scatter albedos.*/
            omega[j+1] = 0.;
        }
        g[j+1] = a3[i+j]*pow(r, b3[i+j]) + c3[i+j];
    }
    tau[0] = tau[1];
    tau[num_spectral_bands+1] = tau[num_spectral_bands];
    omega[0] = omega[1];
    omega[num_spectral_bands+1] = omega[num_spectral_bands];
    g[0] = g[1];
    g[num_spectral_bands+1] = g[num_spectral_bands];

    /*Convert band values to wavenumber [1/cm];*/
    fp_t w[num_spectral_bands+2];
    for (j=0;j<num_spectral_bands+2;++j)
    {
        w[j] = 1.e4/spectral_band[j];
    }

    /*Interpolate from bands to the input spectral grid.*/
    fp_t *x = NULL;
    gmalloc(x, grid.n, HOST_ONLY);
    uint64_t n;
    for (n=0; n<grid.n; ++n)
    {
        x[n] = (fp_t)(grid.w0 + n*grid.dw);
    }
    int e = interpolate2(w, tau, num_spectral_bands+2, x, optical_depth, (size_t)grid.n,
                         angstrom_exponent_sample, NULL);
    if (e == RS_VALUE_ERR)
    {
        gmemset(optical_depth, 0, grid.n, HOST_ONLY);
    }
    else if (e != RS_SUCCESS)
    {
        catch(e);
    }
    catch(interpolate2(w, omega, num_spectral_bands+2, x, single_scatter_albedo,
                       (size_t)grid.n, linear_sample, NULL));
    catch(interpolate2(w, g, num_spectral_bands+2, x, asymmetry_factor,
                       (size_t)grid.n, linear_sample, NULL));
    gfree(x, HOST_ONLY);
    return RS_SUCCESS;
}


EXTERN int cloud_optics(Optics_t * const optics, fp_t * const liquid_water_path,
                        fp_t * const droplet_equivalent_radius)
{
    not_null(optics);
    not_null(liquid_water_path);
    not_null(droplet_equivalent_radius);
    uint64_t n = optics->grid.n;
    fp_t *tau;
    gmalloc(tau, n, HOST_ONLY);
    fp_t *omega;
    gmalloc(omega, n, HOST_ONLY);
    fp_t *g;
    gmalloc(g, n, HOST_ONLY);
    int i;
    for (i=0; i<optics->num_layers; ++i)
    {
        catch(hu_stamnes_1993(liquid_water_path[i], droplet_equivalent_radius[i],
                              optics->grid, tau, omega, g));
        gmemcpy(&(optics->tau[i*n]), tau, n, optics->device, FROM_HOST);
        gmemcpy(&(optics->omega[i*n]), omega, n, optics->device, FROM_HOST);
        gmemcpy(&(optics->g[i*n]), g, n, optics->device, FROM_HOST);
    }
    gfree(tau, HOST_ONLY);
    gfree(omega, HOST_ONLY);
    gfree(g, HOST_ONLY);
    return RS_SUCCESS;
}
