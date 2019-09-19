#include <stdint.h>
#include "curtis_godson.h"
#include "debug.h"
#include "extern.h"
#include "floating_point_type.h"
#include "optics.h"
#include "rayleigh.h"
#include "rs_config.h"


/** @brief Calculate the single-scattering albedo, asymmetric factor, and optical depth
           contribution for Rayleigh scattering of an atmospheric layer at an input
           wavenumber.
    @return RS_SUCCESS or an error code.*/
HOST DEVICE static int rayleigh_kernel(double const w, /**< Wavenumber [1/cm].*/
                                       fp_t const n, /**< Integrated number density [cm^-2].*/
                                       fp_t * const omega, /**< Single-scattering albedo.*/
                                       fp_t * const g, /**< Asymmetric factor.*/
                                       fp_t * const tau /**< Scattering optical depth.*/
                                      )
{
    *omega = 1.;
    *g = 0.;
    fp_t const W = w*1.e-4;
    *tau = (n*1.e-20*W*W*W*W)/(0.268675*1.e5*(9.38076E2 - 10.8426*W*W));
    return RS_SUCCESS;
}


/** @brief Calculate the optical properties due to rayleigh scattering.
    @return RS_SUCCESS or an error code.*/
static int rayleigh(int const num_layers, /**< Number of atmospheric pressure layers.*/
                    double const w0, /**< Spectral grid lower bound [1/cm].*/
                    double const dw, /**< Spectral grid resolution [1/cm].*/
                    uint64_t const num_wpoints, /**< Spectral grid size.*/
                    fp_t const * const n, /**< Integrated number density [cm^-2] of
                                               each atmospheric layer.*/
                    fp_t * const omega, /**< Single-scattering albedo of each
                                             atmospheric layer at each spectral
                                             grid point.*/
                    fp_t * const g, /**< Asymmetric factor of each atmospheric
                                         layer at each spectral grid point.*/
                    fp_t * const tau /**< Optical depth of each atmospheric
                                          layer at each spectral grid point.*/
                   )
{
    int i;
    uint64_t j;
#pragma omp parallel for collapse(2) default(none) private(i,j)
    for (i=0; i<num_layers; ++i)
    {
        for (j=0; j<num_wpoints; ++j)
        {
            double const w = w0 + j*dw;
            uint64_t offset = i*num_wpoints + j;
            rayleigh_kernel(w, n[i], &(omega[offset]), &(g[offset]), &(tau[offset]));
        }
    }
    return RS_SUCCESS;
}


#ifdef __NVCC__
/** @brief Calculate the optical properties due to rayleigh scattering.*/
__global__ static void rayleigh_d(int const num_layers, /**< Number of atmospheric pressure layers.*/
                                  double const w0, /**< Spectral grid lower bound [1/cm].*/
                                  double const dw, /**< Spectral grid resolution [1/cm].*/
                                  uint64_t const num_wpoints, /**< Spectral grid size.*/
                                  fp_t const * const n, /**< Integrated number density [cm^-2] of
                                                             each atmospheric layer.*/
                                  fp_t * const omega, /**< Single-scattering albedo of each
                                                           atmospheric layer at each spectral
                                                           grid point.*/
                                  fp_t * const g, /**< Asymmetric factor of each atmospheric
                                                       layer at each spectral grid point.*/
                                  fp_t * const tau /**< Optical depth of each atmospheric
                                                        layer at each spectral grid point.*/
                                 )
{
    uint64_t i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < num_wpoints)
    {
        double const w = w0 + i*dw;
        int j;
        for (j=0; j<num_layers; ++j)
        {
            uint64_t offset = j*num_wpoints + i;
            rayleigh_kernel(w, n[j], &(omega[offset]), &(g[offset]), &(tau[offset]));
        }
    }
    return;
}
#endif


/*Calculate the optical properties due to Rayleigh scattering.*/
EXTERN int rayleigh_scattering(Optics_t * const optics, fp_t * const pressure)
{
    not_null(optics);
    not_null(pressure);
    fp_t const mbtoatm = 0.000986923f;
    int num_levels = optics->num_layers + 1;
    fp_t p_atm[MAX_NUM_LEVELS];
    int i;
    for (i=0; i<num_levels; ++i)
    {
        p_atm[i] = pressure[i]*mbtoatm;
    }
    fp_t *p;
    fp_t number_density[MAX_NUM_LAYERS];
    fp_t *n;
    if (optics->device == HOST_ONLY)
    {
        p = p_atm;
        n = number_density;
    }
    else
    {
        /*Place data on device.*/
        gmalloc(p, num_levels, optics->device);
        gmemcpy(p, p_atm, num_levels, optics->device, FROM_HOST);
        gmalloc(n, optics->num_layers, optics->device);
    }

    /*Calculate integrated number densities [cm^2] in each atmospheric layer.*/
    glaunch(calc_number_densities, optics->num_layers, optics->device,
            optics->num_layers, p, n);

    /*Calculate optical properties for Rayleigh scattering.*/
    glaunch(rayleigh, optics->grid.n, optics->device, optics->num_layers,
            optics->grid.w0, optics->grid.dw, optics->grid.n, n, optics->omega,
            optics->g, optics->tau);

    /*Clean up.*/
    if (optics->device != HOST_ONLY)
    {
        gfree(p, optics->device);
        gfree(n, optics->device);
    }
    return RS_SUCCESS;
}
