#include "omp.h"
#include "debug.h"
#include "floating_point_type.h"
#include "radiation_solvers.h"
#include "sw_flux.h"


#define sw_flux funcname(sw_flux_,TYPE)


#ifdef __NVCC__
#define MAX_NLEVELS 100

__global__
void calc_sw_flux(int const nlevels,
                  unsigned int const nws,
                  int const w0,
                  double const res,
                  fp_t const * const N,
                  fp_t const mu_dir,
                  fp_t const mu_dif,
                  fp_t const * const tau_gas,
                  fp_t const sfc_alpha_dir,
                  fp_t const sfc_alpha_dif,
                  fp_t const * const solar_flux,
                  fp_t const sol_flux_ratio,
                  fp_t * const flux_up,
                  fp_t * const flux_down,
                  fp_t * const tau_scatter)
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int const nlayers = nlevels - 1;

    fp_t w = w0 + tid*res;
    fp_t tau_total[MAX_NLEVELS-1];
    fp_t omega_avg[MAX_NLEVELS-1];
    fp_t g_avg[MAX_NLEVELS-1];
    int j;
    for (j=0;j<nlayers;j++)
    {
        enum mechanisms
        {
            ABSORB = 0,
            RAYLEIGH,
            NMECHS
        };
        fp_t omega[NMECHS];
        fp_t g[NMECHS];
        fp_t tau[NMECHS];

        /*Absorption.*/
        omega[ABSORB] = 0.;
        g[ABSORB] = 0.;
        tau[ABSORB] = tau_gas[j*nws+tid];

        /*Rayleigh scattering.*/
        omega[RAYLEIGH] = 1.;
        g[RAYLEIGH] = 0.;
        fp_t const W = w*1.e-4;
        tau[RAYLEIGH] = (N[j]*1.e-20*W*W*W*W)/
                        (0.268675*1.e5*(9.38076E2 - 10.8426*W*W));
        tau_scatter[j*nws+tid] += tau[RAYLEIGH];

        /*Get totals/averages.*/
        tau_total[j] = 0.;
        omega_avg[j] = 0.;
        g_avg[j] = 0.;
        int k;
        for (k=0;k<NMECHS;k++)
        {
            tau_total[j] += tau[k];
            omega_avg[j] += omega[k]*tau[k];
            g_avg[j] += g[k]*omega[k]*tau[k];
        }
        g_avg[j] /= omega_avg[j];
        omega_avg[j] /= tau_total[j];
    }

    fp_t flux_up_buf[MAX_NLEVELS];
    fp_t flux_down_buf[MAX_NLEVELS];
    sw_flux(nlevels,
            omega_avg,
            g_avg,
            tau_total,
            mu_dir,
            mu_dif,
            sfc_alpha_dir,
            sfc_alpha_dif,
            solar_flux[tid]*sol_flux_ratio,
            flux_up_buf,
            flux_down_buf);
    for (j=0;j<nlevels;++j)
    {
        flux_up[j*nws+tid] = flux_up_buf[j];
        flux_down[j*nws+tid] = flux_down_buf[j];
    }
}
#endif


int calc_sw_flux_h(int const nlevels,
                   unsigned int const nws,
                   int const w0,
                   double const res,
                   fp_t const * const N,
                   fp_t const mu_dir,
                   fp_t const mu_dif,
                   fp_t const * const tau_gas,
                   fp_t const sfc_alpha_dir,
                   fp_t const sfc_alpha_dif,
                   fp_t const * const solar_flux,
                   fp_t const sol_flux_ratio,
                   fp_t * const flux_up,
                   fp_t * const flux_down,
                   fp_t * const tau_scatter)
{
    int const nlayers = nlevels - 1;
    int i;

/*
#pragma omp parallel for schedule(static) \
                         default(none) \
                         private(i)
*/

    for (i=0;i<nws;++i)
    {
        fp_t w = w0 + i*res;
        fp_t tau_total[nlayers];
        fp_t omega_avg[nlayers];
        fp_t g_avg[nlayers];
        int j;
        for (j=0;j<nlayers;j++)
        {
            enum mechanisms
            {
                ABSORB = 0,
                RAYLEIGH,
                NMECHS
            };
            fp_t omega[NMECHS];
            fp_t g[NMECHS];
            fp_t tau[NMECHS];

            /*Absorption.*/
            omega[ABSORB] = 0.;
            g[ABSORB] = 0.;
            tau[ABSORB] = tau_gas[j*nws+i];

            /*Rayleigh scattering.*/
            omega[RAYLEIGH] = 1.;
            g[RAYLEIGH] = 0.;
            fp_t const W = w*1.e-4;
            tau[RAYLEIGH] = (N[j]*1.e-20*W*W*W*W)/
                            (0.268675*1.e5*(9.38076E2 - 10.8426*W*W));
            tau_scatter[j*nws+i] += tau[RAYLEIGH];

            /*Get totals/averages.*/
            tau_total[j] = 0.;
            omega_avg[j] = 0.;
            g_avg[j] = 0.;
            int k;
            for (k=0;k<NMECHS;k++)
            {
                tau_total[j] += tau[k];
                omega_avg[j] += omega[k]*tau[k];
                g_avg[j] += g[k]*omega[k]*tau[k];
            }
            g_avg[j] /= omega_avg[j];
            omega_avg[j] /= tau_total[j];
        }

        fp_t flux_up_buf[nlevels];
        fp_t flux_down_buf[nlevels];
        rs_check(sw_flux(nlevels,
                         omega_avg,
                         g_avg,
                         tau_total,
                         mu_dir,
                         mu_dif,
                         sfc_alpha_dir,
                         sfc_alpha_dif,
                         solar_flux[i]*sol_flux_ratio,
                         flux_up_buf,
                         flux_down_buf));
        for (j=0;j<nlevels;++j)
        {
            flux_up[j*nws+i] = flux_up_buf[j];
            flux_down[j*nws+i] = flux_down_buf[j];
        }
    }
    return SUCCESS;
}
