#ifndef SW_FLUX_H_
#define SW_FLUX_H_

#include "floating_point_type.h"


#ifdef __NVCC__
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
                  fp_t * const tau_scatter);
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
                   fp_t * const tau_scatter);


#endif
