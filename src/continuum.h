#ifndef CONTINUUM_H_
#define CONTINUUM_H_

#include "floating_point_type.h"

typedef struct ContinuumCoefs
{
    fp_t **coefs; /*Continuum coefficients (wavenumber) [1/cm].*/
} ContinuumCoefs_t;


int get_h2o_continuum_coefs(ContinuumCoefs_t *h2o,
                            unsigned int const nws,
                            int const w,
                            double const res,
                            int put_on_device);


int free_continuum_coeffs(ContinuumCoefs_t *c,
                          int const on_device);


#ifdef __NVCC__
__global__
void calc_ctm_optdepth(unsigned int const nF,
                       int const numLayers,
                       fp_t * const optdepth,
                       fp_t const * const CS,
                       fp_t const * const T,
                       fp_t const * const PS_H2O,
                       fp_t const * const Z,
                       fp_t const * const T0,
                       fp_t const * const CF,
                       fp_t const * const P,
                       fp_t const * const T0F);
#endif


void calc_ctm_optdepth_h(unsigned int const nF,
                         int const numLayers,
                         fp_t * const optdepth,
                         fp_t const * const CS,
                         fp_t const * const T,
                         fp_t const * const PS_H2O,
                         fp_t const * const Z,
                         fp_t const * const T0,
                         fp_t const * const CF,
                         fp_t const * const P,
                         fp_t const * const T0F);


#endif
