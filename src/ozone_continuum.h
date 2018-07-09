#ifndef OZONE_CONTINUUM_H_
#define OZONE_CONTINUUM_H_

#include <stdint.h>
#include "floating_point_type.h"


typedef struct OzoneContinuumCoefs
{
    fp_t *cross_section; /*Ozone continuum cross-section [cm^2].*/
    uint64_t num_wpoints;
} OzoneContinuumCoefs_t;


/*Read in the ozone continuum coefficients.*/
int get_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc,
                              uint64_t const num_wpoints,
                              double const w0,
                              double const res);


int free_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc);


int put_ozone_coefs_on_device(OzoneContinuumCoefs_t const * const in,
                              OzoneContinuumCoefs_t * const out);


int remove_ozone_coefs_from_device(OzoneContinuumCoefs_t * const in);


#ifdef __NVCC__
__global__
void calc_ozone_ctm_optical_depth(uint64_t const num_wpoints,
                                  int const num_layers,
                                  fp_t const * const cross_section,
                                  fp_t const * const N,
                                  fp_t * const tau);
#endif


void calc_ozone_ctm_optical_depth_h(uint64_t const nws,
                                    int const num_layers,
                                    fp_t const * const cross_section,
                                    fp_t const * const N,
                                    fp_t * const tau);


#endif
