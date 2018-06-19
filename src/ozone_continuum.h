#ifndef OZONE_CONTINUUM_H_
#define OZONE_CONTINUUM_H_

#include "floating_point_type.h"


typedef struct OzoneContinuumCoefs
{
    fp_t *cross_section; /*Ozone continuum cross-section [cm^2].*/
    int nws;
} OzoneContinuumCoefs_t;


/*Read in the ozone continuum coefficients.*/
int get_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc,
                              unsigned int const nws,
                              int const w0,
                              double const res);


int free_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc);


int put_ozone_coefs_on_device(OzoneContinuumCoefs_t const * const in,
                              OzoneContinuumCoefs_t * const out);


int remove_ozone_coefs_from_device(OzoneContinuumCoefs_t * const in);


#ifdef __NVCC__
__global__
void calc_ozone_ctm_optdepth(unsigned int const nws,
                             int const nlayers,
                             fp_t const * const cross_section,
                             fp_t const * const N,
                             fp_t * const tau);
#endif


void calc_ozone_ctm_optdepth_h(unsigned int const nws,
                               int const nlayers,
                               fp_t const * const cross_section,
                               fp_t const * const N,
                               fp_t * const tau);


#endif
