#ifndef WATER_VAPOR_CONTINUUM_H_
#define WATER_VAPOR_CONTINUUM_H_

#include "floating_point_type.h"


enum water_vapor_coefs
{
    MTCKD25_F296 = 0,
    MTCKD25_S296,
    CKDF,
    CKDS,
    NUM_COEFS
};


typedef struct WaterVaporContinuumCoefs
{
    fp_t **coefs;
} WaterVaporContinuumCoefs_t;


int get_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc,
                                    unsigned int const nws,
                                    int const w0,
                                    double const res,
                                    int put_on_device);


int free_water_vapor_continuum_coeffs(WaterVaporContinuumCoefs_t *cc,
                                      int const on_device);


#ifdef __NVCC__
__global__
void calc_water_vapor_ctm_optdepth(unsigned int const nF,
                                   int const numLayers,
                                   fp_t * const optdepth,
                                   fp_t const * const CS,
                                   fp_t const * const T,
                                   fp_t const * const Ps,
                                   fp_t const * const N,
                                   fp_t const * const T0,
                                   fp_t const * const CF,
                                   fp_t const * const P,
                                   fp_t const * const T0F);
#endif


void calc_water_vapor_ctm_optdepth_h(unsigned int const nF,
                                     int const numLayers,
                                     fp_t * const optdepth,
                                     fp_t const * const CS,
                                     fp_t const * const T,
                                     fp_t const * const Ps,
                                     fp_t const * const N,
                                     fp_t const * const T0,
                                     fp_t const * const CF,
                                     fp_t const * const P,
                                     fp_t const * const T0F);


#endif
