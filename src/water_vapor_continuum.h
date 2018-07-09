#ifndef WATER_VAPOR_CONTINUUM_H_
#define WATER_VAPOR_CONTINUUM_H_

#include <stdint.h>
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
    uint64_t num_wpoints;
} WaterVaporContinuumCoefs_t;


int get_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc,
                                    uint64_t const num_wpoints,
                                    double const w0,
                                    double const res);


int free_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc);


int put_water_vapor_coefs_on_device(WaterVaporContinuumCoefs_t const * const in,
                                    WaterVaporContinuumCoefs_t * const out);


int remove_water_vapor_coefs_from_device(WaterVaporContinuumCoefs_t * const in);


#ifdef __NVCC__
__global__
void calc_water_vapor_ctm_optical_depth(uint64_t const num_wpoints,
                                        int const num_layers,
                                        fp_t * const tau,
                                        fp_t const * const CS,
                                        fp_t const * const T,
                                        fp_t const * const Ps,
                                        fp_t const * const N,
                                        fp_t const * const T0,
                                        fp_t const * const CF,
                                        fp_t const * const P,
                                        fp_t const * const T0F);
#endif


void calc_water_vapor_ctm_optical_depth_h(uint64_t const num_wpoints,
                                          int const num_layers,
                                          fp_t * const tau,
                                          fp_t const * const CS,
                                          fp_t const * const T,
                                          fp_t const * const Ps,
                                          fp_t const * const N,
                                          fp_t const * const T0,
                                          fp_t const * const CF,
                                          fp_t const * const P,
                                          fp_t const * const T0F);


#endif
