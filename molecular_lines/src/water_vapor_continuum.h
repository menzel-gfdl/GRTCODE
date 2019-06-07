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
    int gpu_id;
} WaterVaporContinuumCoefs_t;


int get_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc,
                                    char const * const h2o_ctm_dir,
                                    uint64_t const num_wpoints,
                                    double const w0,
                                    double const res,
                                    int const gpu_id);


int free_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc);


#endif
