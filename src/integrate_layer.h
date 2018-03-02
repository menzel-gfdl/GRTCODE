#ifndef INTEGRATE_LAYER_H_
#define INTEGRATE_LAYER_H_

#include "floating_point_type.h"


#ifdef __NVCC__
void get_avg_TP(int const nLayer,
                fp_t const * const P,
                fp_t const * const T,
                fp_t * const Pavg,
                fp_t * const Tavg);
#endif


void get_avg_TP_h(int const nLayer,
                  fp_t const * const P,
                  fp_t const * const T,
                  fp_t * const Pavg,
                  fp_t * const Tavg);


#ifdef __NVCC__
__global__
void get_avg_NPs(int const nlayer,
                 fp_t const * const x,
                 fp_t const * const P,
                 fp_t * const N,
                 fp_t * const Psavg);
#endif


void get_avg_NPs_h(int const nlayer,
                   fp_t const * const x,
                   fp_t const * const P,
                   fp_t * const N,
                   fp_t * const Psavg);


#endif
