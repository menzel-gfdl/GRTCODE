#ifndef EVAL_PROFILE_H_
#define EVAL_PROFILE_H_

#include "floating_point_type.h"


#ifdef __NVCC__
__global__ void eval_profile(int const mol_id,
                             unsigned int const num_lines,
                             uint64_t const num_wpoints,
                             double const w0,
                             double const wres,
                             int const num_layers,
                             double const wcutoff,
                             fp_t const * const T,
                             fp_t const * const gamma,
                             fp_t const * const Pshift,
                             fp_t const * const s,
                             fp_t const * const N,
                             fp_t * const tau);
#endif


void eval_profile_h(int const mol_id,
                    unsigned int const num_lines,
                    uint64_t const num_wpoints,
                    double const w0,
                    double const wres,
                    int const num_layers,
                    double const wcutoff,
                    fp_t const * const T,
                    fp_t const * const gamma,
                    fp_t const * const Pshift,
                    fp_t const * const s,
                    fp_t const * const N,
                    fp_t * const tau);


#endif
