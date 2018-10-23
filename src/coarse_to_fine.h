#ifndef COARSE_TO_FINE_H_
#define COARSE_TO_FINE_H_

#include <stdint.h>
#include "floating_point_type.h"


#ifdef __NVCC__
__global__
void coarse_to_fine(int const num_layers,
                    double const w0,
                    double const wres_fine,
                    double const wres_coarse,
                    uint64_t const num_wpoints_fine,
                    uint64_t const num_wpoints_coarse,
                    fp_t * const tau_fine,
                    fp_t const * const tau_coarse);
#endif


void coarse_to_fine_h(int const num_layers,
                      double const w0,
                      double const wres_fine,
                      double const wres_coarse,
                      uint64_t const num_wpoints_fine,
                      uint64_t const num_wpoints_coarse,
                      fp_t * const tau_fine,
                      fp_t const * const tau_coarse);


#endif
