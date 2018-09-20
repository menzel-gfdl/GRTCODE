#ifndef DEVICE_LAUNCH_H_
#define DEVICE_LAUNCH_H_

#include <stdint.h>
#include "floating_point_type.h"
#include "ozone_continuum.h"
#include "parse_HITRAN_file.h"
#include "water_vapor_continuum.h"


int launch(int const num_levels,
           fp_t const * const P,
           fp_t const * const T,
           fp_t const * const x,
           fp_t * const Pavg,
           fp_t * const Tavg,
           fp_t * const N,
           fp_t * const Ns,
           fp_t * const Psavg,
           fp_t * const gamma,
           fp_t * const Pshift,
           fp_t * const s,
           LineParams_t *lines,
           uint64_t const molecule_bit_field,
           LineParams_t ** const line_params,
           double const w0,
           double const wres,
           uint64_t const num_wpoints,
           double const wcutoff,
           int const use_h2o_ctm,
           WaterVaporContinuumCoefs_t * const h2o_cc,
           int const use_o3_ctm,
           OzoneContinuumCoefs_t const * const o3_cc,
           fp_t * const tau);


#endif
