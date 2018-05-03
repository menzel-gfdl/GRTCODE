#ifndef HOST_LAUNCH_H_
#define HOST_LAUNCH_H_

#include "floating_point_type.h"
#include "model_fields.h"
#include "parseHITRANfile.h"
#include "ozone_continuum.h"
#include "output_fields.h"
#include "solar_flux.h"
#include "water_vapor_continuum.h"

typedef struct WorkVars_h
{
    fp_t *P;
    fp_t *T;
    fp_t *x;
    fp_t *Pavg;
    fp_t *Tavg;
    fp_t *N;
    fp_t *Ns;
    fp_t *Psavg;
    fp_t *GAMMA;
    fp_t *PSHIFT;
    fp_t *S;
    line_params_t *LINES;
    fp_t *Snn_ref; /*Only used by host.*/
    fp_t *tau_gas;
    fp_t *tau_scatter;
    fp_t *lw_flux_down_per_w;
    fp_t *lw_flux_up_per_w;
    fp_t *lw_flux_down;
    fp_t *lw_flux_up;
    fp_t *sw_flux_down_per_w;
    fp_t *sw_flux_up_per_w;
    fp_t *sw_flux_down;
    fp_t *sw_flux_up;
} WorkVars_h_t;


int alloc_work_vars_h(WorkVars_h_t *vars,
                      int const nlevels,
                      int const nlines,
                      int const nws);


int free_work_vars_h(WorkVars_h_t *vars);


int launch_h(WorkVars_h_t * const vars,
             req_model_fields_t * const input_data,
             SolarFlux_t const * const solar_flux,
             int const time,
             int const lon,
             int const lat,
             int const nmols,
             line_params_t ** const line_params,
             unsigned int nws,
             fp_t const w,
             double const res,
             int const breadth,
             int const h2o_ctm,
             WaterVaporContinuumCoefs_t * const h2o_continuum,
             int const o3_ctm,
             OzoneContinuumCoefs_t const * const o3_continuum,
             OutputFields_t * const output_data);


#endif
