#ifndef LAUNCH_H_
#define LAUNCH_H_

#include "floating_point_type.h"
#include "model_fields.h"
#include "parseHITRANfile.h"
#include "output_fields.h"


typedef struct WorkVars
{
    fp_t *P;
    fp_t *T;
    fp_t *TSURF;
    fp_t *EMIS;
    fp_t *x;
    fp_t *Pavg;
    fp_t *Tavg;
    fp_t *N;
    fp_t *Psavg;
    fp_t *GAMMA;
    fp_t *PSHIFT;
    fp_t *S;
    line_params_t *LINES;
    fp_t *Snn_ref; /*Only used by host.*/
    fp_t *tau;
    fp_t *lw_flux_down_per_w;
    fp_t *lw_flux_up_per_w;
    fp_t *lw_flux_down;
    fp_t *lw_flux_up;
} WorkVars_t;


int alloc_work_vars(WorkVars_t *vars,
                    int const nlevels,
                    int const nlines,
                    int const nws,
                    int const launch_type);


int free_work_vars(WorkVars_t *vars,
                   int const launch_type);


int launch_host(WorkVars_t * const vars,
                req_model_fields_t * const input_data,
                int const time,
                int const lon,
                int const lat,
                int const nmols,
                line_params_t ** const line_params,
                unsigned int nws,
                fp_t const w,
                double const res,
                int const breadth,
                int const continuum,
                OutputFields_t * const output_data);


#endif
