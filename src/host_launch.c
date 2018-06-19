#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "debug.h"
#include "eval_gamma.h"
#include "eval_profile.h"
#include "eval_pShift.h"
#include "eval_Snn_correction.h"
#include "floating_point_type.h"
#include "host_launch.h"
#include "integrate_layer.h"
#include "lw_flux.h"
#include "model_fields.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "parseHITRANfile.h"
#include "pre_eval_Snn.h"
#include "solar_flux.h"
#include "sw_flux.h"
#include "utils.h"
#include "water_vapor_continuum.h"


int alloc_work_vars_h(WorkVars_h_t *vars,
                      int const nlevels,
                      int const nlines,
                      int const nws)
{
    not_null(vars);
    int const nlayers = nlevels - 1;
    int const l = nlayers*nlines;
    int const m = nlevels*nws;

    /*Will point to already allocated input/output arrays.*/
    vars->P = NULL;
    vars->T = NULL;
    vars->x = NULL;
    vars->LINES = NULL;
    vars->tau_gas = NULL;
    vars->tau_scatter = NULL;
    vars->lw_flux_down = NULL;
    vars->lw_flux_up = NULL;
    vars->sw_flux_down = NULL;
    vars->sw_flux_up = NULL;

    /*Per-column buffers.*/
    check(malloc_ptr((void **)(&(vars->Pavg)),
                     sizeof(*(vars->Pavg))*nlayers));

    check(malloc_ptr((void **)(&(vars->Tavg)),
                     sizeof(*(vars->Tavg))*nlayers));

    check(malloc_ptr((void **)(&(vars->N)),
                     sizeof(*(vars->N))*nlayers));

    check(malloc_ptr((void **)(&(vars->Ns)),
                     sizeof(*(vars->Ns))*nlayers));

    check(malloc_ptr((void **)(&(vars->Psavg)),
                     sizeof(*(vars->Psavg))*nlayers));

    check(malloc_ptr((void **)(&(vars->GAMMA)),
                     sizeof(*(vars->GAMMA))*l));

    check(malloc_ptr((void **)(&(vars->PSHIFT)),
                     sizeof(*(vars->PSHIFT))*l));

    check(malloc_ptr((void **)(&(vars->S)),
                     sizeof(*(vars->S))*l));

    check(malloc_ptr((void **)(&(vars->Snn_ref)),
                     sizeof(*(vars->Snn_ref))*nlines));

    check(malloc_ptr((void **)(&(vars->lw_flux_down_per_w)),
                     sizeof(*(vars->lw_flux_down_per_w))*m));

    check(malloc_ptr((void **)(&(vars->lw_flux_up_per_w)),
                     sizeof(*(vars->lw_flux_up_per_w))*m));

    check(malloc_ptr((void **)(&(vars->sw_flux_down_per_w)),
                     sizeof(*(vars->sw_flux_down_per_w))*m));

    check(malloc_ptr((void **)(&(vars->sw_flux_up_per_w)),
                     sizeof(*(vars->sw_flux_up_per_w))*m));
    return SUCCESS;
}


int free_work_vars_h(WorkVars_h_t *vars)
{
    not_null(vars);
    free(vars->Pavg);
    free(vars->Tavg);
    free(vars->N);
    free(vars->Ns);
    free(vars->Psavg);
    free(vars->GAMMA);
    free(vars->PSHIFT);
    free(vars->S);
    free(vars->Snn_ref);
    free(vars->lw_flux_down_per_w);
    free(vars->lw_flux_up_per_w);
    free(vars->sw_flux_down_per_w);
    free(vars->sw_flux_up_per_w);
    return SUCCESS;
}


int launch_h(WorkVars_h_t * const vars,
             req_model_fields_t * const input_data,
             SolarFlux_t const * const solar_flux,
             int const time,
             int const lon,
             int const lat,
             int const nmols,
             line_params_t ** const line_params,
             unsigned int const nws,
             fp_t const w,
             double const res,
             int const breadth,
             int const h2o_ctm,
             WaterVaporContinuumCoefs_t * const h2o_continuum,
             int const o3_ctm,
             OzoneContinuumCoefs_t const * const o3_continuum,
             OutputFields_t * const output_data)
{
    not_null(vars);
    not_null(input_data);
    not_null(solar_flux);
    not_null(line_params);
    not_null(h2o_continuum);
    not_null(o3_continuum);
    not_null(output_data);

    /*Point to output buffers.*/
    vars->tau_gas = output_data->tau_gas;
    vars->tau_scatter = output_data->tau_scatter;
    vars->lw_flux_down = output_data->lw_flux_down;
    vars->lw_flux_up = output_data->lw_flux_up;
    vars->sw_flux_down = output_data->sw_flux_down;
    vars->sw_flux_up = output_data->sw_flux_up;

    /*Zero out buffers used to accumulate results.*/
    int nlevels = input_data->nlevel;
    int nlayers = nlevels - 1;
    memset(vars->tau_gas,
           0,
           sizeof(*(vars->tau_gas))*nlayers*nws);
    memset(vars->tau_scatter,
           0,
           sizeof(*(vars->tau_scatter))*nlayers*nws);
    memset(vars->lw_flux_down,
           0,
           sizeof(*(vars->lw_flux_down))*nlevels);
    memset(vars->lw_flux_up,
           0,
           sizeof(*(vars->lw_flux_up))*nlevels);
    memset(vars->sw_flux_down,
           0,
           sizeof(*(vars->sw_flux_down))*nlevels);
    memset(vars->sw_flux_up,
           0,
           sizeof(*(vars->sw_flux_up))*nlevels);

    /*Point to the correct column of input data.*/
    unsigned int offset = time*input_data->nlon*input_data->nlat +
                          lon*input_data->nlat + lat;
    fp_t const TSURF = input_data->TSURF[offset];
    fp_t const EMIS = input_data->EMIS[offset];
    fp_t const SFC_DIR_ALB = input_data->SFC_DIR_ALB[offset];
    fp_t const SFC_DIF_ALB = input_data->SFC_DIF_ALB[offset];
    fp_t const MU_DIF = input_data->COS_DIF_BEAM_ANG;
    fp_t const MU_DIR = input_data->COS_SOL_ZEN_ANG[offset];
    fp_t const SOL_FLUX_RATIO = (input_data->TOTAL_SOL_FLUX[offset])/
                                (solar_flux->total_sw_flux);
    offset *= nlevels;
    vars->P = &(input_data->P[offset]);
    vars->T = &(input_data->T[offset]);

    /*Calculate the total number density of air moleucles integrated across
      each layer.*/
    log_mesg("Launching kernel integrated_N_h at point (%d,%d,%d).",
             time,
             lon,
             lat);
    check(integrated_N_h(nlayers,
                         vars->P,
                         vars->N));

    /*Calculate integrated average layer quantities.*/
    log_mesg("Launching kernel Curtis_Godson_PT_h at point (%d,%d,%d).",
             time,
             lon,
             lat);
    check(Curtis_Godson_PT_h(nlayers,
                             vars->P,
                             vars->T,
                             vars->Pavg,
                             vars->Tavg));

    /*Loop over the molecules and calculate the optical depths.*/
    int mol;
    for (mol=0;mol<nmols;++mol)
    {
        /*Point to the line parameters for the current molecule.*/
        vars->LINES = line_params[mol];

        /*Copy the Snn_ref array.*/
        unsigned int nlines = line_params[mol]->nLines;
        memcpy(vars->Snn_ref,
               line_params[mol]->Snn_ref,
               sizeof(*(vars->Snn_ref))*nlines);

        /*Calculate the initial Snn_ref correction.*/
        log_mesg("Launching kernel pre_eval_Snn_h at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        int mol_id = vars->LINES->mol;
        pre_eval_Snn_h(nlines,
                       mol_id,
                       vars->LINES->iso,
                       vars->LINES->Vnn,
                       vars->LINES->En,
                       vars->Snn_ref);

        /*Calculate the integrated average layer partial pressure.*/
        vars->x = &((input_data->x[mol])[offset]);
        log_mesg("Launching kernel Curtis_Godson_PsNs_h at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        Curtis_Godson_PsNs_h(nlayers,
                             vars->P,
                             vars->x,
                             vars->N,
                             vars->Psavg,
                             vars->Ns);

        /*Calcluate the lorentz half-width at half-max (HWHM).*/
        log_mesg("Launching kernel eval_gamma_h at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        eval_gamma_h(nlayers,
                     nlines,
                     vars->Pavg,
                     vars->Tavg,
                     vars->Psavg,
                     vars->LINES->Yself,
                     vars->LINES->Yair,
                     vars->LINES->n,
                     vars->GAMMA);

        /*Calcluate the shift in the line center frequency due to the
          pressure.*/
        log_mesg("Launching kernel eval_pShift_h at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        eval_pShift_h(nlayers,
                      nlines,
                      vars->Pavg,
                      vars->LINES->Vnn,
                      vars->LINES->d,
                      vars->PSHIFT);

        /*Calculate the remainder of the Snn_ref correction.*/
        log_mesg("Launching kernel eval_Snn_correction_h at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        eval_Snn_correction_h(nlayers,
                              nlines,
                              mol_id,
                              vars->Tavg,
                              vars->LINES->iso,
                              vars->LINES->Vnn,
                              vars->LINES->En,
                              vars->Snn_ref,
                              vars->S);

        /*Calculate the molecule's optical depths and add them to existing
          values.*/
        log_mesg("Launching kernel eval_profile_h at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        eval_profile_h(mol_id,
                       nlines,
                       nws,
                       w,
                       res,
                       nlayers,
                       breadth,
                       vars->Tavg,
                       vars->GAMMA,
                       vars->PSHIFT,
                       vars->S,
                       vars->Ns,
                       vars->tau_gas);

        if (h2o_ctm && mol_id == H2O)
        {
            /*Calculate the water vapor continuum optical depths.*/
            log_mesg("Launching kernel calc_water_vapor_ctm_optdetph_h"
                         " at point (%d,%d,%d).",
                     time,
                     lon,
                     lat);
            calc_water_vapor_ctm_optdepth_h(nws,
                                            nlayers,
                                            vars->tau_gas,
                                            h2o_continuum->coefs[MTCKD25_S296],
                                            vars->Tavg,
                                            vars->Psavg,
                                            vars->Ns,
                                            h2o_continuum->coefs[CKDS],
                                            h2o_continuum->coefs[MTCKD25_F296],
                                            vars->Pavg,
                                            h2o_continuum->coefs[CKDF]);
        }
        else if (o3_ctm && mol_id == O3)
        {
            /*Calculate the ozone continuum optical depths.*/
            log_mesg("Launching kernel calc_ozone_ctm_optdetph_h at point"
                         " (%d,%d,%d).",
                     time,
                     lon,
                     lat);
            calc_ozone_ctm_optdepth_h(nws,
                                      nlayers,
                                      o3_continuum->cross_section,
                                      vars->Ns,
                                      vars->tau_gas);
        }
    }

    /*Calculate the longwave fluxes.*/
    log_mesg("Launching kernel calc_lw_flux_h at point (%d,%d,%d)",
             time,
             lon,
             lat);
    calc_lw_flux_h(nws,
                   nlayers,
                   vars->lw_flux_down_per_w,
                   vars->lw_flux_up_per_w,
                   vars->Tavg,
                   TSURF,
                   vars->tau_gas,
                   w,
                   res,
                   EMIS,
                   vars->T);

    /*Integrate the fluxes over wavenumber.*/
    log_mesg("Integrating longwave fluxes across wavenumbers at"
                 " point (%d,%d,%d).",
             time,
             lon,
             lat);
    int i;
    for (i=0;i<nlevels;++i)
    {
        fp_t const M_TO_CM = 100.; /*[cm/m]*/
        reimann_sum(&(vars->lw_flux_up_per_w[i*nws]),
                    nws,
                    res*M_TO_CM,
                    &(vars->lw_flux_up[i]));
        reimann_sum(&(vars->lw_flux_down_per_w[i*nws]),
                    nws,
                    res*M_TO_CM,
                    &(vars->lw_flux_down[i]));
    }

    if (MU_DIR >= 0.)
    {
        /*Calculate the shortwave fluxes.*/
        log_mesg("Launching kernel calc_sw_flux_h at point (%d,%d,%d).",
                 time,
                 lon,
                 lat);
        check(calc_sw_flux_h(nlevels,
                             nws,
                             w,
                             res,
                             vars->N,
                             MU_DIR,
                             MU_DIF,
                             vars->tau_gas,
                             SFC_DIR_ALB,
                             SFC_DIF_ALB,
                             solar_flux->incident_sw_flux,
                             SOL_FLUX_RATIO,
                             vars->sw_flux_up_per_w,
                             vars->sw_flux_down_per_w,
                             vars->tau_scatter));

        /*Integrate the fluxes over wavenumber.*/
        log_mesg("Integrating shortwave fluxes across wavenumbers"
                     " at point (%d,%d,%d).",
                 time,
                 lon,
                 lat);
        for (i=0;i<nlevels;++i)
        {
            reimann_sum(&(vars->sw_flux_up_per_w[i*nws]),
                        nws,
                        res,
                        &(vars->sw_flux_up[i]));
            reimann_sum(&(vars->sw_flux_down_per_w[i*nws]),
                        nws,
                        res,
                        &(vars->sw_flux_down[i]));
        }
    }
    return SUCCESS;
}
