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
#include "integrate_layer.h"
#include "launch.h"
#include "lw_flux.h"
#include "model_fields.h"
#include "parseHITRANfile.h"
#include "pre_eval_Snn.h"
#include "utils.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif


int alloc_work_vars(WorkVars_t *vars,
                    int const nlevels,
                    int const nlines,
                    int const nws,
                    int const launch_type)
{
    not_null(vars);
    check(check_launch_mode(launch_type));
    int const nlayers = nlevels - 1;
    int const l = nlayers*nlines;
    int const m = nlevels*nws;
    if (launch_type == HOST_LAUNCH)
    {
        vars->P = NULL;
        vars->T = NULL;
        vars->TSURF = NULL;
        vars->EMIS = NULL;
        vars->x = NULL;
        malloc_ptr(vars->Pavg,nlayers);
        malloc_ptr(vars->Tavg,nlayers);
        malloc_ptr(vars->N,nlayers);
        malloc_ptr(vars->Psavg,nlayers);
        malloc_ptr(vars->GAMMA,l);
        malloc_ptr(vars->PSHIFT,l);
        malloc_ptr(vars->S,l);
        vars->LINES = NULL;
        malloc_ptr(vars->Snn_ref,nlines);
        vars->tau = NULL;
        malloc_ptr(vars->lw_flux_down_per_w,m);
        malloc_ptr(vars->lw_flux_up_per_w,m);
        vars->lw_flux_down = NULL;
        vars->lw_flux_up = NULL;
    }
    else
    {
        using_gpu();
#ifdef __NVCC__
        HANDLE_ERROR(cudaMalloc(&(vars->P),
                                sizeof(*(vars->P))*nlevels));
        HANDLE_ERROR(cudaMalloc(&(vars->T),
                                sizeof(*(vars->T))*nlevels));
        HANDLE_ERROR(cudaMalloc(&(vars->TSURF),
                                sizeof(*(vars->TSURF))));
        HANDLE_ERROR(cudaMalloc(&(vars->EMIS),
                                sizeof(*(vars->EMIS))));
        HANDLE_ERROR(cudaMalloc(&(vars->x),
                                sizeof(*(vars->x))*nlevels));
        HANDLE_ERROR(cudaMalloc(&(vars->Pavg),
                                sizeof(*(vars->Pavg))*nlayers));
        HANDLE_ERROR(cudaMalloc(&(vars->Tavg),
                                sizeof(*(vars->Tavg))*nlayers));
        HANDLE_ERROR(cudaMalloc(&(vars->N),
                                sizeof(*(vars->N))*nlayers));
        HANDLE_ERROR(cudaMalloc(&(vars->Psavg),
                                sizeof(*(vars->Psavg))*nlayers));
        HANDLE_ERROR(cudaMalloc(&(vars->GAMMA),
                                sizeof(*(vars->GAMMA))*l));
        HANDLE_ERROR(cudaMalloc(&(vars->PSHIFT),
                                sizeof(*(vars->PSHIFT))*l));
        HANDLE_ERROR(cudaMalloc(&(vars->S),
                                sizeof(*(vars->S))*l));
        check(alloc_line_params_device(&(vars->LINES),
                                       nlines));
        int const n = nlayers*nws;
        HANDLE_ERROR(cudaMalloc(&(vars->tau),
                                sizeof(*(vars->tau))*n));
        HANDLE_ERROR(cudaMalloc(&(vars->lw_flux_down_per_w),
                                sizeof(*(vars->lw_flux_down_per_w))*m));
        HANDLE_ERROR(cudaMalloc(&(vars->lw_flux_up_per_w),
                                sizeof(*(vars->lw_flux_up_per_w))*m));
        HANDLE_ERROR(cudaMalloc(&(vars->lw_flux_down),
                                sizeof(*(vars->lw_flux_down))*nlevels));
        HANDLE_ERROR(cudaMalloc(&(vars->lw_flux_up),
                                sizeof(*(vars->lw_flux_up))*nlevels));
#endif
    }
    return SUCCESS;
}


int free_work_vars(WorkVars_t *vars,
                   int const launch_type)
{
    not_null(vars);
    check(check_launch_mode(launch_type));
    if (launch_type == HOST_LAUNCH)
    {
        free(vars->Pavg);
        free(vars->Tavg);
        free(vars->N);
        free(vars->Psavg);
        free(vars->GAMMA);
        free(vars->PSHIFT);
        free(vars->S);
        free(vars->Snn_ref);
        free(vars->lw_flux_down_per_w);
        free(vars->lw_flux_up_per_w);
    }
    else
    {
        using_gpu();
#ifdef __NVCC__
        HANDLE_ERROR(cudaFree(&(vars->P)));
        HANDLE_ERROR(cudaFree(&(vars->T)));
        HANDLE_ERROR(cudaFree(&(vars->TSURF)));
        HANDLE_ERROR(cudaFree(&(vars->EMIS)));
        HANDLE_ERROR(cudaFree(&(vars->x)));
        HANDLE_ERROR(cudaFree(&(vars->Pavg)));
        HANDLE_ERROR(cudaFree(&(vars->Tavg)));
        HANDLE_ERROR(cudaFree(&(vars->N)));
        HANDLE_ERROR(cudaFree(&(vars->Psavg)));
        HANDLE_ERROR(cudaFree(&(vars->GAMMA)));
        HANDLE_ERROR(cudaFree(&(vars->PSHIFT)));
        HANDLE_ERROR(cudaFree(&(vars->S)));
        check(free_line_params_device(&(vars->LINES)));
        HANDLE_ERROR(cudaFree(&(vars->tau)));
        HANDLE_ERROR(cudaFree(&(vars->lw_flux_down_per_w)));
        HANDLE_ERROR(cudaFree(&(vars->lw_flux_up_per_w)));
        HANDLE_ERROR(cudaFree(&(vars->lw_flux_down)));
        HANDLE_ERROR(cudaFree(&(vars->lw_flux_up)));
#endif
    }
    return SUCCESS;
}


int launch_host(WorkVars_t * const vars,
                req_model_fields_t * const input_data,
                int const time,
                int const lon,
                int const lat,
                int const nmols,
                line_params_t ** const line_params,
                unsigned int const nws,
                fp_t const w,
                double const res,
                int const breadth,
                int const continuum,
                OutputFields_t * const output_data)
{
    not_null(vars);
    not_null(input_data);
    not_null(line_params);
    not_null(output_data);

    /*Point to output buffers.*/
    vars->tau = output_data->tau;
    vars->lw_flux_down = output_data->lw_flux_down;
    vars->lw_flux_up = output_data->lw_flux_up;

    /*Zero out buffers used to accumulate results.*/
    int nlayers = input_data->nlevel - 1;
    memset(vars->tau,
           0,
           sizeof(*(vars->tau))*nlayers*nws);
    memset(vars->lw_flux_down,
           0,
           sizeof(*(vars->lw_flux_down))*(input_data->nlevel));
    memset(vars->lw_flux_up,
           0,
           sizeof(*(vars->lw_flux_up))*(input_data->nlevel));

    /*Point to the correct column of input data.*/
    unsigned int offset = time*input_data->nlon*input_data->nlat +
                          lon*input_data->nlat + lat;
    vars->TSURF = &(input_data->TSURF[offset]);
    vars->EMIS = &(input_data->EMIS[offset]);
    offset *= input_data->nlevel;
    vars->P = &(input_data->P[offset]);
    vars->T = &(input_data->T[offset]);

    /*Calculate integrated average layer quantities.*/
    log_mesg("Launching kernel get_avg_NTP_h at point (%d,%d,%d).",
             time,
             lon,
             lat);
    get_avg_TP_h(nlayers,
                 vars->P,
                 vars->T,
                 vars->Pavg,
                 vars->Tavg);

    /*Loop over the molecules and calculate the optical depths.*/
    int mol;
    for (mol=0;mol<nmols;++mol)
    {
        /*Point to the line parameters for the current molecule.*/
        vars->LINES = line_params[mol];

        /*Copy the Snn_ref array.  This only needs to be done by the host.*/
        memcpy(vars->Snn_ref,
               line_params[mol]->Snn_ref,
               sizeof(*(vars->Snn_ref))*(line_params[mol]->nLines));

        /*Calculate the initial Snn_ref correction.*/
        log_mesg("Launching kernel pre_eval_Snn_h at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        pre_eval_Snn_h(vars->LINES->nLines,
                       vars->LINES->mol,
                       vars->LINES->iso,
                       vars->LINES->Vnn,
                       vars->LINES->En,
                       vars->Snn_ref);

        /*Calculate the integrated average layer partial pressure.*/
        vars->x = &((input_data->x[mol])[offset]);
        log_mesg("Launching kernel get_avg_Ps_h at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        get_avg_NPs_h(nlayers,
                      vars->x,
                      vars->P,
                      vars->N,
                      vars->Psavg);

        /*Calcluate the lorentz half-width at half-max (HWHM).*/
        log_mesg("Launching kernel eval_gamma_h at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        eval_gamma_h(nlayers,
                     vars->LINES->nLines,
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
                      vars->LINES->nLines,
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
                              vars->LINES->nLines,
                              vars->LINES->mol,
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

        eval_profile_h(vars->LINES->mol,
                       vars->LINES->nLines,
                       nws,
                       w,
                       res,
                       nlayers,
                       breadth,
                       vars->Tavg,
                       vars->GAMMA,
                       vars->PSHIFT,
                       vars->S,
                       vars->N,
                       vars->tau);
    }

    if (continuum)
    {
        /*Calculate the water vapor continuum optical depths.*/
/*
        log_mesg("Launching kernel calc_ctm_optdetph_h at point (%d,%d,%d).",
                 time,
                 lon,
                 lat);
        calc_ctm_optdepth_h(nF,
                            numLayers,
                            out,
                            CS_h,
                            T,
                            &(PS[H2O*atmosData->npfull]),
                            DELTAZ,
                            T0_h,
                            CF_h,
                            P,
                            T0F_h);
*/
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
                   *(vars->TSURF),
                   vars->tau,
                   w,
                   res,
                   *(vars->EMIS),
                   vars->T);

    /*Integrate the fluxes over wavenumber.*/
    log_mesg("Integrating downward fluxes across wavenumbers at point"
                 " (%d,%d,%d).",
             time,
             lon,
             lat);
    integrate_fluxes(nws,
                     input_data->nlevel,
                     vars->lw_flux_down_per_w,
                     vars->lw_flux_down,
                     res);
    log_mesg("Integrating upward fluxes across wavenumbers at point"
                 " (%d,%d,%d).",
             time,
             lon,
             lat);
    integrate_fluxes(nws,
                     input_data->nlevel,
                     vars->lw_flux_up_per_w,
                     vars->lw_flux_up,
                     res);
    return SUCCESS;
}
