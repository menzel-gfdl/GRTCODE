#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "continuum.h"
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
#include "molecules.h"
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
                    int const put_on_device)
{
    not_null(vars);
    int const nlayers = nlevels - 1;
    int const l = nlayers*nlines;
    int const m = nlevels*nws;
    if (put_on_device)
    {
        using_gpu();
#ifdef __NVCC__
        HANDLE_ERROR(cudaMalloc(&(vars->P),
                                sizeof(*(vars->P))*nlevels));
        HANDLE_ERROR(cudaMalloc(&(vars->T),
                                sizeof(*(vars->T))*nlevels));
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
        vars->LINES = NULL;
        check(alloc_line_params_device(&(vars->LINES),
                                       nlines));
        int const n = nlayers*nws;
        HANDLE_ERROR(cudaMalloc(&(vars->tau),
                                sizeof(*(vars->tau))*n));
        HANDLE_ERROR(cudaMalloc(&(vars->lw_flux_down_per_w),
                                sizeof(*(vars->lw_flux_down_per_w))*m));
        HANDLE_ERROR(cudaMalloc(&(vars->lw_flux_up_per_w),
                                sizeof(*(vars->lw_flux_up_per_w))*m));
#endif
    }
    else
    {
        vars->P = NULL;
        vars->T = NULL;
        vars->x = NULL;
        check(malloc_ptr((void **)(&(vars->Pavg)),
                         sizeof(*(vars->Pavg))*nlayers));
        check(malloc_ptr((void **)(&(vars->Tavg)),
                         sizeof(*(vars->Tavg))*nlayers));
        check(malloc_ptr((void **)(&(vars->N)),
                         sizeof(*(vars->N))*nlayers));
        check(malloc_ptr((void **)(&(vars->Psavg)),
                         sizeof(*(vars->Psavg))*nlayers));
        check(malloc_ptr((void **)(&(vars->GAMMA)),
                         sizeof(*(vars->GAMMA))*l));
        check(malloc_ptr((void **)(&(vars->PSHIFT)),
                         sizeof(*(vars->PSHIFT))*l));
        check(malloc_ptr((void **)(&(vars->S)),
                         sizeof(*(vars->S))*l));
        vars->LINES = NULL;
        check(malloc_ptr((void **)(&(vars->Snn_ref)),
                         sizeof(*(vars->Snn_ref))*nlines));
        vars->tau = NULL;
        check(malloc_ptr((void **)(&(vars->lw_flux_down_per_w)),
                         sizeof(*(vars->lw_flux_down_per_w))*m));
        check(malloc_ptr((void **)(&(vars->lw_flux_up_per_w)),
                         sizeof(*(vars->lw_flux_up_per_w))*m));
        vars->lw_flux_down = NULL;
        vars->lw_flux_up = NULL;
    }
    return SUCCESS;
}


int free_work_vars(WorkVars_t *vars,
                   int const on_device)
{
    not_null(vars);
    if (on_device)
    {
        using_gpu();
#ifdef __NVCC__
        HANDLE_ERROR(cudaFree(vars->P));
        HANDLE_ERROR(cudaFree(vars->T));
        HANDLE_ERROR(cudaFree(vars->x));
        HANDLE_ERROR(cudaFree(vars->Pavg));
        HANDLE_ERROR(cudaFree(vars->Tavg));
        HANDLE_ERROR(cudaFree(vars->N));
        HANDLE_ERROR(cudaFree(vars->Psavg));
        HANDLE_ERROR(cudaFree(vars->GAMMA));
        HANDLE_ERROR(cudaFree(vars->PSHIFT));
        HANDLE_ERROR(cudaFree(vars->S));
        check(free_line_params_device(&(vars->LINES)));
        HANDLE_ERROR(cudaFree(vars->tau));
        HANDLE_ERROR(cudaFree(vars->lw_flux_down_per_w));
        HANDLE_ERROR(cudaFree(vars->lw_flux_up_per_w));
#endif
    }
    else
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
                ContinuumCoefs_t * const h2o_continuum,
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
    int nlevels = input_data->nlevel;
    int nlayers = nlevels - 1;
    memset(vars->tau,
           0,
           sizeof(*(vars->tau))*nlayers*nws);
    memset(vars->lw_flux_down,
           0,
           sizeof(*(vars->lw_flux_down))*nlevels);
    memset(vars->lw_flux_up,
           0,
           sizeof(*(vars->lw_flux_up))*nlevels);

    /*Point to the correct column of input data.*/
    unsigned int offset = time*input_data->nlon*input_data->nlat +
                          lon*input_data->nlat + lat;
    fp_t const TSURF = input_data->TSURF[offset];
    fp_t const EMIS = input_data->EMIS[offset];
    offset *= nlevels;
    vars->P = &(input_data->P[offset]);
    vars->T = &(input_data->T[offset]);

    /*Calculate integrated average layer quantities.*/
    log_mesg("Launching kernel get_avg_TP_h at point (%d,%d,%d).",
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
        log_mesg("Launching kernel get_avg_NPs_h at point (%d,%d,%d)"
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
                       vars->N,
                       vars->tau);

        if (mol_id == H2O && continuum)
        {
            /*Calculate the water vapor continuum optical depths.*/
            log_mesg("Launching kernel calc_ctm_optdetph_h at point"
                         " (%d,%d,%d).",
                     time,
                     lon,
                     lat);
            calc_ctm_optdepth_h(nws,
                                nlayers,
                                vars->tau,
                                h2o_continuum->coefs[CS],
                                vars->Tavg,
                                vars->Psavg,
                                vars->N,
                                h2o_continuum->coefs[T0],
                                h2o_continuum->coefs[CF],
                                vars->Pavg,
                                h2o_continuum->coefs[T0F]);
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
                   vars->tau,
                   w,
                   res,
                   EMIS,
                   vars->T);

    /*Integrate the fluxes over wavenumber.*/
    log_mesg("Integrating downward fluxes across wavenumbers at point"
                 " (%d,%d,%d).",
             time,
             lon,
             lat);
    integrate_fluxes(nws,
                     nlevels,
                     vars->lw_flux_down_per_w,
                     vars->lw_flux_down,
                     res);
    log_mesg("Integrating upward fluxes across wavenumbers at point"
                 " (%d,%d,%d).",
             time,
             lon,
             lat);
    integrate_fluxes(nws,
                     nlevels,
                     vars->lw_flux_up_per_w,
                     vars->lw_flux_up,
                     res);
    return SUCCESS;
}


#ifdef __NVCC__
int launch_device(WorkVars_t * const vars,
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
                  ContinuumCoefs_t * const h2o_continuum,
                  OutputFields_t * const output_data)
{
    not_null(vars);
    not_null(input_data);
    not_null(line_params);
    not_null(output_data);

    /*Zero out buffers used to accumulate results.*/
    int nlevels = input_data->nlevel;
    int nlayers = nlevels - 1;
    HANDLE_ERROR(cudaMemset(vars->tau,
                            0,
                            sizeof(*(vars->tau))*nlayers*nws));
    memset(output_data->lw_flux_down,
           0,
           sizeof(*(output_data->lw_flux_down))*nlevels);
    memset(output_data->lw_flux_up,
           0,
           sizeof(*(output_data->lw_flux_up))*nlevels);

    /*Copy the correct column of input data to the device.*/
    log_mesg("Copying column of input data from host to device at point"
                 " (%d,%d,%d).",
             time,
             lon,
             lat);
    unsigned int offset = time*input_data->nlon*input_data->nlat +
                          lon*input_data->nlat + lat;
    fp_t const TSURF = input_data->TSURF[offset];
    fp_t const EMIS = input_data->EMIS[offset];
    offset *= nlevels;
    HANDLE_ERROR(cudaMemcpy(vars->P,
                            &(input_data->P[offset]),
                            sizeof(*(input_data->P))*nlevels,
                            cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(vars->T,
                            &(input_data->T[offset]),
                            sizeof(*(input_data->T))*nlevels,
                            cudaMemcpyHostToDevice));

    /*Calculate integrated average layer quantities.*/
    log_mesg("Launching kernel get_avg_TP at point (%d,%d,%d).",
             time,
             lon,
             lat);
    get_avg_TP<<<1,nlayers,0,0>>>(nlayers,
                                  vars->P,
                                  vars->T,
                                  vars->Pavg,
                                  vars->Tavg);

    /*Loop over the molecules and calculate the optical depths.*/
    int min_grid_size;
    int dim_block;
    int dim_grid;
    int mol;
    for (mol=0;mol<nmols;++mol)
    {
        /*Copy the molecular abundance from the host to the device.*/
        log_mesg("Copying abundance from host to device for molecule %d.",
                 mol);
        HANDLE_ERROR(cudaMemcpy(vars->x,
                                &((input_data->x[mol])[offset]),
                                sizeof(*(vars->x))*nlevels,
                                cudaMemcpyHostToDevice));

        /*Copy the line parameters for the current molecule to the device.*/
        log_mesg("Copying line parameters from host to device for"
                     " molecule %d.",
                 mol);
        unsigned int nlines = line_params[mol]->nLines;
        HANDLE_ERROR(cudaMemcpy(vars->LINES->iso,
                                line_params[mol]->iso,
                                sizeof(*(vars->LINES->iso))*nlines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(vars->LINES->Vnn,
                                line_params[mol]->Vnn,
                                sizeof(*(vars->LINES->Vnn))*nlines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(vars->LINES->Snn_ref,
                                line_params[mol]->Snn_ref,
                                sizeof(*(vars->LINES->Snn_ref))*nlines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(vars->LINES->Yair,
                                line_params[mol]->Yair,
                                sizeof(*(vars->LINES->Yair))*nlines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(vars->LINES->Yself,
                                line_params[mol]->Yself,
                                sizeof(*(vars->LINES->Yself))*nlines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(vars->LINES->En,
                                line_params[mol]->En,
                                sizeof(*(vars->LINES->En))*nlines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(vars->LINES->n,
                                line_params[mol]->n,
                                sizeof(*(vars->LINES->n))*nlines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(vars->LINES->d,
                                line_params[mol]->d,
                                sizeof(*(vars->LINES->d))*nlines,
                                cudaMemcpyHostToDevice));

        /*Calculate the thread-block size and number of thread blocks that
          maximizes the occupancy on the device.  Round up to make sure
          that all input data is used.  The CUDA API may produce a warning
          that can be safely ignored depending on the sdk version and
          -W flags.*/
        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        pre_eval_Snn,
                                                        0,
                                                        ((int)nlines)));
        dim_grid = (((int)nlines) + dim_block - 1)/dim_block;

        /*Calculate the initial Snn_ref correction.*/
        log_mesg("Launching kernel pre_eval_Snn at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        int mol_id = line_params[mol]->mol;
        pre_eval_Snn<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(nlines,
                                                                                 mol_id,
                                                                                 vars->LINES->iso,
                                                                                 vars->LINES->Vnn,
                                                                                 vars->LINES->En,
                                                                                 vars->LINES->Snn_ref);

        /*Calculate the integrated average layer partial pressure.*/
        log_mesg("Launching kernel get_avg_NPs at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        get_avg_NPs<<<1,nlayers,0,0>>>(nlayers,
                                       vars->x,
                                       vars->P,
                                       vars->N,
                                       vars->Psavg);

        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        eval_gamma,
                                                        0,
                                                        ((int)nlines)));
        dim_grid = (((int)nlines) + dim_block - 1)/dim_block;

        /*Calcluate the lorentz half-width at half-max (HWHM).*/
        log_mesg("Launching kernel eval_gamma at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        eval_gamma<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(nlayers,
                                                                               nlines,
                                                                               vars->Pavg,
                                                                               vars->Tavg,
                                                                               vars->Psavg,
                                                                               vars->LINES->Yself,
                                                                               vars->LINES->Yair,
                                                                               vars->LINES->n,
                                                                               vars->GAMMA);

        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        eval_pShift,
                                                        0,
                                                        ((int)nlines)));
        dim_grid = (((int)nlines) + dim_block - 1)/dim_block;

        /*Calcluate the shift in the line center frequency due to the
          pressure.*/
        log_mesg("Launching kernel eval_pShift at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        eval_pShift<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(nlayers,
                                                                                nlines,
                                                                                vars->Pavg,
                                                                                vars->LINES->Vnn,
                                                                                vars->LINES->d,
                                                                                vars->PSHIFT);

        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        eval_Snn_correction,
                                                        0,
                                                        ((int)nlines)));
        dim_grid = (((int)nlines) + dim_block - 1)/dim_block;

        /*Calculate the remainder of the Snn_ref correction.*/
        log_mesg("Launching kernel eval_Snn_correction at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        eval_Snn_correction<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(nlayers,
                                                                                        nlines,
                                                                                        mol_id,
                                                                                        vars->Tavg,
                                                                                        vars->LINES->iso,
                                                                                        vars->LINES->Vnn,
                                                                                        vars->LINES->En,
                                                                                        vars->LINES->Snn_ref,
                                                                                        vars->S);

        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        eval_profile,
                                                        0,
                                                        ((int)nlines)));
        dim_grid = (((int)nlines) + dim_block - 1)/dim_block;

        /*Calculate the molecule's optical depths and add them to existing
          values.*/
        log_mesg("Launching kernel eval_profile at point (%d,%d,%d)"
                     " for molecule %d.",
                 time,
                 lon,
                 lat,
                 mol);
        eval_profile<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(mol_id,
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
                                                                                 vars->N,
                                                                                 vars->tau);

        if (mol_id == H2O && continuum)
        {
            HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                            &dim_block,
                                                            calc_ctm_optdepth,
                                                            0,
                                                            ((int)nws)));
            dim_grid = (((int)nws) + dim_block - 1)/dim_block;

            /*Calculate the water vapor continuum optical depths.*/
            log_mesg("Launching kernel calc_ctm_optdetph at point"
                         " (%d,%d,%d).",
                     time,
                     lon,
                     lat);
            calc_ctm_optdepth<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(nws,
                                                                                          nlayers,
                                                                                          vars->tau,
                                                                                          h2o_continuum->coefs[CS],
                                                                                          vars->Tavg,
                                                                                          vars->Psavg,
                                                                                          vars->N,
                                                                                          h2o_continuum->coefs[T0],
                                                                                          h2o_continuum->coefs[CF],
                                                                                          vars->Pavg,
                                                                                          h2o_continuum->coefs[T0F]);
        }
    }

    HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                    &dim_block,
                                                    calc_lw_flux,
                                                    0,
                                                    ((int)nws)));
    dim_grid = (((int)nws) + dim_block - 1)/dim_block;

    /*Calculate the longwave fluxes.*/
    log_mesg("Launching kernel calc_lw_flux at point (%d,%d,%d)",
             time,
             lon,
             lat);
    calc_lw_flux<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(nws,
                                                                             nlayers,
                                                                             vars->lw_flux_down_per_w,
                                                                             vars->lw_flux_up_per_w,
                                                                             vars->Tavg,
                                                                             TSURF,
                                                                             vars->tau,
                                                                             w,
                                                                             res,
                                                                             EMIS,
                                                                             vars->T);

    /*Copy the optical depths and fluxes from the device to the host.*/
    log_mesg("Copying optical depths and fluxes from device to host"
                 " at point (%d,%d,%d).",
             time,
             lon,
             lat);
    HANDLE_ERROR(cudaMemcpy(output_data->tau,
                            vars->tau,
                            sizeof(*(vars->tau))*nws*nlayers,
                            cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaMemcpy(output_data->lw_flux_down_per_w,
                            vars->lw_flux_down_per_w,
                            sizeof(*(vars->lw_flux_down_per_w))*nws*nlevels,
                            cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaMemcpy(output_data->lw_flux_up_per_w,
                            vars->lw_flux_up_per_w,
                            sizeof(*(vars->lw_flux_up_per_w))*nws*nlevels,
                            cudaMemcpyDeviceToHost));

    /*Integrate the fluxes over wavenumber.*/
    log_mesg("Integrating downward fluxes across wavenumbers at point"
                 " (%d,%d,%d).",
             time,
             lon,
             lat);
    integrate_fluxes(nws,
                     nlevels,
                     output_data->lw_flux_down_per_w,
                     output_data->lw_flux_down,
                     res);
    log_mesg("Integrating upward fluxes across wavenumbers at point"
                 " (%d,%d,%d).",
             time,
             lon,
             lat);
    integrate_fluxes(nws,
                     nlevels,
                     output_data->lw_flux_up_per_w,
                     output_data->lw_flux_up,
                     res);
    return SUCCESS;
}
#endif
