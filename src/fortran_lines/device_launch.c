#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef __NVCC__
#include "cuda_helpers.cuh"
#else
#error
#endif
#include "debug.h"
#include "device_launch.h"
#include "eval_gamma.h"
#include "eval_profile.h"
#include "eval_pshift.h"
#include "eval_snn_correction.h"
#include "floating_point_type.h"
#include "integrate_layer.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "parse_HITRAN_file.h"
#include "pre_eval_snn.h"
#include "utils.h"
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
           int const num_molecules,
           LineParams_t ** const line_params,
           double const w0,
           double const wres,
           uint64_t const num_wpoints,
           double const wcutoff,
           int const use_h2o_ctm,
           WaterVaporContinuumCoefs_t * const h2o_cc,
           int const use_o3_ctm,
           OzoneContinuumCoefs_t const * const o3_cc,
           fp_t * const tau)
{
    /*Zero out buffers used to accumulate results.*/
    int num_layers = num_levels - 1;
    HANDLE_ERROR(cudaMemset(tau,
                            0,
                            sizeof(*tau)*num_layers*num_wpoints));

    /*Calculate the total number density of air moleucles integrated across
      each layer.*/
    log_mesg("Integration total number density across %d layers.",
             num_layers);
    integrated_N<<<1,num_layers,0,0>>>(num_layers,
                                       P,
                                       N);

    /*Calculate integrated average layer quantities.*/
    log_mesg("Calculating Curtis-Godson pressure and temperature across"
                 " %d layers.",
             num_layers);
    Curtis_Godson_PT<<<1,num_layers,0,0>>>(num_layers,
                                           P,
                                           T,
                                           Pavg,
                                           Tavg);

    /*Loop over the molecules and calculate the optical depths.*/
    int min_grid_size;
    int dim_block;
    int dim_grid;
    int mol;
    for (mol=0;mol<num_molecules;++mol)
    {
        /*Copy line parameters to device.*/
        unsigned int num_lines = line_params[mol]->num_lines;
        int mol_id = line_params[mol]->mol;
        HANDLE_ERROR(cudaMemcpy(lines->iso,
                                line_params[mol]->iso,
                                sizeof(*(lines->iso))*num_lines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(lines->vnn,
                                line_params[mol]->vnn,
                                sizeof(*(lines->vnn))*num_lines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(lines->snn_ref,
                                line_params[mol]->snn_ref,
                                sizeof(*(lines->snn_ref))*num_lines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(lines->yair,
                                line_params[mol]->yair,
                                sizeof(*(lines->yair))*num_lines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(lines->yself,
                                line_params[mol]->yself,
                                sizeof(*(lines->yself))*num_lines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(lines->en,
                                line_params[mol]->en,
                                sizeof(*(lines->en))*num_lines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(lines->n,
                                line_params[mol]->n,
                                sizeof(*(lines->n))*num_lines,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(lines->d,
                                line_params[mol]->d,
                                sizeof(*(lines->d))*num_lines,
                                cudaMemcpyHostToDevice));

        /*Calculate the initial Snn_ref correction.*/
        log_mesg("Launching kernel pre_eval_snn across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        pre_eval_snn,
                                                        0,
                                                        ((int)num_lines)));
        dim_grid = (((int)num_lines) + dim_block - 1)/dim_block;
        pre_eval_snn<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(num_lines,
                                                                                 mol_id,
                                                                                 lines->iso,
                                                                                 lines->vnn,
                                                                                 lines->en,
                                                                                 lines->snn_ref);

        /*Calculate the integrated average layer partial pressure.*/
        fp_t const *xp = &(x[mol*num_levels]);
        log_mesg("Calculating Curtis-Godson partial pressure and abundance"
                     " across %d layers for molecule %d.",
                 num_layers,
                 mol);
        Curtis_Godson_PsNs<<<1,num_layers,0,0>>>(num_layers,
                                                 P,
                                                 xp,
                                                 N,
                                                 Psavg,
                                                 Ns);

        /*Calcluate the lorentz half-width at half-max (HWHM).*/
        log_mesg("Launching kernel eval_gamma across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        eval_gamma,
                                                        0,
                                                        ((int)num_lines)));
        dim_grid = (((int)num_lines) + dim_block - 1)/dim_block;
        eval_gamma<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(num_layers,
                                                                               num_lines,
                                                                               Pavg,
                                                                               Tavg,
                                                                               Psavg,
                                                                               lines->yself,
                                                                               lines->yair,
                                                                               lines->n,
                                                                               gamma);

        /*Calcluate the shift in the line center frequency due to the
          pressure.*/
        log_mesg("Launching kernel eval_pshift across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        eval_pshift,
                                                        0,
                                                        ((int)num_lines)));
        dim_grid = (((int)num_lines) + dim_block - 1)/dim_block;
        eval_pshift<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(num_layers,
                                                                                num_lines,
                                                                                Pavg,
                                                                                lines->vnn,
                                                                                lines->d,
                                                                                Pshift);

        /*Calculate the remainder of the Snn_ref correction.*/
        log_mesg("Launching kernel eval_snn_correction across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        eval_snn_correction,
                                                        0,
                                                        ((int)num_lines)));
        dim_grid = (((int)num_lines) + dim_block - 1)/dim_block;
        eval_snn_correction<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(num_layers,
                                                                                        num_lines,
                                                                                        mol_id,
                                                                                        Tavg,
                                                                                        lines->iso,
                                                                                        lines->vnn,
                                                                                        lines->en,
                                                                                        lines->snn_ref,
                                                                                        s);

        /*Calculate the molecule's optical depths and add them to existing
          values.*/
        log_mesg("Launching kernel eval_profile across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                        &dim_block,
                                                        eval_profile,
                                                        0,
                                                        ((int)num_lines)));
        dim_grid = (((int)num_lines) + dim_block - 1)/dim_block;
        eval_profile<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(mol_id,
                                                                                 num_lines,
                                                                                 num_wpoints,
                                                                                 w0,
                                                                                 wres,
                                                                                 num_layers,
                                                                                 wcutoff,
                                                                                 Tavg,
                                                                                 gamma,
                                                                                 Pshift,
                                                                                 s,
                                                                                 Ns,
                                                                                 tau);

        if (use_h2o_ctm && mol_id == H2O)
        {
            /*Calculate the water vapor continuum optical depths.*/
            log_mesg("Calculating optical depth due to the water vapor"
                         " continuum across %d layers.",
                     num_layers);
            HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                            &dim_block,
                                                            calc_water_vapor_ctm_optical_depth,
                                                            0,
                                                            ((int)num_lines)));
            dim_grid = (((int)num_lines) + dim_block - 1)/dim_block;
            calc_water_vapor_ctm_optical_depth<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(num_wpoints,
                                                                                                           num_layers,
                                                                                                           tau,
                                                                                                           h2o_cc->coefs[MTCKD25_S296],
                                                                                                           Tavg,
                                                                                                           Psavg,
                                                                                                           Ns,
                                                                                                           h2o_cc->coefs[CKDS],
                                                                                                           h2o_cc->coefs[MTCKD25_F296],
                                                                                                           Pavg,
                                                                                                           h2o_cc->coefs[CKDF]);
        }
        else if (use_o3_ctm && mol_id == O3)
        {
            /*Calculate the ozone continuum optical depths.*/
            log_mesg("Calculating optical depth due to the ozone"
                         " continuum across %d layers.",
                     num_layers);
            HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&min_grid_size,
                                                            &dim_block,
                                                            calc_ozone_ctm_optical_depth,
                                                            0,
                                                            ((int)num_lines)));
            dim_grid = (((int)num_lines) + dim_block - 1)/dim_block;
            calc_ozone_ctm_optical_depth<<<((unsigned int)dim_grid),((unsigned int)dim_block),0,0>>>(num_wpoints,
                                                                                                     num_layers,
                                                                                                     o3_cc->cross_section,
                                                                                                     Ns,
                                                                                                     tau);
        }
    }
    return SUCCESS;
}
