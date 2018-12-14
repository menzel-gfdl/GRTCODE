#include <string.h>
#include "debug.h"
#include "floating_point_type.h"
#ifdef __NVCC__
#include "kernels.cuh"
#endif
#include "kernels.h"
#include "launch.h"
#include "molecular_lines.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "spectral_bin.h"
#include "water_vapor_continuum.h"


int launch(GrtContext_t * const context,
           fp_t *p,
           fp_t *t,
           fp_t * const tau
          )
{
    not_null(context);
    not_null(p);
    not_null(t);
    not_null(tau);

    /*Set pointers to input data.*/
    if (context->gpu_id == HOST_ONLY)
    {
        context->p = p;
        context->t = t;
        context->tau = tau;
    }
    else
    {
        gmemcpy(context->p,p,context->num_levels,context->gpu_id,FROM_HOST);
        gmemcpy(context->t,t,context->num_levels,context->gpu_id,FROM_HOST);
    }

    /*Zero out buffers used to accumulate results.*/
    gmemset(context->tau,0,context->num_layers*context->num_wpoints,context->gpu_id);
    gmemset(context->bins.tau,0,context->bins.isize,context->bins.gpu_id);

    /*Calculate the total number density of air molecules integrated across
      each layer.*/
    log_info("Integrating total number density across %d layers.",
             context->num_layers);
    glaunch(calc_number_densities,
            context->num_layers,
            context->gpu_id,
            context->num_layers,
            context->p,
            context->n);

    /*Calculate integrated average layer quantities.*/
    log_info("Calculating Curtis-Godson pressures and temperatures across"
                 " %d layers.",
             context->num_layers);
    glaunch(calc_pressures_and_temperatures,
            context->num_layers,
            context->gpu_id,
            context->num_layers,
            context->p,
            context->t,
            context->pavg,
            context->tavg);

    /*Loop over the molecules and calculate the optical depths.*/
    int m;
    for (m=0;m<context->num_molecules;++m)
    {
        Molecule_t *mol = &(context->mols[m]);
        int index;
        throw(molecule_hash(mol->id,
                            &index));

        /*Calculate the integrated average layer partial pressure.*/
        fp_t const *xp = &(context->x[index*context->num_levels]);
        log_info("Calculating Curtis-Godson partial pressures and abundances"
                     " across %d layers for molecule %s.",
                 context->num_layers,
                 mol->name);
        glaunch(calc_partial_pressures_and_number_densities,
                context->num_layers,
                context->gpu_id,
                context->num_layers,
                context->p,
                xp,
                context->n,
                context->psavg,
                context->ns);

        /*Calculate pressure shifted line center positions.*/
        log_info("Calculating pressure-shifted line center positions"
                     " across %d layers for molecule %s.",
                 context->num_layers,
                 mol->name);
        glaunch(calc_line_centers,
                mol->line_params.num_lines,
                context->gpu_id,
                mol->line_params.num_lines,
                context->num_layers,
                mol->line_params.vnn,
                mol->line_params.d,
                context->pavg,
                context->linecenter);

        /*Calculate total partition functions.*/
        log_info("Calculating total partition functions across %d layers"
                     " for molecule %s.",
                 context->num_layers,
                 mol->name);
        glaunch(calc_partition_functions,
                mol->num_isotopologues,
                context->gpu_id,
                context->num_layers,
                mol->id,
                mol->num_isotopologues,
                context->tavg,
                mol->q);

        /*Calculate temperature-corrected line strengths.*/
        log_info("Calculating temperature-corrected line strengths"
                     " across %d layers for molecule %s.",
                 context->num_layers,
                 mol->name);
        glaunch(calc_line_strengths,
                mol->line_params.num_lines,
                context->gpu_id,
                mol->line_params.num_lines,
                context->num_layers,
                mol->num_isotopologues,
                mol->line_params.iso,
                mol->line_params.snn,
                mol->line_params.vnn,
                mol->line_params.en,
                context->tavg,
                mol->q,
                context->snn);

        /*Calcluate temperature and pressure corrected lorentz half-widths.*/
        log_info("Calculating temperature- and pressure-corrected lorentz"
                     " half-widths across %d layers for molecule %s.",
                 context->num_layers,
                 mol->name);
        glaunch(calc_lorentz_hw,
                mol->line_params.num_lines,
                context->gpu_id,
                mol->line_params.num_lines,
                context->num_layers,
                mol->line_params.n,
                mol->line_params.yair,
                mol->line_params.yself,
                context->tavg,
                context->pavg,
                context->psavg,
                context->gamma);

        /*Calculate doppler half-widths.*/
        log_info("Calculating doppler half-widths across %d layers for"
                     " molecule %s.",
                 context->num_layers,
                 mol->name);
        glaunch(calc_doppler_hw,
                mol->line_params.num_lines,
                context->gpu_id,
                mol->line_params.num_lines,
                context->num_layers,
                mol->mass,
                context->linecenter,
                context->tavg,
                context->alpha);

        /*Calculate the molecule's optical depths and add them to existing
          values.*/
        log_info("Calculating optical depths across %d layers for molecule"
                     " %s.",
                 context->num_layers,
                 mol->name);
        switch (context->optical_depth_method)
        {
            case wavenumber_sweep:
                glaunch(sort_lines,
                        context->num_layers,
                        context->gpu_id,
                        mol->line_params.num_lines,
                        context->num_layers,
                        context->linecenter,
                        context->snn,
                        context->gamma,
                        context->alpha);
                glaunch(calc_optical_depth_bin_sweep,
                        context->bins.n,
                        context->gpu_id,
                        mol->line_params.num_lines,
                        context->num_layers,
                        context->linecenter,
                        context->snn,
                        context->gamma,
                        context->alpha,
                        context->ns,
                        context->bins,
                        context->tau);
                break;
            case line_sweep:
                glaunch(calc_optical_depth_line_sweep,
                        context->bins.n,
                        context->gpu_id,
                        mol->line_params.num_lines,
                        context->num_layers,
                        context->linecenter,
                        context->snn,
                        context->gamma,
                        context->alpha,
                        context->ns,
                        context->bins,
                        context->tau);
                break;
            case line_sample:
                glaunch(calc_optical_depth_line_sample,
                        mol->line_params.num_lines,
                        context->gpu_id,
                        mol->line_params.num_lines,
                        context->num_layers,
                        context->linecenter,
                        context->snn,
                        context->gamma,
                        context->alpha,
                        context->ns,
                        context->bins,
                        context->tau);
                break;
        }

        if (context->use_h2o_ctm && mol->id == H2O)
        {
            /*Calculate the water vapor continuum optical depths.*/
            log_info("Calculating optical depth contribution due to the"
                         " water vapor continuum across %d layers.",
                     context->num_layers);
            glaunch(calc_water_vapor_ctm_optical_depth,
                    context->bins.num_wpoints,
                    context->gpu_id,
                    context->bins.num_wpoints,
                    context->num_layers,
                    context->tau,
                    context->h2o_cc.coefs[MTCKD25_S296],
                    context->tavg,
                    context->psavg,
                    context->ns,
                    context->h2o_cc.coefs[CKDS],
                    context->h2o_cc.coefs[MTCKD25_F296],
                    context->pavg,
                    context->h2o_cc.coefs[CKDF]);
        }
        else if (context->use_o3_ctm && mol->id == O3)
        {
            /*Calculate the ozone continuum optical depths.*/
            log_info("Calculating optical depth contribution due to the"
                         " ozone continuum across %d layers.",
                     context->num_layers);
            glaunch(calc_ozone_ctm_optical_depth,
                    context->bins.num_wpoints,
                    context->gpu_id,
                    context->bins.num_wpoints,
                    context->num_layers,
                    context->o3_cc.cross_section,
                    context->ns,
                    context->tau);
        }
    }

    if (context->optical_depth_method != line_sample)
    {
        /*Interpolate line wing optical depth contributions.*/
        log_info("Interpolating line wing optical depth contributions across"
                     " %d layers.",
                 context->num_layers);
        glaunch(interpolate,
                context->bins.n-1,
                context->gpu_id,
                context->bins,
                context->tau);
        glaunch(interpolate_last_bin,
                context->num_layers,
                context->gpu_id,
                context->bins,
                context->tau);
    }

    if (context->gpu_id != HOST_ONLY)
    {
        gmemcpy(tau,context->tau,context->num_layers*context->num_wpoints,context->gpu_id,FROM_DEVICE);
    }
    return SUCCESS;
}
