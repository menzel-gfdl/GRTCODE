#include <string.h>
#include "calc_optical_depth.h"
#include "debug.h"
#include "floating_point_type.h"
#include "host_launch.h"
#include "molecular_lines.h"
#include "molecules.h"
#include "spectral_bin.h"

#ifdef FOO
#include "ozone_continuum.h"
#include "water_vapor_continuum.h"
#endif


int launch_h(GrtContext_t * const context,
             fp_t const * const p,
             fp_t const * const t,
             fp_t * const tau)

{
    not_null(context);
    not_null(p);
    not_null(t);
    not_null(tau);

    /*Zero out buffers used to accumulate results.*/
    memset(tau,
           0,
           sizeof(*tau)*context->num_layers*context->num_wpoints);
    memset(context->bins.tau,
           0,
           sizeof(*(context->bins.tau))*context->bins.isize);

    /*Calculate the total number density of air molecules integrated across
      each layer.*/
    log_info("Integrating total number density across %d layers.",
             context->num_layers);
    check(calc_number_densities(context->num_layers,
                                p,
                                context->n));

    /*Calculate integrated average layer quantities.*/
    log_info("Calculating Curtis-Godson pressures and temperatures across"
                 " %d layers.",
             context->num_layers);
    check(calc_pressures_and_temperatures(context->num_layers,
                                          p,
                                          t,
                                          context->pavg,
                                          context->tavg));

    /*Loop over the molecules and calculate the optical depths.*/
    int m;
    for (m=0;m<context->num_molecules;++m)
    {
        Molecule_t *mol = &(context->mols[m]);
        int index;
        check(molecule_hash(mol->id,
                            &index));

        /*Calculate the integrated average layer partial pressure.*/
        fp_t const *xp = &(context->x[index*context->num_levels]);
        log_info("Calculating Curtis-Godson partial pressures and abundances"
                     " across %d layers for molecule %s.",
                 context->num_layers,
                 mol->name);
        check(calc_partial_pressures_and_number_densities(context->num_layers,
                                                          p,
                                                          xp,
                                                          context->n,
                                                          context->psavg,
                                                          context->ns));

        /*Calculate pressure shifted line center positions.*/
        log_info("Calculating pressure-shifted line center positions"
                     " across %d layers for molecule %s.",
                 context->num_layers,
                 mol->name);
        check(calc_line_centers(mol->line_params.num_lines,
                                context->num_layers,
                                mol->line_params.vnn,
                                mol->line_params.d,
                                context->pavg,
                                context->linecenter));

        /*Calculate temperature-corrected line strengths.*/
        log_info("Calculating temperature-corrected line strengths"
                     " across %d layers for molecule %s.",
                 context->num_layers,
                 mol->name);
        check(calc_line_strengths(mol->line_params.num_lines,
                                  context->num_layers,
                                  mol->id,
                                  mol->line_params.iso,
                                  mol->line_params.snn,
                                  mol->line_params.vnn,
                                  mol->line_params.en,
                                  context->tavg,
                                  context->snn));

        /*Calcluate temperature and pressure corrected lorentz half-widths.*/
        log_info("Calculating temperature- and pressure-corrected lorentz"
                     " half-widths across %d layers for molecule %s.",
                 context->num_layers,
                 mol->name);
        check(calc_lorentz_hw(mol->line_params.num_lines,
                              context->num_layers,
                              mol->line_params.n,
                              mol->line_params.yair,
                              mol->line_params.yself,
                              context->tavg,
                              context->pavg,
                              context->psavg,
                              context->gamma));

        /*Calculate doppler half-widths.*/
        log_info("Calculating doppler half-widths across %d layers for"
                     " molecule %s.",
                 context->num_layers,
                 mol->name);
        check(calc_doppler_hw(mol->line_params.num_lines,
                              context->num_layers,
                              mol->mass,
                              context->linecenter,
                              context->tavg,
                              context->alpha));

        /*Calculate the molecule's optical depths and add them to existing
          values.*/
        log_info("Calculating optical depths across %d layers for molecule"
                     " %s.",
                 context->num_layers,
                 mol->name);
        check(calc_optical_depth(mol->line_params.num_lines,
                                 context->num_layers,
                                 context->linecenter,
                                 context->snn,
                                 context->gamma,
                                 context->alpha,
                                 context->ns,
                                 &(context->bins),
                                 tau));

/*
        check(calc_optical_depth_old(mol->line_params.num_lines,
                                     context->num_layers,
                                     context->linecenter,
                                     context->snn,
                                     context->gamma,
                                     context->alpha,
                                     context->ns,
                                     &(context->bins),
                                     tau));
*/

        if (context->use_h2o_ctm && mol->id == H2O)
        {
            /*Calculate the water vapor continuum optical depths.*/
            log_info("Calculating optical depth contribution due to the"
                         " water vapor continuum across %d layers.",
                     context->num_layers);
            calc_water_vapor_ctm_optical_depth_h(context->bins.num_wpoints,
                                                 context->num_layers,
                                                 tau,
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
            calc_ozone_ctm_optical_depth_h(context->bins.num_wpoints,
                                           context->num_layers,
                                           context->o3_cc.cross_section,
                                           context->ns,
                                           tau);
        }
    }

    /*Interpolate line wing optical depth contributions.*/
    log_info("Interpolating line wing optical depth contributions across"
                 " %d layers.",
             context->num_layers);
    check(interpolate(&(context->bins),
                      tau));
    return SUCCESS;
}
