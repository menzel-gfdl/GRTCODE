#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "eval_gamma.h"
#include "eval_profile.h"
#include "eval_pshift.h"
#include "eval_snn_correction.h"
#include "floating_point_type.h"
#include "host_launch.h"
#include "integrate_layer.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "parse_HITRAN_file.h"
#include "pre_eval_snn.h"
#include "utils.h"
#include "water_vapor_continuum.h"


int launch_h(int const num_levels,
             fp_t const * const P,
             fp_t const * const T,
             fp_t const * const x,
             fp_t * const Pavg,
             fp_t * const Tavg,
             fp_t * const N,
             fp_t * const Ns,
             fp_t * const Psavg,
             fp_t * const snn_ref,
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
    not_null(P);
    not_null(T);
    not_null(x);
    not_null(Pavg);
    not_null(Tavg);
    not_null(N);
    not_null(Ns);
    not_null(Psavg);
    not_null(snn_ref);
    not_null(gamma);
    not_null(Pshift);
    not_null(s);
    not_null(line_params);
    not_null(tau);

    /*Zero out buffers used to accumulate results.*/
    int num_layers = num_levels - 1;
    memset(tau,
           0,
           sizeof(*tau)*num_layers*num_wpoints);

    /*Calculate the total number density of air moleucles integrated across
      each layer.*/
    log_info("Integrating total number density across %d layers.",
             num_layers);
    check(integrated_N_h(num_layers,
                         P,
                         N));

    /*Calculate integrated average layer quantities.*/
    log_info("Calculating Curtis-Godson pressure and temperature across"
                 " %d layers.",
             num_layers);
    check(Curtis_Godson_PT_h(num_layers,
                             P,
                             T,
                             Pavg,
                             Tavg));

    /*Loop over the molecules and calculate the optical depths.*/
    int mol;
    for (mol=0;mol<num_molecules;++mol)
    {
        /*Point to the molecules line parameter structure.*/
        lines = line_params[mol];

        /*Copy the snn_ref array.*/
        unsigned int num_lines = lines->num_lines;
        memcpy(snn_ref,
               lines->snn_ref,
               sizeof(*snn_ref)*num_lines);

        /*Calculate the initial Snn_ref correction.*/
        log_info("Launching kernel pre_eval_snn_h across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        int mol_id = lines->mol;
        pre_eval_snn_h(num_lines,
                       mol_id,
                       lines->iso,
                       lines->vnn,
                       lines->en,
                       snn_ref);

        /*Calculate the integrated average layer partial pressure.*/
        fp_t const *xp = &(x[mol*num_levels]);
        log_info("Calculating Curtis-Godson partial pressure and abundance"
                     " across %d layers for molecule %d.",
                 num_layers,
                 mol);
        Curtis_Godson_PsNs_h(num_layers,
                             P,
                             xp,
                             N,
                             Psavg,
                             Ns);

        /*Calcluate the lorentz half-width at half-max (HWHM).*/
        log_info("Launching kernel eval_gamma_h across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        eval_gamma_h(num_layers,
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
        log_info("Launching kernel eval_pshift_h across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        eval_pshift_h(num_layers,
                      num_lines,
                      Pavg,
                      lines->vnn,
                      lines->d,
                      Pshift);

        /*Calculate the remainder of the Snn_ref correction.*/
        log_info("Launching kernel eval_snn_correction_h across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        eval_snn_correction_h(num_layers,
                              num_lines,
                              mol_id,
                              Tavg,
                              lines->iso,
                              lines->vnn,
                              lines->en,
                              snn_ref,
                              s);

        /*Calculate the molecule's optical depths and add them to existing
          values.*/
        log_info("Launching kernel eval_profile_h across %d layers"
                     " for molecule %d.",
                 num_layers,
                 mol);
        eval_profile_h(mol_id,
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
            not_null(h2o_cc);
            log_info("Calculating optical depth due to the water vapor"
                         " continuum across %d layers.",
                     num_layers);
            calc_water_vapor_ctm_optical_depth_h(num_wpoints,
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
            not_null(o3_cc);
            log_info("Calculating optical depth due to the ozone"
                         " continuum across %d layers.",
                     num_layers);
            calc_ozone_ctm_optical_depth_h(num_wpoints,
                                           num_layers,
                                           o3_cc->cross_section,
                                           Ns,
                                           tau);
        }
    }
    return SUCCESS;
}
