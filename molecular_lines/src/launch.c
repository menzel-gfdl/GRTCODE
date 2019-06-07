#include <string.h>
#include "cfcs.h"
#include "collision_induced_absorption.h"
#include "curtis_godson.h"
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


/** @brief Driver for optical depth calculation.
    @return RS_SUCCESS or an error code.*/
int launch(MolecularLines_t * const ml, /**< Molecular lines object.*/
           fp_t *p, /**< Pressure [atm] (level).*/
           fp_t *t, /**< Temperature [K] (level).*/
           fp_t * const tau /**< Optical depth (level, wavenumber).*/
          )
{
    not_null(ml);
    not_null(p);
    not_null(t);
    not_null(tau);

    /*Set pointers to input data.*/
    if (ml->device == HOST_ONLY)
    {
        ml->p = p;
        ml->t = t;
        ml->tau = tau;
    }
    else
    {
        gmemcpy(ml->p, p, ml->num_levels, ml->device, FROM_HOST);
        gmemcpy(ml->t, t, ml->num_levels, ml->device, FROM_HOST);
    }

    /*Zero out buffers used to accumulate results.*/
    gmemset(ml->tau, 0 ,ml->num_layers*ml->grid.n, ml->device);
    gmemset(ml->bins.tau, 0, ml->bins.isize*ml->bins.num_layers,
            ml->bins.gpu_id);

    /*Calculate the total number density of air molecules integrated across
      each layer.*/
    glaunch(calc_number_densities, ml->num_layers, ml->device, ml->num_layers,
            ml->p, ml->n);

    /*Calculate integrated average layer quantities.*/
    glaunch(calc_pressures_and_temperatures, ml->num_layers, ml->device, ml->num_layers,
            ml->p, ml->t, ml->pavg, ml->tavg);

    /*Loop over the molecules and calculate the optical depths.*/
    int m;
    for (m=0; m<ml->num_molecules; ++m)
    {
        Molecule_t *mol = &(ml->mols[m]);
        int index;
        catch(molecule_hash(mol->id, &index));
        char *mesg = "Calculating spectra for %s.";
        log_info(mesg, mol->name);

        /*Calculate the integrated average layer partial pressure.*/
        fp_t const *xp = &(ml->x[index*ml->num_levels]);
        glaunch(calc_partial_pressures_and_number_densities, ml->num_layers, ml->device,
                ml->num_layers, ml->p, xp, ml->n, ml->psavg, ml->ns);

        /*Calculate pressure shifted line center positions.*/
        glaunch(calc_line_centers, mol->line_params.num_lines, ml->device,
                mol->line_params.num_lines, ml->num_layers, mol->line_params.vnn,
                mol->line_params.d, ml->pavg, ml->linecenter);

        /*Calculate total partition functions.*/
        glaunch(calc_partition_functions, mol->num_isotopologues, ml->device,
                ml->num_layers, mol->id, mol->num_isotopologues, ml->tavg, mol->q);

        /*Calculate temperature-corrected line strengths.*/
        glaunch(calc_line_strengths, mol->line_params.num_lines, ml->device,
                mol->line_params.num_lines, ml->num_layers, mol->num_isotopologues,
                mol->line_params.iso, mol->line_params.snn, mol->line_params.vnn,
                mol->line_params.en, ml->tavg, mol->q, ml->snn);

        /*Calcluate temperature and pressure corrected lorentz half-widths.*/
        glaunch(calc_lorentz_hw, mol->line_params.num_lines, ml->device,
                mol->line_params.num_lines, ml->num_layers, mol->line_params.n,
                mol->line_params.yair, mol->line_params.yself, ml->tavg, ml->pavg,
                ml->psavg, ml->gamma);

        /*Calculate doppler half-widths.*/
        glaunch(calc_doppler_hw, mol->line_params.num_lines, ml->device,
                mol->line_params.num_lines, ml->num_layers, mol->mass,
                ml->linecenter, ml->tavg, ml->alpha);

        /*Calculate the molecule's optical depths and add them to existing
          values.*/
        switch (ml->optical_depth_method)
        {
            case wavenumber_sweep:
                glaunch(sort_lines, ml->num_layers, ml->device, mol->line_params.num_lines,
                        ml->num_layers, ml->linecenter, ml->snn, ml->gamma, ml->alpha);
                glaunch(calc_optical_depth_bin_sweep, ml->bins.n, ml->device,
                        mol->line_params.num_lines, ml->num_layers, ml->linecenter,
                        ml->snn, ml->gamma, ml->alpha, ml->ns, ml->bins, ml->tau);
                break;
            case line_sweep:
                glaunch(calc_optical_depth_line_sweep, mol->line_params.num_lines,
                        ml->device, mol->line_params.num_lines, ml->num_layers,
                        ml->linecenter, ml->snn, ml->gamma, ml->alpha, ml->ns,
                        ml->bins, ml->tau);
                break;
            case line_sample:
                glaunch(calc_optical_depth_line_sample, mol->line_params.num_lines,
                        ml->device, mol->line_params.num_lines, ml->num_layers,
                        ml->linecenter, ml->snn, ml->gamma, ml->alpha, ml->ns,
                        ml->bins, ml->tau);
                break;
        }

        if (ml->use_h2o_ctm && mol->id == H2O)
        {
            /*Calculate the water vapor continuum optical depths.*/
            glaunch(calc_water_vapor_ctm_optical_depth, ml->bins.num_wpoints,
                    ml->device, ml->bins.num_wpoints, ml->num_layers,
                    ml->tau, ml->h2o_cc.coefs[MTCKD25_S296], ml->tavg,
                    ml->psavg, ml->ns, ml->h2o_cc.coefs[CKDS],
                    ml->h2o_cc.coefs[MTCKD25_F296], ml->pavg, ml->h2o_cc.coefs[CKDF]);
        }
        else if (ml->use_o3_ctm && mol->id == O3)
        {
            /*Calculate the ozone continuum optical depths.*/
            glaunch(calc_ozone_ctm_optical_depth, ml->bins.num_wpoints,
                    ml->device, ml->bins.num_wpoints, ml->num_layers,
                    ml->o3_cc.cross_section, ml->ns, ml->tau);
        }
    }

    for (m=0; m<ml->num_cfcs; ++m)
    {
        CfcCrossSection_t *cfc = &(ml->cfcs[m]);
        int index = cfc->id;
        fp_t const *xp = &(ml->x_cfc[index*ml->num_levels]);
        char *mesg = "Calculating spectra for %s.";
        log_info(mesg, cfc->name);

        /*Calculate CFC optical depths.*/
        glaunch(calc_cfc_optical_depth, ml->bins.num_wpoints, ml->device,
                ml->bins.num_wpoints, ml->num_layers, ml->n, xp,
                cfc->cross_section, ml->tau);
    }

    for (m=0; m<ml->num_cias; ++m)
    {
        CollisionInducedAbsorption_t *cia = &(ml->cia[m]);
        int index1 = cia->id[0];
        fp_t const *xp1 = &(ml->x_cia[index1*ml->num_levels]);
        int index2 = cia->id[1];
        fp_t const *xp2 = &(ml->x_cia[index2*ml->num_levels]);
        char *mesg = "Calculating CIA spectra for %s - %s.";
        log_info(mesg, cia->name[0], cia->name[1]);

        /*Calculate collision-induced absorption optical depths.*/
        glaunch(calc_cia_optical_depth, ml->bins.num_wpoints, ml->device,
                ml->bins.num_wpoints, ml->num_layers, ml->p, ml->tavg, xp1, xp2,
                cia->cross_section, ml->tau);
    }

    if (ml->optical_depth_method != line_sample)
    {
        /*Interpolate line wing optical depth contributions.*/
        glaunch(interpolate, ml->bins.n-1, ml->device, ml->bins, ml->tau);
        glaunch(interpolate_last_bin, ml->num_layers, ml->device, ml->bins, ml->tau);
    }

    if (ml->device != HOST_ONLY)
    {
        gmemcpy(tau, ml->tau, ml->num_layers*ml->grid.n, ml->device, FROM_DEVICE);
    }
    return RS_SUCCESS;
}
