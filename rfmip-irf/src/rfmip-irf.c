#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "atmosphere.h"
#include "argparse.h"
#include "device.h"
#include "floating_point_type.h"
#include "longwave.h"
#include "molecular_lines.h"
#include "optics.h"
#include "rayleigh.h"
#include "return_codes.h"
#include "shortwave.h"
#include "solar_flux.h"
#include "spectral_grid.h"
#include "verbosity.h"


#define catch(e) { \
    if (e != RS_SUCCESS) { \
        fprintf(stderr, "[%s, %d] Error.\n", __FILE__, __LINE__); \
        return EXIT_FAILURE; \
    }}
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/** @brief Integrate using simple trapezoids on a uniform grid.*/
static void integrate(fp_t const * const in, /**< Data to be integrated.*/
                      uint64_t const n, /**< Size of input data array.*/
                      fp_t const dx, /**< Grid spacing.*/
                      fp_t * const out /**< Result of integral.*/
                     )
{
    *out = 0.;
    uint64_t i;
    for (i=0; i<n-1; ++i)
    {
        *out += dx*0.5*(in[i] + in[i+1]);
    }
    return;
}


/** @brief Set a species as active.*/
static void activate_species(Parser_t const parser, /**< Parser object.*/
                             char const *const arg, /**< Arg to check for.*/
                             int * const array, /**< Array.*/
                             int const tag, /**< Value to store in array.*/
                             int * spot, /**< Current spot in array.*/
                             int const max_size, /**< Size of input array.*/
                             char **values /**< Array to store argument values.*/
                            )
{
    char buffer[valuelen];
    if (get_argument(parser, arg, buffer))
    {
        array[*spot] = tag;
        if (values != NULL)
        {
            snprintf(values[*spot], valuelen, "%s", buffer);
        }
        *spot += 1;
        if (*spot >= max_size)
        {
            fprintf(stderr, "Array is too small, increase size.\n");
            exit(EXIT_FAILURE);
        }
    }
    return;
}


/*Calculate the radiative fluxes for the CMIP6 RFMIP-IRF test cases.*/
int main(int argc, char **argv)
{
    /*Add/parse command line arguments.*/
    char *description = "Calculates the radiative fluxes for the CMIP6 RFMIP-IRF test"
                        " cases.";
    Parser_t parser = create_parser(argc, argv, description);
    add_argument(&parser, "input_file", NULL, "Input data file.", NULL);
    add_argument(&parser, "experiment", NULL, "Experiment number.", NULL);
    add_argument(&parser, "hitran_file", NULL, "HITRAN database file.", NULL);
    add_argument(&parser, "solar_flux", NULL, "Solar flux CSV file.", NULL);
    int one = 1;
    add_argument(&parser, "-CH4", NULL, "Include CH4.", NULL);
    add_argument(&parser, "-CO", NULL, "Include CO.", NULL);
    add_argument(&parser, "-CO2", NULL, "Include CO2.", NULL);
    add_argument(&parser, "-F11", NULL, "CSV file with F11 cross sections.", &one);
    add_argument(&parser, "-F12", NULL, "CSV file with F12 cross sections.", &one);
    add_argument(&parser, "-H2O", NULL, "Include H2O.", NULL);
    add_argument(&parser, "-N2O", NULL, "Include N2O.", NULL);
    add_argument(&parser, "-O2", NULL, "Include O2.", NULL);
    add_argument(&parser, "-O3", NULL, "Include O3.", NULL);
    add_argument(&parser, "-c", "--line-cutoff", "Cutoff [1/cm] from line center.", &one);
    add_argument(&parser, "-h2o-ctm", NULL, "Directory containing H2O continuum files", &one);
    add_argument(&parser, "-o", NULL, "Name of output file.", &one);
    add_argument(&parser, "-o3-ctm", NULL, "Directory containing O3 continuum files", &one);
    add_argument(&parser, "-r", "--spectral-resolution", "Spectral resolution [1/cm].", &one);
    add_argument(&parser, "-v", "--verbose", "Increase verbosity.", NULL);
    add_argument(&parser, "-w", "--spectral-lower-bound", "Spectral lower bound [1/cm].", &one);
    add_argument(&parser, "-W", "--spectral-upper-bound", "Spectral upper bound [1/cm].", &one);
    add_argument(&parser, "-x", "--column-lower-bound", "Starting column index.", &one);
    add_argument(&parser, "-X", "--column-upper-bound", "Ending column index.", &one);
    add_argument(&parser, "-z", "--level-lower-bound", "Starting level index.", &one);
    add_argument(&parser, "-Z", "--level-upper-bound", "Ending level index.", &one);
    parse_args(parser);

    /*Set verbosity.*/
    char buffer[valuelen];
    if (get_argument(parser, "-v", NULL))
    {
        rs_set_verbosity(RS_INFO);
    }
    else
    {
        rs_set_verbosity(RS_WARN);
    }

    /*Set device.*/
    Device_t device;
    catch(create_device(&device, NULL));

    /*Create a spectral grid.*/
    double const w0 = 1.;
    double const wn = 50000.;
    double const dw = 0.1;
    SpectralGrid_t grid;
    catch(create_spectral_grid(&grid, w0, wn, dw));

    /*Determine which molecules to use.*/
    int molecules[32];
    int num_molecules = 0;
    activate_species(parser, "-CH4", molecules, CH4, &num_molecules, 32, NULL);
    activate_species(parser, "-CO", molecules, CO, &num_molecules, 32, NULL);
    activate_species(parser, "-CO2", molecules, CO2, &num_molecules, 32, NULL);
    activate_species(parser, "-H2O", molecules, H2O, &num_molecules, 32, NULL);
    activate_species(parser, "-N2O", molecules, N2O, &num_molecules, 32, NULL);
    activate_species(parser, "-O2", molecules, O2, &num_molecules, 32, NULL);
    activate_species(parser, "-O3", molecules, O3, &num_molecules, 32, NULL);

    /*Determine which CFCs to use.*/
    int cfcs[32];
    char *cfc_paths[32];
    int i;
    for (i=0; i<32; ++i)
    {
        cfc_paths[i] = malloc(sizeof(*(cfc_paths[i]))*valuelen);
    }
    int num_cfcs = 0;
    activate_species(parser, "-F11", cfcs, F11, &num_cfcs, 32, cfc_paths);
    activate_species(parser, "-F12", cfcs, F12, &num_cfcs, 32, cfc_paths);

    /*Read in the atmospheric input data.*/
    Atmosphere_t atm;
    atm.num_wavenumber = grid.n;
    atm.x = 0;
    atm.num_columns = 1;
    atm.z = 0;
    atm.num_levels = 61;
    atm.num_layers = atm.num_levels - 1;
    get_argument(parser, "experiment", buffer);
    int experiment = atoi(buffer);
    get_argument(parser, "input_file", buffer);
    create_atmosphere(&atm, buffer, experiment, molecules, num_molecules, cfcs,
                      num_cfcs);

    /*Read in the incident solar flux.*/
    SolarFlux_t solar_flux;
    fp_t total_solar_irradiance = 1407.679; /*[W/m^2]*/
    get_argument(parser, "solar_flux", buffer);
    catch(create_solar_flux(&solar_flux, &grid, buffer, total_solar_irradiance));

    /*Initialize a molecular lines object.*/
    char hitran_path[valuelen];
    get_argument(parser, "hitran_file", hitran_path);
    char h2o_ctm[valuelen];
    if (!get_argument(parser, "-h2o-ctm", h2o_ctm))
    {
        snprintf(h2o_ctm, valuelen, "%s", "none");
    }
    char o3_ctm[valuelen];
    if (!get_argument(parser, "-o3-ctm", o3_ctm))
    {
        snprintf(o3_ctm, valuelen, "%s", "none");
    }
    MolecularLines_t molecular_lines;
    catch(create_molecular_lines(&molecular_lines, atm.num_levels, &grid, &device,
                                 hitran_path, h2o_ctm, o3_ctm, NULL, NULL));

    /*Add molecules and CFCs.*/
    for (i=0; i<num_molecules; ++i)
    {
        catch(grt_add_molecule(&molecular_lines, molecules[i], NULL, NULL));
    }
    for (i=0; i<num_cfcs; ++i)
    {
        catch(grt_add_cfc(&molecular_lines, cfcs[i], cfc_paths[i]));
    }

    /*Initialize an optics object.*/
    Optics_t optics_ml;
    catch(create_optics(&optics_ml, atm.num_layers, &grid, &device));
    Optics_t optics_rayleigh;
    catch(create_optics(&optics_rayleigh, atm.num_layers, &grid, &device));

    /*Initialize a longwave object.*/
    Longwave_t longwave;
    catch(create_longwave(&longwave, atm.num_levels, &grid, &device));

    /*Initialize a shortwave object.*/
    Shortwave_t shortwave;
    catch(create_shortwave(&shortwave, atm.num_levels, &grid, &device));

    /*Initialize the output file.*/
    if (!get_argument(parser, "-o", buffer))
    {
        snprintf(buffer, valuelen, "%s", "rfmip-irf.output.nc");
    }
    Output_t output = create_flux_file(buffer, &atm);

    /*Loop through the columns.*/
    fp_t *flux_up = malloc(sizeof(*flux_up)*atm.num_levels*grid.n);
    fp_t *flux_down = malloc(sizeof(*flux_down)*atm.num_levels*grid.n);
    for (i=0; i<atm.num_columns; ++i)
    {
        /*Calculate molecular spectra.*/
        fp_t *level_pressure = &(atm.level_pressure[i*atm.num_levels]);
        fp_t *level_temperature = &(atm.level_temperature[i*atm.num_levels]);
        int j;
        for (j=0; j<num_molecules; ++j)
        {
            fp_t *ppmv = atm.ppmv[j];
            ppmv = &(ppmv[i*atm.num_levels]);
            catch(grt_set_molecule_ppmv(&molecular_lines, molecules[j], ppmv));
        }
        for (j=0; j<num_cfcs; ++j)
        {
            fp_t *ppmv = atm.cfc_ppmv[j];
            ppmv = &(ppmv[i*atm.num_levels]);
            catch(grt_set_cfc_ppmv(&molecular_lines, cfcs[j], ppmv));
        }
        catch(grt_calculate_optical_depth(&molecular_lines, level_pressure,
                                          level_temperature, &optics_ml));

        /*Calculate longwave fluxes.*/
        fp_t surface_temperature = atm.surface_temperature[i];
        fp_t *layer_temperature = &(atm.layer_temperature[i*atm.num_layers]);
        fp_t *surface_emissivity = &(atm.surface_emissivity[i*atm.num_wavenumber]);
        double lw_solver_w0 = 1.;
        double lw_solver_wn = 3250.;
        SpectralGrid_t lw_solver_grid;
        catch(create_spectral_grid(&lw_solver_grid, lw_solver_w0, lw_solver_wn, grid.dw));
        catch(calculate_lw_fluxes(&longwave, &optics_ml, surface_temperature,
                                  layer_temperature, level_temperature,
                                  surface_emissivity, flux_up, flux_down,
                                  &(lw_solver_grid.w0), &(lw_solver_grid.wn)));

        /*Integrate fluxes and write them to the output file.*/
        fp_t flux_up_total[atm.num_levels];
        fp_t flux_down_total[atm.num_levels];
        for (j=0; j<atm.num_levels; ++j)
        {
            integrate(&(flux_up[j*lw_solver_grid.n]), lw_solver_grid.n, lw_solver_grid.dw,
                      &(flux_up_total[j]));
            integrate(&(flux_down[j*lw_solver_grid.n]), lw_solver_grid.n, lw_solver_grid.dw,
                      &(flux_down_total[j]));
        }
        write_fluxes(&output, RLU, i+atm.x, flux_up_total);
        write_fluxes(&output, RLD, i+atm.x, flux_down_total);

        /*Calculate the optical properities of a column.*/
        catch(rayleigh_scattering(&optics_rayleigh, level_pressure));

        /*Calculate the combined optical properties.*/
        Optics_t const * const optics_mech[2] = {&optics_ml, &optics_rayleigh};
        Optics_t optics_combined;
        catch(add_optics(optics_mech, 2, &optics_combined));

        /*Calculate shortwave fluxes.*/
        fp_t const zen_dir = atm.solar_zenith_angle[i];
        fp_t const zen_dif = 0.5;
        fp_t const albedo_dir = atm.surface_albedo[i];
        fp_t const albedo_dif = albedo_dir;
        catch(calculate_sw_fluxes(&shortwave, &optics_combined, zen_dir, zen_dif, albedo_dir,
                                  albedo_dif, solar_flux.incident_flux, flux_up, flux_down));
        catch(destroy_optics(&optics_combined));

        /*Integrate fluxes and write them to the output file.*/
        for (j=0; j<atm.num_levels; ++j)
        {
            integrate(&(flux_up[j*grid.n]), grid.n, grid.dw, &(flux_up_total[j]));
            integrate(&(flux_down[j*grid.n]), grid.n, grid.dw, &(flux_down_total[j]));
        }
        write_fluxes(&output, RSU, i+atm.x, flux_up_total);
        write_fluxes(&output, RSD, i+atm.x, flux_down_total);
    }

    /*Clean up.*/
    close_flux_file(&output);
    free(flux_up);
    free(flux_down);
    catch(destroy_shortwave(&shortwave));
    catch(destroy_longwave(&longwave));
    catch(destroy_solar_flux(&solar_flux));
    catch(destroy_optics(&optics_ml));
    catch(destroy_optics(&optics_rayleigh));
    catch(destroy_molecular_lines(&molecular_lines));
    destroy_parser(&parser);
    for (i=0; i<32; ++i)
    {
        free(cfc_paths[i]);
    }
    return EXIT_SUCCESS;
}
