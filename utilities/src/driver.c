/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "argparse.h"
#include "driver.h"
#include "gas_optics.h"
#include "grtcode_utilities.h"
#include "longwave.h"
#include "rayleigh.h"
#include "shortwave.h"
#include "solar_flux.h"


#define catch(e) { \
    int e_ = e; \
    if (e_ != GRTCODE_SUCCESS) { \
        fprintf(stderr, "[%s, %d] Error:\n", __FILE__, __LINE__); \
        char b_[1024]; \
        grtcode_errstr(e_, b_, 1024); \
        fprintf(stderr, "%s", b_); \
        return EXIT_FAILURE; \
    }}


/*Integrate using simple trapezoids on a uniform grid.*/
static void integrate(fp_t const * const in, /*Data to be integrated.*/
                      uint64_t const n, /*Size of input data array.*/
                      fp_t const dx, /*Grid spacing.*/
                      fp_t * const out /*Result of integral.*/
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


/*Add molecules, CFCs, and collision-induced absoprtion.*/
static int add_molecules(GasOptics_t * const lbl, /*Gas optics object.*/
                         Atmosphere_t const atm /*Atmospheric state.*/
                        )
{
    int i;
    for (i=0; i<atm.num_molecules; ++i)
    {
        catch(add_molecule(lbl, atm.molecules[i], NULL, NULL));
    }
    for (i=0; i<atm.num_cfcs; ++i)
    {
        catch(add_cfc(lbl, atm.cfc[i].id, atm.cfc[i].path));
    }
    for (i=0; i<atm.num_cias; ++i)
    {
        catch(add_cia(lbl, atm.cia[i].id[0], atm.cia[i].id[1], atm.cia[i].path));
    }
    return GRTCODE_SUCCESS;
}


/*Calculate optics.*/
static int calculate_optics(GasOptics_t * const lbl, /*Gas optics object.*/
                            Atmosphere_t const atm, /*Atmospheric state.*/
                            int const column, /*Column index.*/
                            int const offset, /*Additional index offset.*/
                            Optics_t * const gas, /*Optics due to molecular lines.*/
                            Optics_t * const rayleigh, /*Optics due to rayleigh scattering.*/
                            Optics_t * const aerosol, /*Optics due to aerosols.*/
                            Optics_t * const cloud, /*Optics due to clouds.*/
                            Optics_t * total /*Sum of all optics.*/
                           )
{
    int i = column + offset;
    fp_t *level_pressure = &(atm.level_pressure[i*atm.num_levels]);
    fp_t *level_temperature = &(atm.level_temperature[i*atm.num_levels]);
    int j;
    for (j=0; j<atm.num_molecules; ++j)
    {
        fp_t const *ppmv = atm.ppmv[j];
        ppmv = &(ppmv[i*atm.num_levels]);
        catch(set_molecule_ppmv(lbl, atm.molecules[j], ppmv));
    }
    for (j=0; j<atm.num_cfcs; ++j)
    {
        fp_t const *ppmv = atm.cfc_ppmv[j];
        ppmv = &(ppmv[i*atm.num_levels]);
        catch(set_cfc_ppmv(lbl, atm.cfc[j].id, ppmv));
    }
    for (j=0; j<atm.num_cia_species; ++j)
    {
        fp_t const *ppmv = atm.cia_ppmv[j];
        ppmv = &(ppmv[i*atm.num_levels]);
        catch(set_cia_ppmv(lbl, atm.cia_species[j], ppmv));
    }
    catch(calculate_optical_depth(lbl, level_pressure, level_temperature, gas));
    catch(rayleigh_scattering(rayleigh, level_pressure));
    Optics_t *optics[4] = {gas, rayleigh};
    int num_optics = 2;
    if (!atm.clean)
    {
        int n = atm.num_layers*aerosol->grid.n;
        fp_t *tau = (fp_t *)malloc(sizeof(*tau)*n);
        fp_t *omega = (fp_t *)malloc(sizeof(*omega)*n);
        fp_t *g = (fp_t *)malloc(sizeof(*g)*n);
        for (j=0; j<atm.num_layers; ++j)
        {
            n = i*atm.num_layers*atm.aerosol_grid_size + j*atm.aerosol_grid_size;
            catch(interpolate_to_grid(aerosol->grid, atm.aerosol_grid, &(atm.aerosol_optical_depth[n]),
                                      atm.aerosol_grid_size, &(tau[j*aerosol->grid.n]),
                                      linear_sample, NULL));
            catch(interpolate_to_grid(aerosol->grid, atm.aerosol_grid, &(atm.aerosol_single_scatter_albedo[n]),
                                      atm.aerosol_grid_size, &(omega[j*aerosol->grid.n]),
                                      linear_sample, NULL));
            catch(interpolate_to_grid(aerosol->grid, atm.aerosol_grid, &(atm.aerosol_asymmetry_factor[n]),
                                      atm.aerosol_grid_size, &(g[j*aerosol->grid.n]),
                                      linear_sample, NULL));
        }
        catch(update_optics(aerosol, tau, omega, g));
        optics[num_optics] = aerosol;
        num_optics++;
        free(tau);
        free(omega);
        free(g);
    }
    if (!atm.clear)
    {
        catch(cloud_optics(cloud, atm.liquid_water_path,
                           atm.liquid_water_droplet_radius, atm.cloud_param));
        optics[num_optics] = cloud;
        num_optics++;
    }
    catch(add_optics(optics, num_optics, total));
    return GRTCODE_SUCCESS;
}


static int driver(Atmosphere_t const atm, /*Atmospheric state.*/
                  char const * const hitran_path, /*Path to HITRAN database ascii file.*/
                  char const * const solar_flux_path, /*Path to solar flux ascii file.*/
                  SpectralGrid_t const lw_grid, /*Longwave spectral grid.*/
                  SpectralGrid_t const sw_grid, /*Shortwave spectral grid.*/
                  Device_t const device, /*Device object.*/
                  Output_t * const output /*Output object.*/
                 )
{
    /*Create a grid that spans the entire wavenumber range.*/
    SpectralGrid_t grid;
    double w0 = lw_grid.w0 < sw_grid.w0 ? lw_grid.w0 : sw_grid.w0;
    double wn = sw_grid.wn > lw_grid.wn ? sw_grid.wn : lw_grid.wn;
    double res = sw_grid.dw > lw_grid.dw ? sw_grid.dw : lw_grid.dw;
    catch(create_spectral_grid(&grid, w0, wn, res));

    /*Intialize gas optics objects.*/
    GasOptics_t lbl_lw;
    int method = line_sweep;
    catch(create_gas_optics(&lbl_lw, atm.num_levels, &lw_grid, &device,
                            hitran_path, atm.h2o_ctm, atm.o3_ctm, NULL, &method));
    catch(add_molecules(&lbl_lw, atm));
    GasOptics_t lbl_sw;
    catch(create_gas_optics(&lbl_sw, atm.num_levels, &sw_grid, &device,
                            hitran_path, atm.h2o_ctm, atm.o3_ctm, NULL, &method));
    catch(add_molecules(&lbl_sw, atm));

    /*Initialize optics objects.*/
    Optics_t optics_lw_gas;
    catch(create_optics(&optics_lw_gas, atm.num_layers, &lw_grid, &device));
    Optics_t optics_sw_gas;
    catch(create_optics(&optics_sw_gas, atm.num_layers, &sw_grid, &device));
    Optics_t optics_lw_rayleigh;
    catch(create_optics(&optics_lw_rayleigh, atm.num_layers, &lw_grid, &device));
    Optics_t optics_sw_rayleigh;
    catch(create_optics(&optics_sw_rayleigh, atm.num_layers, &sw_grid, &device));
    Optics_t optics_lw_aerosol;
    Optics_t optics_sw_aerosol;
    if (!atm.clean)
    {
        catch(create_optics(&optics_lw_aerosol, atm.num_layers, &lw_grid, &device));
        catch(create_optics(&optics_sw_aerosol, atm.num_layers, &sw_grid, &device));
    }
    Optics_t optics_lw_cloud;
    Optics_t optics_sw_cloud;
    if (!atm.clear)
    {
        catch(create_optics(&optics_lw_cloud, atm.num_layers, &lw_grid, &device));
        catch(create_optics(&optics_sw_cloud, atm.num_layers, &sw_grid, &device));
    }
    Optics_t optics;
    catch(create_optics(&optics, atm.num_layers, &grid, &device));

    /*Initialize the solar fluxe object.*/
    SolarFlux_t solar_flux;
    catch(create_solar_flux(&solar_flux, &grid, solar_flux_path));

    /*Initialize solver objects.*/
    Longwave_t longwave;
    catch(create_longwave(&longwave, atm.num_levels, &lw_grid, &device));
    Shortwave_t shortwave;
    catch(create_shortwave(&shortwave, atm.num_levels, &grid, &device));

    /*Initialize other buffers.*/
    fp_t *albedo_dir = (fp_t *)malloc(sizeof(*albedo_dir)*grid.n);
    fp_t *surface_emissivity = (fp_t *)malloc(sizeof(*surface_emissivity)*lw_grid.n);
    fp_t *lw_flux_up = (fp_t *)malloc(sizeof(*lw_flux_up)*atm.num_levels*lw_grid.n);
    fp_t *lw_flux_down = (fp_t *)malloc(sizeof(*lw_flux_down)*atm.num_levels*lw_grid.n);
    fp_t *sw_flux_up = (fp_t *)malloc(sizeof(*sw_flux_up)*atm.num_levels*grid.n);
    fp_t *sw_flux_down = (fp_t *)malloc(sizeof(*sw_flux_down)*atm.num_levels*grid.n);

    /*Loop through the columns.*/
    int t;
    for (t=0; t<atm.num_times; ++t)
    {
        int offset = t*atm.num_columns;
        int i;
        for (i=0; i<atm.num_columns; ++i)
        {
            /*Longwave.*/
            Optics_t optics_lw_total;
            catch(calculate_optics(&lbl_lw, atm, i, offset, &optics_lw_gas,
                                   &optics_lw_rayleigh, &optics_lw_aerosol,
                                   &optics_lw_cloud, &optics_lw_total));
            fp_t surface_temperature = atm.surface_temperature[offset+i];
            fp_t *layer_temperature = &(atm.layer_temperature[(offset+i)*atm.num_layers]);
            fp_t *level_temperature = &(atm.level_temperature[(offset+i)*atm.num_levels]);
            int n = (offset + i)*atm.emissivity_grid_size;
            catch(interpolate_to_grid(lw_grid, atm.emissivity_grid, &(atm.surface_emissivity[n]),
                                      atm.emissivity_grid_size, surface_emissivity,
                                      linear_sample, constant_extrapolation));
            catch(calculate_lw_fluxes(&longwave, &optics_lw_total, surface_temperature,
                                      layer_temperature, level_temperature,
                                      surface_emissivity, lw_flux_up, lw_flux_down));
            fp_t flux_up_total[atm.num_levels];
            fp_t flux_down_total[atm.num_levels];
            int j;
            for (j=0; j<atm.num_levels; ++j)
            {
                integrate(&(lw_flux_up[j*lw_grid.n]), lw_grid.n, lw_grid.dw, &(flux_up_total[j]));
                integrate(&(lw_flux_down[j*lw_grid.n]), lw_grid.n, lw_grid.dw, &(flux_down_total[j]));
            }
            write_output(output, RLU, flux_up_total, t, i);
            write_output(output, RLD, flux_down_total, t, i);

            /*Shortwave.*/
            fp_t const zen_dir = atm.solar_zenith_angle[offset+i];
            if (zen_dir > 0.)
            {
                Optics_t optics_sw_total;
                catch(calculate_optics(&lbl_sw, atm, i, offset, &optics_sw_gas,
                                       &optics_sw_rayleigh, &optics_sw_aerosol,
                                       &optics_sw_cloud, &optics_sw_total));
                /*Combine the longwave and shortwave optics.*/
                catch(sample_optics(&optics, &optics_lw_total, &lw_grid.w0, &lw_grid.wn));
                catch(sample_optics(&optics, &optics_sw_total, &sw_grid.w0, &sw_grid.wn));
                fp_t const zen_dif = 0.5;
                int n = (offset + i)*atm.albedo_grid_size;
                catch(interpolate_to_grid(grid, atm.albedo_grid, &(atm.surface_albedo[n]),
                                          atm.albedo_grid_size, albedo_dir, linear_sample,
                                          constant_extrapolation));
                fp_t *albedo_dif = albedo_dir;
                catch(calculate_sw_fluxes(&shortwave, &optics, zen_dir, zen_dif,
                                          albedo_dir, albedo_dif, atm.total_solar_irradiance[offset+i],
                                          solar_flux.incident_flux, sw_flux_up, sw_flux_down));
                catch(destroy_optics(&optics_sw_total));

                /*Integrate fluxes and write them to the output file.*/
                for (j=0; j<atm.num_levels; ++j)
                {
                    integrate(&(sw_flux_up[j*grid.n]), grid.n, grid.dw, &(flux_up_total[j]));
                    integrate(&(sw_flux_down[j*grid.n]), grid.n, grid.dw, &(flux_down_total[j]));
                }
                write_output(output, RSU, flux_up_total, t, i);
                write_output(output, RSD, flux_down_total, t, i);
            }
            catch(destroy_optics(&optics_lw_total));
        }
    }

    /*Release memory.*/
    free(albedo_dir);
    free(surface_emissivity);
    free(lw_flux_up);
    free(lw_flux_down);
    free(sw_flux_up);
    free(sw_flux_down);
    catch(destroy_shortwave(&shortwave));
    catch(destroy_longwave(&longwave));
    catch(destroy_solar_flux(&solar_flux));
    catch(destroy_optics(&optics_lw_gas));
    catch(destroy_optics(&optics_sw_gas));
    catch(destroy_optics(&optics_lw_rayleigh));
    catch(destroy_optics(&optics_sw_rayleigh));
    if (!atm.clean)
    {
        catch(destroy_optics(&optics_lw_aerosol));
        catch(destroy_optics(&optics_sw_aerosol));
    }
    Optics_t optics_clouds;
    if (!atm.clear)
    {
        catch(destroy_optics(&optics_lw_cloud));
        catch(destroy_optics(&optics_sw_cloud));
    }
    catch(destroy_optics(&optics));
    catch(destroy_gas_optics(&lbl_lw));
    catch(destroy_gas_optics(&lbl_sw));
    return GRTCODE_SUCCESS;
}


/*Main driver program.  When linking an executable, the use must provide a c file that
  includes driver.h in this directory and provides implementations for:
   - struct output;
   - void close_flux_file(Output_t *output);
   - Atmosphere_t create_atmosphere(Parser_t parser);
   - void destroy_atmosphere(Atmosphere_t *atm);
   - Output_t create_flux_file(char *path, Atmosphere_t *atm);
   - void write_output(Output_t *output, enum varid id, fp_t *data, int time, int column);*/
int main(int argc, char **argv)
{
    /*Add/parse command line arguments.*/
    char *description = "Calculates radiative fluxes using the line-by-line method.";
    Parser_t parser = create_parser(argc, argv, description);
    add_argument(&parser, "hitran_file", NULL, "HITRAN database file.", NULL);
    add_argument(&parser, "solar_flux", NULL, "Solar flux file.", NULL);
    int one = 1;
    add_argument(&parser, "-c", "--line-cutoff", "Cutoff [1/cm] from line center.", &one);
    add_argument(&parser, "-d", "--device", "GPU id", &one);
    add_argument(&parser, "-o", NULL, "Name of output file.", &one);
    add_argument(&parser, "-r-lw", "--lw-resolution", "Longwave spectral resolution [1/cm].", &one);
    add_argument(&parser, "-r-sw", "--sw-resolution", "Shortwave spectral resolution [1/cm].", &one);
    add_argument(&parser, "-s", "--solver", "Shortwave solver.", &one);
    add_argument(&parser, "-v", "--verbose", "Increase verbosity.", NULL);
    add_argument(&parser, "-w-lw", "--lw-lower-bound", "Longwave spectral lower bound [1/cm].", &one);
    add_argument(&parser, "-w-sw", "--sw-lower-bound", "Shortwave spectral lower bound [1/cm].", &one);
    add_argument(&parser, "-W-lw", "--lw-upper-bound", "Longwave spectral upper bound [1/cm].", &one);
    add_argument(&parser, "-W-sw", "--sw-upper-bound", "Shortwave spectral upper bound [1/cm].", &one);

    /*Get the atmospheric data.*/
    Atmosphere_t atm = create_atmosphere(&parser);

    /*Set GRTCODE verbosity.*/
    int verbosity_level = get_argument(parser, "-v", NULL) ? GRTCODE_INFO : GRTCODE_WARN;
    grtcode_set_verbosity(verbosity_level);

    /*Get paths to required input files.*/
    char hitran_path[valuelen];
    get_argument(parser, "hitran_file", hitran_path);
    char solar_flux_path[valuelen];
    get_argument(parser, "solar_flux", solar_flux_path);

    /*Create a spectral grids.*/
    char buffer[valuelen];
    double w0 = get_argument(parser, "-w-lw", buffer) ? atof(buffer) : 1.;
    double wn = get_argument(parser, "-W-lw", buffer) ? atof(buffer) : 3250.;
    double dw = get_argument(parser, "-r-lw", buffer) ? atof(buffer) : 0.1;
    SpectralGrid_t lw_grid;
    catch(create_spectral_grid(&lw_grid, w0, wn, dw));
    dw = get_argument(parser, "-r-sw", buffer) ? atof(buffer) : 1.;
    w0 = get_argument(parser, "-w-sw", buffer) ? atof(buffer) : dw + wn;
    wn = get_argument(parser, "-W-sw", buffer) ? atof(buffer) : 50000.;
    SpectralGrid_t sw_grid;
    catch(create_spectral_grid(&sw_grid, w0, wn, dw));

    /*Set the device to run on.*/
    Device_t device;
    int *device_id = NULL;
    if (get_argument(parser, "-d", buffer))
    {
        int d = atoi(buffer);
        device_id = &d;
    }
    catch(create_device(&device, device_id));

    /*Initialize the output file.*/
    if (!get_argument(parser, "-o", buffer))
    {
        snprintf(buffer, valuelen, "%s", "output.nc");
    }
    Output_t *output;
    create_flux_file(&output, buffer, &atm);

    /*Calculate the fluxes over all the columns.*/
    catch(driver(atm, hitran_path, solar_flux_path, lw_grid, sw_grid, device, output));

    /*Clean up.*/
    close_flux_file(output);
    destroy_atmosphere(&atm);
    destroy_parser(&parser);
    return EXIT_SUCCESS;
}
