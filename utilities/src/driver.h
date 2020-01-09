#ifndef DRIVER_H_
#define DRIVER_H_

#include "argparse.h"
#include "cloud_optics.h"
#include "grtcode_utilities.h"


typedef enum VarId
{
    RLD = 0,
    RLU,
    RSU,
    RSD,
    NUM_VARS
} VarId_t;


typedef struct Cfc
{
    int id; /*GRTCODE CFC id.*/
    char path[valuelen]; /*Path to CFC cross-section file.*/
} Cfc_t;


typedef struct Cia
{
    int id[2]; /*GRTCODE CIA id of species.*/
    char path[valuelen]; /*Path to CIA cross-section file.*/
} Cia_t;


typedef struct Atmosphere
{
    fp_t *aerosol_grid; /*Wavenumber [cm-1] grid for the aerosol optical properties (wavenumber).*/
    int aerosol_grid_size; /*Length of aerosol_grid array.*/
    fp_t *aerosol_optical_depth; /*Aerosol optical depth (time, column, layer, wavenumber).*/
    fp_t *aerosol_single_scatter_albedo; /*Aerosol single-scatter albedo (time, column, layer, wavenumber).*/
    fp_t *aerosol_asymmetry_factor; /*Aerosol asymmetry factor (time, column, layer, wavenumber).*/
    fp_t *albedo_grid; /*Wavenumber [cm-1] grid for the surface albedo (wavenumber).*/
    size_t albedo_grid_size; /*Length of albedo_grid array.*/
    Cfc_t *cfc; /*GRTCODE cfc ids and paths.*/
    fp_t **cfc_ppmv; /*CFC abundance [ppmv] (CFC, time, column, level).*/
    LiquidCloud_t cloud_param; /*Cloud parameterization.*/
    Cia_t *cia; /*GRTCODE CIA ids and paths.*/
    fp_t **cia_ppmv; /*CIA abundance [ppmv] (molecule, time, column, level).*/
    int *cia_species; /*GRTCODE CIA ids.*/
    int clean; /*Flag for clean (aerosol-free) sky.*/
    int clear; /*Flag for clear (cloud-free) sky.*/
    fp_t *emissivity_grid; /*Wavenumber [cm-1] grid for the surface emissivity (wavenumber).*/
    size_t emissivity_grid_size; /*Length of emissivity_grid array.*/
    char h2o_ctm[valuelen]; /*Path to water vapor continuum directory.*/
    fp_t *layer_pressure; /*Pressure [atm] (time, column, layer).*/
    fp_t *layer_temperature; /*Temperature [K] (time, column, layer).*/
    fp_t *level_pressure; /*Pressure [atm] (time, column, level).*/
    fp_t *level_temperature; /*Temperature [K] (time, column, level).*/
    fp_t *liquid_water_droplet_radius; /*Liquid water equivalent radius [microns] (time, column, layer).*/
    fp_t *liquid_water_path; /*Liquid water path [g/m2] (time, column, layer).*/
    int *molecules; /*GRTCODE molecule ids.*/
    int num_cfcs; /*Length of cfc array.*/
    int num_cia_species; /*Length of cia_species array.*/
    int num_cias; /*Length of cia array.*/
    int num_columns; /*Number of columns.*/
    int num_layers; /*Number of atmospheric layers.*/
    int num_levels; /*Number of atmospheric levels.*/
    int num_molecules; /*Length of molecules array.*/
    int num_times; /*Number of times.*/
    char o3_ctm[valuelen]; /*Path to ozone continuum file.*/
    fp_t **ppmv; /*Molecular abundance [ppmv] (molecule, time, column, level).*/
    fp_t *solar_zenith_angle; /*Cosine of solar zenith angle (time, column).*/
    fp_t *surface_albedo; /*Surface albedo (time, column, wavenumber).*/
    fp_t *surface_emissivity; /*Surface emissivity (time, column, wavenumber).*/
    fp_t *surface_temperature; /*Surface temperature [K] (time, column).*/
    fp_t *total_solar_irradiance; /*Total solar irradiance [W m-2] at TOA (time, column).*/
} Atmosphere_t;


typedef struct Output Output_t;


void close_flux_file(Output_t * const output);


void create_flux_file(Output_t **output, char const * const path, Atmosphere_t const * const atm);


Atmosphere_t create_atmosphere(Parser_t * const parser);


void destroy_atmosphere(Atmosphere_t *atm);


void write_output(Output_t *output, VarId_t id, fp_t *data, int time, int column);


#endif
