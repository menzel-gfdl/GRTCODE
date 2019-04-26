#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "floating_point_type.h"


/** @brief Atmospheric properties.*/
typedef struct Atmosphere
{
    int x; /**< Column lower bound index.*/
    int z; /**< Level lower bound index.*/
    int num_columns; /**< Number of atmospheric columns.*/
    int num_levels; /**< Number of atmopsheric levels.*/
    int num_layers; /**< Number of atmophseric layers.*/
    uint64_t num_wavenumber; /**< Number of wavenumber points.*/
    int num_molecules; /**< Number of molecules.*/
    int num_cfcs; /**< Number of CFCs.*/
    int num_cias; /**< Number of CIAs.*/
    fp_t *level_pressure; /**< Pressure [atm] (column, level).*/
    fp_t *layer_pressure; /**< Pressure [atm] (column, layer).*/
    fp_t *level_temperature; /**< Temperature [K] (column, level).*/
    fp_t *layer_temperature; /**< Temperature [K] (column, layer).*/
    fp_t *surface_temperature; /**< Surface temperature [K] (column).*/
    fp_t *total_solar_irradiance; /**< Total solar irradiance at TOA [W/m^2] (column).*/
    fp_t *solar_zenith_angle; /**< Cosine of solar zenith angle (column).*/
    fp_t *surface_albedo; /**< Surface albedo (column).*/
    fp_t *surface_emissivity; /**< Surface emissivity (column, wavenumber).*/
    fp_t **ppmv; /**< Molecular abundance [ppmv] (molecule, column, level).*/
    fp_t **cfc_ppmv; /**< CFC abundance [ppmv] (CFC, column, level).*/
    fp_t **cia_ppmv; /**< CIA abindance [ppmv] (molecule, column, level).*/
} Atmosphere_t;


/**@ brief Dimension indices.*/
enum dimensions
{
    EXPT = 0,
    SITE,
    LEVEL,
    NUM_DIMS
};


/** @brief Flux variable indices.*/
enum variables
{
    RLU = 0,
    RLD,
    RSU,
    RSD,
    NUM_VARS
};


/** @brief Output file object.*/
typedef struct Output
{
    int ncid;
    int dimid[NUM_DIMS];
    int varid[NUM_VARS];
} Output_t;


/**@ brief Reserve memory and read in atmospheric data.*/
void create_atmosphere(Atmosphere_t * const atm, /**< Atmosphere object.*/
                       char const * const filepath, /**< Input data file.*/
                       int const experiment, /**< Experiment index.*/
                       int const * const molecules, /**< Array of molecule ids.*/
                       int const num_molecules, /**< Number of molecules.*/
                       int const * const cfcs, /**< Array of CFC ids.*/
                       int const num_cfcs, /**< Number of CFCs.*/
                       int const * const cias, /**< Array of CIA ids.*/
                       int const num_cias /**< NUmber of CIAs.*/
                      );


/** @brief Free memory for atmosphere.*/
void destroy_atmosphere(Atmosphere_t * const atm /*Atmosphere object.*/
                       );


/** @brief Create an output file and write metadata.*/
Output_t create_flux_file(char const * const filepath, /**< File path.*/
                          Atmosphere_t const * const atm /**< Atmosphere object.*/
                         );


/** @brief Close output file.*/
void close_flux_file(Output_t const * const o /**< Output object.*/
                    );


/** @brief Write a column of fluxes to the output file.*/
void write_fluxes(Output_t const * const o, /**< Output object.*/
                  int const index, /**< Variable index.*/
                  int const column, /**< Column index.*/
                  fp_t const * const flux /**< Fluxes [W/m^2] (level).*/
                 );


#endif
