/** @file */
#ifndef MOLECULAR_LINES_H_
#define MOLECULAR_LINES_H_

#include <stdint.h>
#include "cfcs.h"
#include "collision_induced_absorption.h"
#include "device.h"
#include "extern.h"
#include "floating_point_type.h"
#include "molecules.h"
#include "optics.h"
#include "ozone_continuum.h"
#include "spectral_bin.h"
#include "spectral_grid.h"
#include "water_vapor_continuum.h"


/*Macros.*/
#define DIR_PATH_LEN 1024


/**
    @defgroup capi C API
    @section Overview
        Given an atmospheric column made up of at least one layer,
        this code calculates the total optical depth of each layer
        at each point on an input spectral grid.  To use this API,

            #include "molecular_lines.h"

        and follow these steps:
            -# Declare library context pointer(s).  These variables must be
               pointers of type GrtContext_t, as in:

                   GrtContext_t *x;

               Each context pointer will
               hold the address of a struct that contains data that is
               required by the library.
               @note Each context can be associated with a single GPU.  If
                   you wish to use multiple GPUs, you must create a context
                   pointer for each GPU you wish to use and pass
                   in the appropriate device id when initializing the context
                   pointer (i.e., by calling the @ref grt_context_init
                   function).

            -# Initialize the context pointer(s) by calling the
               @ref grt_context_init function.  Here you must provide the
               number of levels per atmospheric column and parameters
               (lower bound, upper bound, and resolution) for the
               spectral grid on which the optical depths will be calculated.
               @attention You must call this function before calling any
                   other function included in this library.  Failure to do
                   so will result in undefined behavior.
               @note Each atmospheric level corresponds to an interface
                   between adjacent atmospheric layers or the lower/upper
                   edge of the atmosphere.  Thus, the number of atmospheric
                   levels = the number of atmospheric layers plus one.
                   Since at least one atmospheric layer is required,
                   the number of atmospheric levels must be greater than or
                   equal to two.

               In addition, the following parameters may be set (passing
               in NULL implies using the default values):
               - A cutoff value for the molecular lines.  This
                 value denotes how far (in terms of wavenumber) from the line
                 center each molecular line is calculated out to.  By default
                 this value is set to 25 [1/cm], so if for example a line
                 center lies at the wavenumber 150 [1/cm], then it contributes
                 to the optical depth at all spectral grid points in the range
                 125 <= w <= 175 [1/cm].
               - The id of the GPU device you wish to associate with this
                 context (see the readme for a simple way to determine the ids of
                 any GPUs you may have on your system).  If you do not specify a
                 GPU id when initializing the context (i.e., by passing in NULL
                 for this argument), the library will query the system for available
                 GPUs.  If any are found, then the first device (device 0) will be used.
                 If you wish to force the code to run on your host CPU, pass in a value
                 of -1.
               - If you wish to run with either the water vapor
                 or ozone continua, you must provide a path to a directory
                 containing the necessary input files.  The required input
                 files are included with this library in directories named
                 "water_vapor_continuum" and "ozone_continuum" respectively,
                 and are located in the base directory of this repository.
                 If you wish to run without either continuum, instead pass
                 in the value NULL.
            -# Add each molecule that you want included in the optical
               depth calculation by calling the @ref grt_add_molecule function.
               Each added molecule requires an ascii [HITRAN](http://hitran.org)
               database file containing the necessary molecular line parameters.  The
               format of these files must match that described in Table 1 of
               [Rothman et al. 2013, Journal of Quantitative Spectroscopy & Radiative
               Transfer, 130](http://dx.doi.org/10.1016/j.jqsrt.2013.07.002).
               Example HITRAN database files for a select set of molecules are included
               with this library in a directory labeled HITRAN_files in the base of this
               repository.
               @attention Ozone and water vapor continua will only be
                   included in the optical depth calculation if the ozone
                   and water vapor molecules are added.

            -# Set the abundance [ppmv] of each added molecule by calling
               the @ref grt_set_molecule_ppmv function.
               @attention The input abundance array must be contiguous and
                   its size (number of elements) must be equal to the number
                   of atmospheric levels passed into the @ref grt_context_init
                   function, or else the behavior is undefined.

            -# Calculate the optical depth for each layer in the column at each
               spectral grid point by calling the @ref grt_calculate_optical_depth
               function.
               @attention All input arrays must be contiguous.  In addition,
                   the number of elements in the input pressure [mb] and
                   temperature [K] arrays must be equal to the number of
                   atmospheric levels.  The number of elements in the input
                   optical depth array must be equal to the number of
                   atmospheric layers (number of atmospheric levels minus one)
                   times the number of spectral grid points (returned by
                   the @ref grt_get_spectral_grid_size function).  If any of
                   these arrays are not contiguous or have an incorrect size,
                   the behavior is undefined.

            -# Release the memory allocated by the context(s) by calling the
               @ref grt_context_free function.
    @section Example
    Here is a simple example demonstrating how to use this library.
    @include example.c
    In order to build this code, copy this code into a file, modify
    the paths to input files (as needed), and (assuming you have gcc
    installed) run:

        $ gcc <file> -fopenmp -o example.x -I<path to library include directory> \
              -L<path to library lib directory> -lmolecular_lines \
              -Wl,-rpath=<path to library lib directory>

    To run this example on your GPU, make sure that you have compiled
    the library using the NVCC compiler (i.e., by using the provided
    Makefile.nvcc).
*/


/** @ingroup capi
    @brief Molecular lines object.
*/
typedef struct MolecularLines
{
    Device_t device; /**< Id of the device associated with this object.*/
    int num_levels; /**< Number of atmospheric levels.*/
    int num_layers; /**< Number of atmospheric layers (num_levels-1).*/

    int num_molecules; /**< Number of molecules.*/
    uint64_t molecule_bit_field; /**< Bit field used to determine which molecules are currently in use.*/
    Molecule_t mols[NUM_MOLS]; /**< Array of molecule structures.*/

    int num_cfcs; /**< Number of cfcs.*/
    uint64_t cfc_bit_field; /**< Bit field used to determine which cfcs are currently in use.*/
    CfcCrossSection_t cfcs[NUM_CFCS]; /**< CFC cross section data structures.*/
    fp_t *x_cfc; /**< CFC abundance (CFC, levels).*/

    int num_cias; /**< Number of collision-induced absorption continua.*/
    uint64_t cia_bit_field; /**< Bit field used to determine with species are currently in use.*/
    CollisionInducedAbsorption_t cia[MAX_NUM_CIAS]; /**< Collision-induce absorption structures.*/
    fp_t *x_cia; /**< Collision-induced absorption abundances (molecule, level).*/

    char h2o_ctm_dir[DIR_PATH_LEN]; /**< Path to the water vapor continuum directory.*/
    int use_h2o_ctm; /**< Flag indicating if using the water vapor continuum is used.*/
    WaterVaporContinuumCoefs_t h2o_cc; /**< Water vapor continuum coefficients.*/

    char o3_ctm_dir[DIR_PATH_LEN]; /**< Path to the ozone continuum directory.*/
    int use_o3_ctm; /**< Flag indicating if the ozone continuum is used.*/
    OzoneContinuumCoefs_t o3_cc; /**< Oone continuum coefficients.*/

    SpectralGrid_t grid; /**< Spectral grid.*/
    SpectralBins_t bins; /**< Spectral bins.*/
    char hitran_path[DIR_PATH_LEN]; /**< Path to the HITRAN database file.*/
    double wcutoff; /**< Cutoff from spectral line center [1/cm].*/
    int optical_depth_method; /**< Flag specifying which method will be used to calculate the optical depths.*/
    fp_t *x; /**< Abundance (molecule, level).*/
    fp_t *n; /**< Integrated number density [1/cm^2] (layer).*/
    fp_t *pavg; /**< Pressure [atm] (layer).*/
    fp_t *tavg; /**< Temperature [K] (layer).*/
    fp_t *psavg; /**< Molecular partial pressure [atm] (layer).*/
    fp_t *ns; /**< Integrated number density [1/cm^2] of a particular species (layer).*/
    fp_t *linecenter; /**< Pressure-shifted line center position [1/cm] (layer, line).*/
    fp_t *snn; /**< Line strength [1/cm] (layer, line).*/
    fp_t *gamma; /**< Temperature- and pressure-corrected lorentz half-width [1/cm] (layer, line).*/
    fp_t *alpha; /**< Doppler half-width [1/cm] (layer, line).*/
    fp_t *p; /**< Pressure [atm] (level).*/
    fp_t *t; /**< Temperature [K] (level).*/
    fp_t *tau; /**< Optical depth (layer, wavenumber).*/
} MolecularLines_t;


/** @ingroup capi
    @brief Flags used to specifiy which method is used to calculate the
           optical depths.*/
enum OpticalDepthMethod
{
    wavenumber_sweep,
    line_sweep,
    line_sample
};


/** @brief Reserve memory for molecular lines.
    @return RS_SUCCESS or an error code.*/
EXTERN int create_molecular_lines(MolecularLines_t * const ml, /**< Molecular lines object.*/
                                  int const num_levels, /**< Number of atmospheric levels.*/
                                  SpectralGrid_t const * const grid, /**< Spectral grid.*/
                                  Device_t const * const device, /**< Device.*/
                                  char const * const hitran_path, /**< Path to HITRAN database file.*/
                                  char const * const h2o_ctm_dir, /**< Path to water vapor continuum directory.*/
                                  char const * const o3_ctm_dir, /**< Path to ozone continuum directory.*/
                                  double const * const wcutoff, /**< Cutoff from line center [1/cm].*/
                                  int const * const optical_depth_method /**< Method to use to calculate optical depths.*/
                                 );


/** @brief Free memory for the molecular lines.
    @return RS_SUCCESS or an error code.*/
EXTERN int destroy_molecular_lines(MolecularLines_t * const ml /**< Molecular lines object.*/
                                  );


/** @brief Add a molecule.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_add_molecule(MolecularLines_t * const ml, /**< Molecular lines object.*/
                            int const molecule_id, /**< Molecule id.*/
                            double const * const min_line_center, /**< Lower bound [1/cm]
                                                                       for spectral line
                                                                       centers.*/
                            double const * const max_line_center /**< Upper bound [1/cm]
                                                                      for spectral line
                                                                      centers.*/
                           );


/** @brief Update a molecule's ppmv.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_set_molecule_ppmv(MolecularLines_t * const ml, /**< Molecular lines object.*/
                                 int const molecule_id, /**< Molecule id.*/
                                 fp_t const * const ppmv /**< Abundance [ppmv] (level).*/
                                );


/** @brief Add a CFC.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_add_cfc(MolecularLines_t * const ml, /**< Molecular lines object.*/
                       int const cfc_id, /**< CFC id.*/
                       char const * const filepath /**< Path to CFC cross section csv file.*/
                      );


/** @brief Update a CFC's ppmv.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_set_cfc_ppmv(MolecularLines_t * const ml, /**< Molecular lines object.*/
                            int const cfc_id, /**< CFC id.*/
                            fp_t const * const ppmv /**< Abundance [ppmv] (level).*/
                           );


/** @brief Activate collision-induced absorption between two species.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_add_cia(MolecularLines_t * const ml, /**< Molecular lines object.*/
                       int const species1, /**< Id of species.*/
                       int const species2, /**< Id of species.*/
                       char const * const filepath /**< Path to cross section csv file.*/
                      );


/** @brief Update a CIA species' ppmv.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_set_cia_ppmv(MolecularLines_t * const ml, /**< Molecularlines object.*/
                            int const cia_id, /**< CIA species id.*/
                            fp_t const * const ppmv /**< Abundance [ppmv] (level).*/
                           );


/** @brief Calcluate the total optical depth in each layer at each spectral grid point.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_calculate_optical_depth(MolecularLines_t * const ml, /**< Molecular lines object.*/
                                       fp_t * const pressure, /**< Pressure [mb] (level).*/
                                       fp_t * const temperature, /**< Temperature [K] (level).*/
                                       Optics_t * const optics /**< Optics object.*/
                                      );


/** @brief Get the number of molecules.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_get_num_molecules(MolecularLines_t const * const ml, /**< Molecular lines object.*/
                                 int * const n /**< Number of molecules.*/
                                );


/** @brief Return a message for an input return code.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_errstr(int const code, /**< Error code.*/
                      char * const buf, /**< Buffer to hold error message.*/
                      int const buf_size /**< Size of input buffer.*/
                     );


#endif
