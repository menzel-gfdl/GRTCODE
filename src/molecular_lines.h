/** @file */
#ifndef MOLECULAR_LINES_H_
#define MOLECULAR_LINES_H_

#include <stdint.h>
#include "floating_point_type.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "spectral_bin.h"
#include "water_vapor_continuum.h"


/*Macros.*/
#define DIR_PATH_LEN 1024
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif


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
    @brief Library context.
*/
typedef struct GrtContext
{
    /*--- Parameters directly supplied by the user. ---*/
    int gpu_id; /**< Id of the GPU that is associated with this context.*/
    int num_threads; /**< Number of CPU threads that will be used to calculate
                          the lines (if not using a GPU).*/
    int num_levels; /**< Number of atmospheric levels.*/
    double w0; /**< First point of the spectral grid [1/cm].*/
    double wn; /**< Last point of the spectral grid [1/cm].*/
    double wres; /**< Spectral grid resolution [1/cm].*/
    double wcutoff; /**< Cutoff from spectral line center [1/cm].*/
    char hitran_path[DIR_PATH_LEN]; /**< Path to the HITRAN database file.*/
    char h2o_ctm_dir[DIR_PATH_LEN]; /**< Path to the water vapor continuum
                                         directory.*/
    char o3_ctm_dir[DIR_PATH_LEN]; /**< Path to the ozone continuum
                                        directory.*/
    int optical_depth_method; /**< Flag specifying which method will be
                                   used to calculate the optical depths.*/

    /*--- Parameters implicitly defined by the code. ---*/
    int num_layers; /**< Number of atmospheric layers (num_levels-1).*/
    int num_molecules; /**< Number of molecules.*/
    uint64_t molecule_bit_field; /**< Bit field used to determine which
                                      molecules are currently in use.*/
    Molecule_t mols[NUM_MOLS]; /**< Array of molecule structures.*/
    uint64_t num_wpoints; /**< Number of spectral grid points.*/
    int use_h2o_ctm; /**< Flag indicating if using the water vapor
                          continuum is used.*/
    int use_o3_ctm; /**< Flag indicating if the ozone continuum is used.*/
    SpectralBins_t bins; /**< Spectral bins.*/

    fp_t *x; /**< Abundance (molecules,levels).*/
    fp_t *n; /**< Total number of molecules [1/cm^2] (layers).*/
    fp_t *pavg; /**< Pressure [atm] (layers).*/
    fp_t *tavg; /**< Temperature [K] (layers).*/
    fp_t *psavg; /**< Molecular partial pressure [atm] (layers).*/
    fp_t *ns; /**< Number of molecules [1/cm^2] of a particular species
                   (layers).*/
    fp_t *linecenter; /**< Pressure-shifted line center positions [1/cm]
                           (layers,lines).*/
    fp_t *snn; /**< Line strength [1/cm] (layers,lines).*/
    fp_t *gamma; /**< Temperature- and pressure-corrected lorentz
                       half-width [1/cm] (layers,lines).*/
    fp_t *alpha; /**< Doppler half-width [1/cm] (layers,lines).*/
    WaterVaporContinuumCoefs_t h2o_cc; /**< Structure containing water vapor
                                            continuum coefficients.*/
    OzoneContinuumCoefs_t o3_cc; /**< Structure containing ozone continuum
                                      coefficients.*/
    fp_t *p; /**< Pressure [atm] (levels).*/
    fp_t *t; /**< Temperatture [K] (levels).*/
    fp_t *tau; /**< Optical depths (layer,wavenumber).*/


#ifdef FOO
    LineParams_t *lines; /**< Molecular line parameters.*/
#endif
} GrtContext_t;


/** @ingroup capi
    @brief Flags used to specifiy which method is used to calculate the
           optical depths.*/
enum OpticalDepthMethod
{
    wavenumber_sweep,
    line_sweep,
    line_sample
};


/**
    @ingroup capi
    @brief Initialize a context.
    @return 0 if completed successfully, or else an error code.
*/
EXTERN int grt_context_init(GrtContext_t **context, /**< Library context.*/
                            int const num_levels, /**< Number of atmospheric levels.*/
                            double const w0, /**< Lowest wavenumber [1/cm] on spectral grid.*/
                            double const wn, /**< Highest wavenumber [1/cm] on spectral grid.*/
                            double const wres, /**< Spectral grid resolution [1/cm].*/
                            char const * const hitran_path, /**< Path to the HITRAN
                                                                 database file.*/
                            char const * const h2o_ctm_dir, /**< Directory containing the
                                                                 provided water vapor continuum
                                                                 input files.  If NULL, then
                                                                 the water vapor continuum
                                                                 is not included in the optical
                                                                 depth calculation.*/
                            char const * const o3_ctm_dir, /**< Directory containing the
                                                                provided ozone continuum
                                                                input files.  If NULL, then
                                                                the ozone continuum is not
                                                                included in the optical depth
                                                                calculation.*/
                            double const * const wcutoff, /**< Cutoff [1/cm] from spectral
                                                               line center.  If NULL, this
                                                               defaults to 25 [1/cm].*/
                            int const * const gpu_id, /**< Id of the GPU that will be associated
                                                           with this context.  If NULL, then
                                                           use GPU 0 if at least one GPU
                                                           exists on the system, or else
                                                           set to -1 (corresponding to
                                                           a host only run.*/
                            int const * const num_threads, /**< If running on the host CPU,
                                                                determines the maximum number
                                                                of OpenMP threads that will
                                                                be used.  If NULL, default to
                                                                omp_get_max_threads (or one
                                                                if not build with OpenMP).*/
                            int const * const optical_depth_method /**< Flag specifying which
                                                                        method will be used to
                                                                        calculate the optical
                                                                        depths.  Defaults to
                                                                        wavenumber_sweep.*/
                           );


/**
    @ingroup capi
    @brief Release memory allocated by the context.
    @return 0 if completed successfully, or else an error code.
*/
EXTERN int grt_context_free(GrtContext_t **context /**< Library context.*/
                           );


/**
    @ingroup capi
    @brief Add a molecule.  The optical depths of all added molecules
           will be computed and summed to give the total optical optical
           depth of each atmospheric layer at each spectral grid point.
    @return 0 if completed successfully, or else an error code.
*/
EXTERN int grt_add_molecule(GrtContext_t *context, /**< Library context.*/
                            int const molecule_id, /**< Id that is associated with the molecule. */
                            double const * const min_line_center, /**< Lower bound [1/cm] of spectral range.
                                                                       Only lines with line center wavenumbers
                                                                       greater than or equal to this will be
                                                                       computed.  If the value NULL is passed
                                                                       in, this defaults to 1 [1/cm].*/
                            double const * const max_line_center /**< Upper bound [1/cm] of spectral range.
                                                                      Only lines with line center wavenumbers
                                                                      less than or equal to this will be
                                                                      computed.  If the value NULL is passed
                                                                      in, this defaults to 3250 [1/cm].*/
                           );


/**
    @ingroup capi
    @brief Update the abundances [ppmv] for a molecule.
    @return 0 if completed successfully, or else an error code.
*/
EXTERN int grt_set_molecule_ppmv(GrtContext_t *context, /**< Library context.*/
                                 int const molecule_id, /**< Molecule id returned by
                                                             @ref add_molecule.*/
                                 fp_t const * const ppmv /**< Array of molecular abundances [ppmv].
                                                              The size of this array must be
                                                              equal to the number of atmospheric
                                                              levels.*/
                                );


/**
    @ingroup capi
    @brief Calculate the total optical depth of each atmospheric layer at
           each spectral grid point.
    @return 0 if completed successfully, or else an error code.
*/
EXTERN int grt_calculate_optical_depth(GrtContext_t *context, /**< Library context.*/
                                       fp_t const * const pressure, /**< Array of atmospheric pressures [mb].
                                                                         The size of this array must be
                                                                         equal to the number of atmospheric
                                                                         levels.*/
                                       fp_t const * const temperature, /**< Array of atmospheric temperatures [K].
                                                                            The size of this array must be
                                                                            equal to the number of atmospheric
                                                                            levels.*/
                                       fp_t *optical_depth /**< Array of atmospheric optical depths.
                                                                The size of this array must be equal
                                                                to the number of atmospheric layers
                                                                times the number of spectral grid
                                                                points.  Memory is layed out as
                                                                (layer,wavenumber) (i.e., the
                                                                fastest changing dimension is the
                                                                one corresponding to the spectral
                                                                grid.)*/
                                      );


/**
    @ingroup capi
    @brief Get the number of molecules that have been added to the context.
    @return 0 if completed successfully, or else an error code.
*/
EXTERN int grt_get_num_molecules(GrtContext_t const * const context, /**< Library context.*/
                                 int * const n /**< Number of molecules.*/
                                );


/**
    @ingroup capi
    @brief Get the number of spectral grid points for the input context.
    @return 0 if completed successfully, or else an error code.
*/
EXTERN int grt_get_spectral_grid_size(GrtContext_t const * const context, /**< Library context.*/
                                      uint64_t * const n /**< Spectral grid size.*/
                                     );


/**
    @ingroup capi
    @brief Return a message describing the input return code.
    @return 0 if completed successfully, or else an error code.
*/
EXTERN int grt_errstr(int const code, /**< Code returned from one of the
                                           GRT functions.*/
                      char * const buf, /**< Buffer where message will be stored.*/
                      int const buf_size /**< Size of the input message buffer.*/
                     );


/**
    @ingroup capi
    @brief Set the verbosity level for the library.
*/
EXTERN void grt_set_verbosity(int const level /**< Verbosity level.  Levels range
                                                   from 0 (least verbose) to 3
                                                   (most verbose).  The default
                                                   level is 0.*/
                             );


/**
    @ingroup capi
    @brief Get the verbosity level for the library.
    @return Current verbosity level.
*/
EXTERN int grt_get_verbosity(void);


#endif
