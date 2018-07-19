/** @file */
#ifndef MOLECULAR_LINES_H_
#define MOLECULAR_LINES_H_

#include <stdint.h>
#include "floating_point_type.h"
#include "ozone_continuum.h"
#include "parse_HITRAN_file.h"
#include "water_vapor_continuum.h"


/**
    \defgroup capi C API
    \brief foobar
    \section Overview
        Given an atmospheric column made up of at least one layer,
        this code calculates the total optical depth of each layer
        at each point on an input spectral grid.  Typical usage of this
        library includes the following steps:\n\n
            -# Declare library context pointer(s).  Each context pointer will
               hold the address of a struct that contains data that is
               required by the library.\n\n
               \note Each context can be associated with a single GPU.  If
                   you wish to use multiple GPUs, you must create a context
                   pointer for each GPU you wish to use and pass
                   in the appropriate device id when initalizing the context
                   pointer (i.e., by calling the @ref grt_context_init
                   function).

               \n
            -# Initialize the context pointer(s) by calling the
               @ref grt_context_init function.  Here you must provide the
               number of levels per atmospheric column and parameters
               (lower bound, upper bound, and resolution) for the
               spectral grid on which the optical depths will be calculated.
               \n\n
               \attention You must call this function before calling any
                   other function included in this library.  Failure to do
                   so will result in undefined behavior.

               \n
               \note Each atmospheric level corresponds to an interface
                   between adjacent atmospheric layers or the lower/upper
                   edge of the atmosphere.  Thus, the number of atmospheric
                   levels = the number of atmospheric layers plus one.
                   Since at least one atmospheric layer is required,
                   the number of atmospheric levels must be greater than or
                   equal to two.

               \n In addition, the following parameters may be set (passing
               in NULL implies using the default values):\n\n
               - A cutoff value for the molecular lines.  This
                 value denotes how far (in terms of wavenumber) from the line
                 center each molecular line is calculated out to.  By default
                 this value is set to 25 [1/cm], so if for example a line
                 center lies at the wavenumber 150 [1/cm], then it contributes
                 to the optical depth at all spectral grid points in the range
                 125 <= w <= 175 [1/cm].\n\n
               - The id of the GPU device you wish to associate with this
                 context.  A list of NVIDIA GPUs on your system can be
                 found by running:\n\n
                 ```$ nvidia-smi --list-gpus```\n\n
                 If you do not specify a GPU id when initializing the context,
                 the library will query the system for available GPUs.  If any
                 are found, then the first device (device 0) will be used.
                 If you wish to run only on your host CPU, pass in a value
                 of -1.\n\n
                 \note If you wish to run on a GPU, you must have a NVIDIA GPU,
                     install CUDA, and compile with the NVCC compiler (see
                     Requirements below).

                 \n
               - If you wish to run with either the water vapor
                 or ozone continua, you must provide a path to a directory
                 containing the necessary input files.  The required input
                 files are included with this library in directories named
                 "water_vapor_continuum" and "ozone_continuum" respectively,
                 and are located in the base directory of this repository.
                 If you wish to run without either continuum, instead pass
                 in the value NULL.\n\n
            -# Add each molecule that you want included in the optical
               depth calculation by calling the @ref grt_add_molecule function.
               Each added molecule requires an ascii [HITRAN]
               (http://hitran.org) database file
               containing the necessary molecular line parameters.  The
               format of these files must match that described in Table 1 of
               [Rothman et al. 2013, Journal of Quantitative Spectroscopy
               & Radiative Transfer, 130]
               (http://dx.doi.org/10.1016/j.jqsrt.2013.07.002).
               Example HITRAN database files
               for a select set of molecules are included with this library in
               a directory labeled HITRAN_files in the base of this
               repository.\n\n
               \attention Ozone and water vapor continua will only be
                   included in the optical depth calculation if the ozone
                   and water vapor molecules are added.

               \n
            -# Set the abundance [ppmv] of each added molecule by calling
               the @ref grt_set_molecule_ppmv function.\n\n
               \attention The input abundance array must be contiguous and
                   its size (number of elements) must be equal to the number
                   of atmospheric levels passed into the @ref grt_context_init
                   function, or else the behavior is undefined.

               \n
            -# Calculate the optical depth for each layer in the column
               at each spectral grid point by calling the
               @ref grt_calculate_optical_depth function.\n\n
               \attention All input arrays must be contiguous.  In addition,
                   the number of elements in the input pressure [atm] and
                   temperature [K] arrays must be equal to the number of
                   atmospheric levels.  The number of elements in the input
                   optical depth array must be equal to the number of
                   atmospheric layers (number of atmospheric levels minus one)
                   times the number of spectral grid points (returned by
                   the @ref grt_context_init function).  If any of
                   these arrays are not contiguous or have an incorrect size,
                   the behavior is undefined.

               \n
            -# Release the memory allocated by the context(s) by calling the
               @ref finalize_grt function.\n\n

    \section Example
        Here is a simple example demonstrating how to use this library.
        \include example.c
        In order to build this code, copy this code into a file and
        (assuming you have gcc installed), run:\n\n
        ```gcc <file> -o example.x -I<path to library include directory>
           -L<path to library lib directory> -lmolecular_lines```\n\n
        To run this example on your GPU, make sure that you have compiled
        the library using the NVCC compiler (i.e., by using the provided
        Makefile.nvcc).
*/


/** @ingroup capi
    @brief Library context.*/
typedef struct GrtContext
{
    int num_levels; /**< Number of atmospheric levels.*/
    int num_layers; /**< Number of atmospheric layers (num_levels-1).*/
    int num_molecules; /**< Number of molecules.*/
    LineParams_t **line_params; /**< Array of structures containg molecular
                                     line parameters (from HITRAN database
                                     files.)*/
    double w0; /**< First point of the spectral grid [1/cm].*/
    double wn; /**< Last point of the spectral grid [1/cm].*/
    double wres; /**< Spectral resolution [1/cm].*/
    uint64_t num_wpoints; /**< Number of spectral grid points.*/
    double wcutoff; /**< Cutoff from spectral line center [1/cm].*/
    int gpu_id; /**< Id of the GPU that is associated with this context.*/
    int num_threads; /**< Number of CPU threads that will be used to calculate
                          the lines (if not using a GPU).*/
    int use_h2o_ctm; /**< Flag for using the water vapor continuum.*/
    WaterVaporContinuumCoefs_t *h2o_cc; /**< Structure containing water vapor
                                             continuum coefficients.*/
    int use_o3_ctm; /**< Flag for using the ozone continuum.*/
    OzoneContinuumCoefs_t *o3_cc; /**< Structure containing ozone continuum
                                       coefficients.*/
    fp_t *P; /**< Pressure [atm] at each level.*/
    fp_t *T; /**< Temperatture [K] at each level.*/
    fp_t *x; /**< Molecular abundance at each level.*/
    fp_t *Pavg; /**< Pressure [atm] in each layer.*/
    fp_t *Tavg; /**< Temperature [K] in each layer.*/
    fp_t *N; /**< Total number of molecules [1/cm^2] in each layer.*/
    fp_t *Ns; /**< Number of molecules [1/cm^2] (of a particular species) in
                   each layer.*/
    fp_t *Psavg; /**< Molecular partial pressure [atm] in each layer.*/
    fp_t *snn_ref; /**< */
    fp_t *gamma; /**< */
    fp_t *Pshift; /**< */
    fp_t *s; /**< */
    fp_t *tau; /**< Optical depths (layer,wavenumber).*/
    LineParams_t *lines; /**< Molecular line parameters.*/
} GrtContext_t;


/**
    @ingroup capi
    @brief Initialize a context.
    @return 0 if completed successfully, or else an error code.
*/
#ifdef __NVCC__
extern "C"
#endif
int grt_context_init(GrtContext_t **context, /**< Library context.*/
                     int const num_levels, /**< Number of atmospheric levels.*/
                     double const w0, /**< Lowest wavenumber [1/cm] on spectral grid.*/
                     double const wn, /**< Highest wavenumber [1/cm] on spectral grid.*/
                     double const wres, /**< Spectral grid resolution [1/cm].*/
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
                                                         be used.  Defaults to
                                                         omp_get_max_threads (or one
                                                         if not build with OpenMP).*/
                     char const * const h2o_ctm_dir, /**< Directory containing the
                                                          provided water vapor continuum
                                                          input files.  If NULL, then
                                                          the water vapor continuum
                                                          is not included in the optical
                                                          depth calculation.*/
                     char const * const o3_ctm_dir /**< Directory containing the
                                                        provided ozone continuum
                                                        input files.  If NULL, then
                                                        the ozone continuum is not
                                                        included in the optical depth
                                                        calculation.*/
                    );


/**
    @ingroup capi
    @brief Release memory allocated by the context.
    @return 0 if completed successfully, or else an error code.
*/
#ifdef __NVCC__
extern "C"
#endif
int grt_context_free(GrtContext_t **context /**< Library context.*/
                    );


/**
    @ingroup capi
    @brief Add a molecule.  The optical depths of all added molecules
           will be computed and summed to give the total optical optical
           depth of each atmospheric layer at each spectral grid point.
    @return 0 if completed successfully, or else an error code.
*/
#ifdef __NVCC__
extern "C"
#endif
int grt_add_molecule(GrtContext_t *context, /**< Library context.*/
                     char const * const hitran_filepath, /**< Path to HITRAN ascii file containing
                                                              molecular line parameters.*/
                     int * const molecule_id, /**< Id that is associated with the molecule. */
                     double const * const min_line_center_wavenumber, /**< Lower bound [1/cm] of spectral range.
                                                                           Only lines with line center wavenumbers
                                                                           greater than or eqaul to this will be
                                                                           computed.  Defaults to 1 [1/cm].*/
                     double const * const max_line_center_wavenumber /**< Upper bound [1/cm] of spectral range.
                                                                          Only lines with line center wavenumbers
                                                                          less than or equal to this will be
                                                                          computed.  Defaults to 50,000 [1/cm].*/
                    );


/**
    @ingroup capi
    @brief Update the abundances [ppmv] for a molecule.
    @return 0 if completed successfully, or else an error code.
*/
#ifdef __NVCC__
extern "C"
#endif
int grt_set_molecule_ppmv(GrtContext_t *context, /**< Library context.*/
                          int const molecule_id, /**< Molecule id returned by
                                                      @ref add_molecule.*/
                          fp_t const * const ppmv /**< Array of molecular abundances [ppmv].
                                                       The size of this array must be
                                                       eqaul to the number of atmospheric
                                                       levels.*/
                         );


/**
    @ingroup capi
    @brief Calculate the total optical depth of each atmospheric layer at
           each spectral grid point.
    @return 0 if completed successfully, or else an error code.
*/
#ifdef __NVCC__
extern "C"
#endif
int grt_calculate_optical_depth(GrtContext_t *context, /**< Library context.*/
                                fp_t const * const pressure, /**< Array of atmospheric pressures [atm].
                                                                  The size of this array must be
                                                                  eqaul to the number of atmospheric
                                                                  levels.*/
                                fp_t const * const temperature, /**< Array of atmospheric temperatures [K].
                                                                     The size of this array must be
                                                                     eqaul to the number of atmospheric
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
#ifdef __NVCC__
extern "C"
#endif
int grt_get_num_molecules(GrtContext_t const * const context, /**< Library context.*/
                          int * const n /**< Number of molecules.*/
                         );


/**
    @ingroup capi
    @brief Get the number of spectral grid points for the input context.
    @return 0 if completed successfully, or else an error code.
*/
#ifdef __NVCC__
extern "C"
#endif
int grt_get_spectral_grid_size(GrtContext_t const * const context, /**< Library context.*/
                               uint64_t * const n /**< Spectral grid size.*/
                              );


/**
    @ingroup capi
    @brief Return a message describing the input return code.
    @return 0 if completed successfully, or else an error code.
*/
#ifdef __NVCC__
extern "C"
#endif
int grt_errstr(int const code, /**< Code returned from one of the
                                    GRT functions.*/
               char * const buf, /**< Buffer where message will be stored.*/
               int const buf_size /**< Size of the input message buffer.*/
              );


#endif
