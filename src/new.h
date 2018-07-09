/** @file */
#ifndef NEW_H_
#define NEW_H_

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
            -# Declare a library context pointer.  This pointer will hold
               the address of a struct that contains all data that is
               stored/used by the library.\n\n
            -# Set the number of atmospheric levels per column.  Each
               atmospheric level corresponds to an interface  between
               adjacent atmospheric layers.  Thus, the number of atmospheric
               levels = the number of atmospheric layers plus one.  Since
               at least one atmospheric layer is required, the number of
               atmospheric levels must be greater than or equal to two.\n\n
            -# Set the upper and lower bounds and resolution of the spectral
               grid on which the optical depths will be calculated.\n\n
            -# Optional: Set a cutoff value for the molecular lines.  This
               value denotes how far (in terms of wavenumber) from the line
               center each molecular line is calculated out to.  By default
               this value is set to 25 [1/cm], so if for example a line
               center lies at the wavenumber 150 [1/cm], then it contributes
               to the optical depth at all spectral grid points in the range
               125 <= w <= 175 [1/cm].\n\n
            -# Optional: Set a flag denoting the architecture you wish to
               run on.  Here a value of 0 corresponds to running on your host
               CPU, while any other value corresponds to running on your GPU.
               If this flag is not set, the library will attempt to use your
               default GPU (if you have one and you build with the NVIDIA
               NVCC compiler, see Requirements below).  If no NVIDIA GPU is
               detected, then only the host CPU will be used.  Please note
               that if you wish to run on your GPU, it must be an NVIDIA GPU
               and you must have CUDA installed and compile with the NVCC
               compiler.\n\n
            -# Optional: Set flags determining whether or not the continuum
               for water vapor and ozone will be included in the optical
               depth calculation.  A value of 0 corresponds to running
               without the specified continuum, while another value
               corresponds to running with the continuum.  Please note that
               running with the continua requires specific input files.
               By default the library will attempt to include the continua
               if the required inputs are found, otherwise they will not
               be included in the optical depth calculation.\n\n
            -# Initialize the library by calling the @ref initialize_grt
               function with the necessary arguments.\n\n
            -# Add each molecule that you want included in the optical
               depth calculation by calling the @ref add_molecule function.
               Each added molecule requires an ascii [HITRAN]
               (http://hitran.org) database file
               containing the necessary molecular line parameters.  The
               format of these files must match that described in Table 1 of
               [Rothman et al. 2013, Journal of Quantitative Spectroscopy
               & Radiative Transfer, 130]
               (http://dx.doi.org/10.1016/j.jqsrt.2013.07.002).
               Example HITRAN database files
               for a select set of molecules is included with this library in
               a directory labeled HITRAN_FILES in the base of this
               repository. Please note that the ozone and water vapor
               continua will only be included if the ozone and water vapor
               molecules are added.\n\n
            -# Set the abundance [ppmv] of each added molecule by calling
               the @ref set_molecule_ppmv function.  The number of elements
               in the input abundance array must be
               equal to the number of atmospheric levels passed into the
               @ref initialize_grt function, or else the behavior is
               undefined.\n\n
            -# Calculate the optical depth for each layer in the column
               at each spectral grid point by calling the
               @ref calculate_optical_depth function.  The number of elements
               in the input pressure [atm] and temperature [K] arrays must
               be equal to the number of atmospheric levels.  The number
               of elements in the input optical depth array must be equal
               to the number of atmospheric layers (number of atmospheric
               levels minus one) times the number of spectral grid points
               (returned by the @ref initialize_grt function).  If any of
               these arrays have an incorrect size, the behavior is
               undefined.\n\n
            -# Finalize the library by calling the @ref finalize_grt function.
               This function frees all memory allocated by the library.\n\n

    \section Requirements
        This library requires a c and fortran compiler, such as the freely
        available gcc and gfortran.  The c compiler must support the c99
        standard, and the fortran compiler must support this [Fortran 2012
        Technical Specification] (https://www.iso.org/standard/45136.html).
        In order to run on a NVIDIA GPU,
        CUDA must also be installed and the library must be build using the
        included NVCC compiler.  This library also requires make and install.

    \section Code
        The source code currently resides in this [Gitlab repository]
        (https://gitlab.gfdl.noaa.gov/Raymond.Menzel/grtcodev2), on branch
        modular_lines.  To obtain the code, run:\n\n
        ```$ git clone
             https://gitlab.gfdl.noaa.gov/Raymond.Menzel/grtcodev2.git```\n
        ```$ git checkout modular_lines```

    \section Building
        The library provides two Makefiles, one called Makefile for building
        for CPU only runs and one call Makefile.nvcc for building for GPU
        runs.  The CPU only makefile assumes gcc and gfortran as the default
        compilers, but those can be overwritten by setting CC= your c compiler
        and FC= your fortran compiler when running make.  Note that the
        CFLAGS= and FFLAGS= will most likely also need to be overridden.
        For example, to build for CPU only runs using the default settings,
        simply run:\n\n
        ```$ make```\n\n
        or to build with the intel compilers, run:\n\n
        ```$ make CC=icc CFLAGS=-O3 FC=ifort FFLAGS=-O3```\n\n
        To build for GPU runs with CUDA installed and the NVCC compiler
        in your path, simply run:\n\n
        ```$ make -f Makefile.nvcc```\n\n
        After building the library, run: \n\n
        ```$ make test```\n\n and/or\n\n ```$ make -f Makefile.nvcc test```\n\n
        to run some tests to make sure everything is working.  Lastly,
        run:\n\n
        ```make install```\n\n or
        ```make -f Makefile.nvcc install```\n
        to install the library.  PREFIX= can be used to select the directory
        where the library will be installed, or else the library will be
        installed in the current directory.

    \section Example
        Here is a simple example demonstrating how to use this library.
        \include example.c
        In order to build this code, copy this code into a file and
        (assuming you have gcc installed), run:\n\n
        ```gcc <file> -o example.x -I<path to library include directory>
           -L<path to library lib directory> -lmolecular_lines```\n\n
        To run this example on your GPU, switch the use_gpu flag to any
        value other than 0, and recompile the library using nvcc instead
        of gcc.
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
    int use_gpu; /**< Flag telling if the lines will calculated on a GPU.*/
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
    @brief Initialize library parameters.
    @return 0 if completed successfully, or else an error code.
*/
#ifdef __NVCC__
extern "C"
#endif
int initialize_grt(GrtContext_t **context, /**< Library context.*/
                   int const num_levels, /**< Number of atmospheric levels.*/
                   double const w0, /**< Lowest wavenumber [1/cm] on spectral grid.*/
                   double const wn, /**< Highest wavenumber [1/cm] on spectral grid.*/
                   double const wres, /**< Spectral grid resolution [1/cm].*/
                   uint64_t * const num_wpoints, /**< Number of points on spectral grid.*/
                   double const * const wcutoff, /**< Cutoff [1/cm] from spectral line center.
                                                      Default value is 25.*/
                   int const * const use_gpu, /**< Flag to determine architecture where:\n
                                                   0 implies run on host CPU\n
                                                   != 0 implies run on GPU\n
                                                   Defaults to running on the host CPU.*/
                   int const * const num_threads, /**< If running on the host CPU,
                                                       determines the maximum number
                                                       of OpenMP threads that will
                                                       be used.  Defaults to
                                                       omp_get_max_threads (or one
                                                       if not build with OpenMP).*/
                   int const * const use_h2o_ctm, /**< Flag to determine if the water
                                                       vapor continumm will be included where:\n
                                                       0 implies no continuum\n
                                                       != 0 implies use the continuum\n
                                                       Defaults to running with the continuum.*/
                   int const * const use_o3_ctm /**< Flag to determine if the ozone
                                                     continumm will be included where:\n
                                                     0 implies no continuum\n
                                                     != 0 impiles use the continuum\n
                                                     Defaults to running with the continuum.*/
                  );


/**
    @ingroup capi
    @brief Release memory allocated by the library.
    @return 0 if completed successfully, or else an error code.
*/
#ifdef __NVCC__
extern "C"
#endif
int finalize_grt(GrtContext_t **context /**< Library context.*/
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
int add_molecule(GrtContext_t *context, /**< Library context.*/
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
int set_molecule_ppmv(GrtContext_t *context, /**< Library context.*/
                      int const molecule_id, /**< molecule_id !>Molecule id returned by
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
int calculate_optical_depth(GrtContext_t *context, /**< Library context.*/
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
