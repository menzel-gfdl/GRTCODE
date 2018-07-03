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
    \details
    \section Overview
        Given an atmospheric column made up of at least one layer,
        this code calculates the total optical depth of each layer
        at each point on an input spectral grid.
    \section Example
        \include example.c
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


#endif
