#ifndef EVAL_PROFILE_H_
#define EVAL_PROFILE_H_

#include <stdint.h>
#include "floating_point_type.h"


#ifdef __NVCC__
__global__
void eval_profile(int const mol_id, /**< Molecule id.*/
                  unsigned int const num_lines, /**< Number of molecular
                                                     lines.*/
                  uint64_t const num_wpoints_fine, /**< Size of the fine
                                                        mesh.*/
                  uint64_t const num_wpoints_coarse, /**< Size of the coarse
                                                          mesh.*/
                  double const w0, /**< Lowest allowed wavenumber [1/cm].*/
                  double const wres_fine, /**< Resolution of the fine
                                               mesh [1/cm].*/
                  double const wres_coarse, /**< Resoluion of the coarse
                                                 mesh [1/cm].*/
                  int const num_layers, /**< Number of atmospheric layers.*/
                  double const wcutoff, /**< Cutoff from line center [1/cm].*/
                  fp_t const * const T, /**< Layer temperatures [K].*/
                  fp_t const * const gamma, /**< Pressure broadened line
                                                 halfwidths [1/cm].*/
                  fp_t const * const Pshift, /**< Pressure shifted line
                                                  centers [1/cm].*/
                  fp_t const * const s, /**< Line strengths [cm^2].*/
                  fp_t const * const N, /**< Integrated Layer number
                                             densities [cm^-2].*/
                  fp_t * const tau_fine, /**< Optical depths on the fine
                                              mesh.*/
                  fp_t * const tau_coarse, /**< Optical depths on the
                                                coarse mesh.*/
                  fp_t const fine_factor /**< Cutoff from the line center
                                              for the fine mesh in units
                                              of Pressure broadened line
                                              halwidths.*/
                 );
#endif


void eval_profile_h(int const mol_id,
                    unsigned int const num_lines,
                    uint64_t const num_wpoints_fine,
                    uint64_t const num_wpoints_coarse,
                    double const w0,
                    double const wres_fine,
                    double const wres_coarse,
                    int const num_layers,
                    double const wcutoff,
                    fp_t const * const T,
                    fp_t const * const gamma,
                    fp_t const * const Pshift,
                    fp_t const * const s,
                    fp_t const * const N,
                    fp_t * const tau_fine,
                    fp_t * const tau_coarse,
                    fp_t const fine_factor);


#endif
