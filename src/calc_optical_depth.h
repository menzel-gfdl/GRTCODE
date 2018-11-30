#ifndef CALC_OPTICAL_DEPTH_H_
#define CALC_OPTICAL_DEPTH_H_

#include <stdint.h>
#include "floating_point_type.h"
#include "spectral_bin.h"


/** @brief Calculate integrated number densities.
    @return SUCCESS or an error code.*/
int calc_number_densities(int const num_layers, /*Number of atmospheric layers.*/
                          fp_t const * const p, /*Pressure [atm] (levels).*/
                          fp_t * const n /*Integrated number densities
                                           [cm^-2] (layers).*/
                         );


/** @brief Calculate layer pressures and temperatures.
    @return SUCCESS or an error code.*/
int calc_pressures_and_temperatures(int const num_layers, /*Number of atmospheric
                                                            layers.*/
                                    fp_t const * const p, /*Pressure [atm] (levels).*/
                                    fp_t const * const t, /*Temperature [K] (levels).*/
                                    fp_t * const pavg, /*Pressure [atm] (layers).*/
                                    fp_t * const tavg /*Pressure [atm] (layers).*/
                                   );


/** @brief Calculate partial pressures and number densities.
    @return SUCCESS or an error code.*/
int calc_partial_pressures_and_number_densities(int const num_layers, /*Number of
                                                                        atmospheric layers.*/
                                                fp_t const * const p, /*Pressure [atm] (levels).*/
                                                fp_t const * const x, /*Abundance (levels).*/
                                                fp_t const * const n, /*Integrated number densities
                                                                        [cm^-2] (layers).*/
                                                fp_t * const ps, /*Partial pressure [atm]
                                                                   (layers).*/
                                                fp_t * const ns /*Integrated molecular number
                                                                  densities [cm^-2] (layers).*/
                                               );


/** @brief Calculate pressure-shifted line center positions.
    @return SUCCESS or an error code.*/
int calc_line_centers(uint64_t const num_lines, /*Number of molecular lines.*/
                      int const num_layers, /*Number of atmospheric layers.*/
                      fp_t const * const v0, /*Unshifted line center
                                               positions [1/cm] (lines).*/
                      fp_t const * const delta, /*Air-broadened pressure
                                                  shift [1/(cm*atm)] (lines).*/
                      fp_t const * const p, /*Pressure [atm] (layers).*/
                      fp_t * const vnn /*Pressure-shifted line center
                                         positions [1/cm] (layers,lines).*/
                     );


/** @brief Calculate temperature-corrected line intensities.
    @return SUCCESS or an error code.*/
int calc_line_strengths(uint64_t const num_lines, /*Number of molecular lines.*/
                        int const num_layers, /*Number of atmospheric layers.*/
                        int const mol_id, /*Molecule id.*/
                        int const num_iso, /*Number of molecular
                                             isotopologues.*/
                        int const * const iso, /*Isotopologue id (lines).*/
                        fp_t const * const s0, /*Uncorrected line strengths
                                                 [1/cm] (lines).*/
                        fp_t const * const vnn, /*Line center position [1/cm]
                                                  (lines).*/
                        fp_t const * const en, /*Lower state energies [1/cm]
                                                 (lines).*/
                        fp_t const * const t, /*Temperature [K] (layers).*/
                        fp_t * const snn /*Temperature-corrected line
                                           strengths [1/cm] (layers,lines).*/
                       );


/** @brief Calculate lorentz halfwidths.
    @return SUCCESS or an error code.*/
int calc_lorentz_hw(uint64_t const num_lines, /*Number of molecular lines.*/
                    int const num_layers, /*Number of atmospheric layers.*/
                    fp_t const * const n, /*Coefficient of temperature
                                            dependence of air-broadened
                                            halfwidths (lines).*/
                    fp_t const * const yair, /*Air-broadended halfwidths
                                               [1/(cm*atm)] at 296K and 1atm.*/
                    fp_t const * const yself, /*Self-broadended halfwidths
                                               [1/(cm*atm)] at 296K and 1atm.*/
                    fp_t const * const t, /*Temperature [K] (layers).*/
                    fp_t const * const p, /*Pressure [atm] (layers).*/
                    fp_t const * const ps, /*Partial pressure [atm] (layers.*/
                    fp_t * const gamma /*Temperature and pressure corrected
                                         lorentz halfwidths [1/cm]
                                         (layers,lines).*/
                   );


/** @brief Calculate doppler halfwidths.
    @return SUCCESS or an error code.*/
int calc_doppler_hw(uint64_t const num_lines, /*Number of molecular lines.*/
                    int const num_layers, /*Number of atmospheric layers.*/
                    fp_t const m, /*Molecular mass [g].*/
                    fp_t const * const vnn, /*Pressure-shifted line center
                                              position [1/cm] (layers,lines).*/
                    fp_t const * const t, /*Temperature [K] (layers).*/
                    fp_t * const alpha /*Doppler halfwidths [1/cm]
                                         (layers,lines).*/
                   );


/** @brief Calculate optical depths.
    @return SUCCESS or an error code.*/
int calc_optical_depth(uint64_t const num_lines, /*Number of molecular lines.*/
                       int const num_layers, /*Number of atmospheric layers.*/
                       fp_t * const vnn, /*Pressure-shifted line
                                           center positions [1/cm].
                                           (layers,lines).*/
                       fp_t * const snn, /*Line strength [1/cm]
                                           (layers,lines).*/
                       fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                             (layers,lines).*/
                       fp_t * const alpha, /*Doppler halfwidth [1/cm]
                                             (layers,lines).*/
                       fp_t const * const n, /*Integrated number density
                                               [cm^-2] (layers).*/
                       SpectralBins_t * const bins, /*Spectral bins.*/
                       fp_t * const tau /*Optical depth (layer,wavenumber).*/
                      );


/** @brief Calculate optical depths.
    @return SUCCESS or an error code.*/
int calc_optical_depth_2(uint64_t const num_lines, /*Number of molecular lines.*/
                         int const num_layers, /*Number of atmospheric layers.*/
                         fp_t * const vnn, /*Pressure-shifted line
                                             center positions [1/cm].
                                             (layers,lines).*/
                         fp_t * const snn, /*Line strength [1/cm]
                                             (layers,lines).*/
                         fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                               (layers,lines).*/
                         fp_t * const alpha, /*Doppler halfwidth [1/cm]
                                               (layers,lines).*/
                         fp_t const * const n, /*Integrated number density
                                                 [cm^-2] (layers).*/
                         SpectralBins_t * const bins, /*Spectral bins.*/
                         fp_t * const tau /*Optical depth (layer,wavenumber).*/
                        );


/** @brief Calculate optical depths.
    @return SUCCESS or an error code.*/
int calc_optical_depth_old(uint64_t const num_lines, /*Number of molecular lines.*/
                           int const num_layers, /*Number of atmospheric layers.*/
                           fp_t * const vnn, /*Pressure-shifted line
                                               center positions [1/cm].
                                               (layers,lines).*/
                           fp_t * const snn, /*Line strength [1/cm]
                                               (layers,lines).*/
                           fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                                 (layers,lines).*/
                           fp_t * const alpha, /*Doppler halfwidth [1/cm]
                                                 (layers,lines).*/
                           fp_t const * const n, /*Integrated number density
                                                   [cm^-2] (layers).*/
                           SpectralBins_t const * const bins, /*Spectral bins.*/
                           fp_t * const tau /*Optical depth (layer,wavenumber).*/
                          );


#endif
