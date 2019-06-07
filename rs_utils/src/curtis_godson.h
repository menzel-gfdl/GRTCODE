#ifndef CURTIS_GODSON_H_
#define CURTIS_GODSON_H_

#include "floating_point_type.h"


/** @brief Calculate integrated number densities.
    @return RS_SUCCESS or an error code.*/
int calc_number_densities(int const num_layers, /**< Number of atmospheric layers.*/
                          fp_t const * const p, /**< Pressure [atm] (levels).*/
                          fp_t * const n /**< Integrated number densities [cm^-2] (layers).*/
                         );


#ifdef __NVCC__
/** @brief Calculate integrated number densities.*/
__global__ void calc_number_densities_d(int const num_layers, /**< Number of atmospheric layers.*/
                                        fp_t const * const p, /**< Pressure [atm] (levels).*/
                                        fp_t * const n /**< Integrated number densities [cm^-2] (layers).*/
                                       );
#endif


/** @brief Calculate layer pressures and temperatures.
    @return RS_SUCCESS or an error code.*/
int calc_pressures_and_temperatures(int const num_layers, /**< Number of atmospheric layers.*/
                                    fp_t const * const p, /**< Pressure [atm] (levels).*/
                                    fp_t const * const t, /**< Temperature [K] (levels).*/
                                    fp_t * const pavg, /**< Pressure [atm] (layers).*/
                                    fp_t * const tavg /**< Pressure [atm] (layers).*/
                                   );


#ifdef __NVCC__
/** @brief Calculate layer pressures and temperatures.*/
__global__ void calc_pressures_and_temperatures_d(int const num_layers, /**< Number of atmospheric layers.*/
                                                  fp_t const * const p, /**< Pressure [atm] (levels).*/
                                                  fp_t const * const t, /**< Temperature [K] (levels).*/
                                                  fp_t * const pavg, /**< Pressure [atm] (layers).*/
                                                  fp_t * const tavg /**< Pressure [atm] (layers).*/
                                                 );
#endif


/** @brief Calculate partial pressures and number densities.
    @return RS_SUCCESS or an error code.*/
int calc_partial_pressures_and_number_densities(int const num_layers, /**< Number of atmospheric layers.*/
                                                fp_t const * const p, /**< Pressure [atm] (levels).*/
                                                fp_t const * const x, /**< Abundance (levels).*/
                                                fp_t const * const n, /**< Integrated number densities [cm^-2] (layers).*/
                                                fp_t * const ps, /**< Partial pressure [atm] (layers).*/
                                                fp_t * const ns /**< Integrated molecular number densities [cm^-2] (layers).*/
                                               );


#ifdef __NVCC__
/** @brief Calculate partial pressures and number densities.*/
__global__ void calc_partial_pressures_and_number_densities_d(int const num_layers, /**< Number of atmospheric layers.*/
                                                              fp_t const * const p, /**< Pressure [atm] (levels).*/
                                                              fp_t const * const x, /**< Abundance (levels).*/
                                                              fp_t const * const n, /**< Integrated number densities [cm^-2] (layers).*/
                                                              fp_t * const ps, /**< Partial pressure [atm] (layers).*/
                                                              fp_t * const ns /**< Integrated molecular number densities [cm^-2] (layers).*/
                                                             );
#endif


#endif
