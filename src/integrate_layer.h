#ifndef INTEGRATE_LAYER_H_
#define INTEGRATE_LAYER_H_

#include "floating_point_type.h"


/*Integrate the total number density of the air across each atmospheric layer,
  assuming that:
    - each layer is hydrostatic.
    - accerlation due to gravity is constant across each layer.
*/
#ifdef __NVCC__
__global__
void integrated_N(int const num_layers, /**<Number of atmospheric
                                            layers.*/
                  fp_t const * const P, /**<Pressure [atm] at
                                            atmospheric layer edges.*/
                  fp_t * const N); /**<Integrated layer number density
                                       [cm^-2].*/
#endif


int integrated_N_h(int const num_layers,
                   fp_t const * const P,
                   fp_t * const N);


/*Calculate the Curtis-Godson integrals for pressure and temperature
  across each atmospheric layer, assuming that:
    - each layer is hydrostatic.
    - accerlation due to gravity is constant across each layer.
    - temperature is a linear function of pressure.
*/
#ifdef __NVCC__
__global__
void Curtis_Godson_PT(int const num_layers, /**<Number of atmospheric
                                                layers.*/
                      fp_t const * const P, /**<Pressure [atm] at
                                                atmospheric layer edges.*/
                      fp_t const * const T, /**<Temperature [K] at
                                                atmospheric layer edges.*/
                      fp_t * const Pavg, /**<Average layer pressure [atm].*/
                      fp_t * const Tavg); /**<Average layer temperature [K].*/
#endif


int Curtis_Godson_PT_h(int const num_layers,
                       fp_t const * const P,
                       fp_t const * const T,
                       fp_t * const Pavg,
                       fp_t * const Tavg);


/*Calculate the Curtis-Godson integrals for partial pressure and
  number density of a specific molecular species across each atmospheric
  layer, assuming that:
    - each layer is hydrostatic.
    - accerlation due to gravity is constant across each layer.
    - the molecular abundance is a linear function of pressure.
*/
#ifdef __NVCC__
__global__
void Curtis_Godson_PsNs(int const num_layers, /**<Number of atmospheric
                                                  layers.*/
                       fp_t const * const P, /**<Pressure [atm] at
                                                 atmospheric layer edges.*/
                       fp_t const * const x, /**<Molecular abundance at
                                                 atmospheric layer edges.*/
                       fp_t const * const N, /**<Total integrated layer
                                                 number density [cm^-2].*/
                       fp_t * const Psavg, /**<Average layer partial
                                               pressure [atm].*/
                       fp_t * const Ns); /**<Integrated number density [cm^-2]
                                             for the input molecule.*/
#endif


int Curtis_Godson_PsNs_h(int const num_layers,
                         fp_t const * const P,
                         fp_t const * const x,
                         fp_t const * const N,
                         fp_t * const Psavg,
                         fp_t * const Ns);


#endif
