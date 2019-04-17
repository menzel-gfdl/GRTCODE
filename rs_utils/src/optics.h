#ifndef OPTICS_H_
#define OPTICS_H_

#include "floating_point_type.h"
#include "spectral_grid.h"


/** @brief Container for optical properties.*/
typedef struct Optics
{
    int num_layers; /**< Number of atmospheric pressure layers.*/
    SpectralGrid_t grid; /**< Spectral grid object.*/
    fp_t *g; /**< Asymmetric factor (layers,wavenumber).*/
    fp_t *omega; /**< Single-scattering albedo (layers,wavenumber).*/
    fp_t *tau; /**< Optical depth (layers,wavenumber).*/
} Optics_t;


/** @brief Reserve memory for the optics.
    @return RS_SUCCESS or an error code.*/
int create_optics(Optics_t * const optics, /**< Optics object.*/
                  int const num_layers, /**< Number of atmospheric layers.*/
                  SpectralGrid_t const * const grid /**< Spectral grid object.*/
                 );


/** @brief Free memory for the optics.
    @return RS_SUCCESS or an error code.*/
int destroy_optics(Optics_t * const optics /**< Optics object.*/
                  );


/** @brief Determine if two optics objects are compatible.
    @return RS_SUCCESS or an error code.*/
int optics_compatible(Optics_t const * const one, /**< Optics object.*/
                      Optics_t const * const two, /**< Optics object.*/
                      int * const result /**< 1 if compatable, 0 if not.*/
                     );


/** @brief Add optical properties together.
    @return RS_SUCCESS or an error code.*/
int add_optics(Optics_t const * const * const optics, /**< Array of optics objects.*/
               int const num_optics, /**< Size of the input array of optics objects.*/
               Optics_t * const result /**< Resulting optics objects.*/
              );


#endif
