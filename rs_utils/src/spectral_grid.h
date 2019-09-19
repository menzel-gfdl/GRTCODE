#ifndef SPECTRAL_GRID_H_
#define SPECTRAL_GRID_H_

#include <stdint.h>
#include "extern.h"


/** @brief Container for spectral grid properties.*/
typedef struct SpectralGrid
{
    double w0; /**< Lower bound [1/cm].*/
    double wn; /**< Upper bound [1/cm].*/
    double dw; /**< Grid spacing [1/cm.*/
    uint64_t n; /**< Number of grid points.*/
} SpectralGrid_t;


/** @brief Initialize a spectral grid.
    @return RS_SUCCESS or an error code.*/
EXTERN int create_spectral_grid(SpectralGrid_t * const grid, /**< Spectral grid object.*/
                                double const w0, /**< Lower bound [1/cm].*/
                                double const wn, /**< Upper bound [1/cm].*/
                                double const dw /**< Grid spacing [1/cm.*/
                               );


/** @brief Determine if two spectral grids are the same.
    @return RS_SUCCESS or an error code.*/
EXTERN int compare_spectral_grids(SpectralGrid_t const * const one, /**< Spectral grid object.*/
                                  SpectralGrid_t const * const two, /**< Spectral grid object.*/
                                  int * const result /**< 1 if they match, 0 if not.*/
                                 );


#endif
