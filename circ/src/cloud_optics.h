#ifndef CLOUD_OPTICS_H_
#define CLOUD_OPTICS_H_

#include "extern.h"
#include "floating_point_type.h"
#include "optics.h"
#include "spectral_grid.h"


typedef int(*LiquidCloud_t)(fp_t const, fp_t const, SpectralGrid_t const,
                            fp_t * const, fp_t * const, fp_t * const);


/** @brief Calculate cloud optics using the parameterization described in
           https://doi.org/10.1175/1520-0442(1993)006<0728:AAPOTR>2.0.CO;2
    @return RS_SUCCESS or an error code.*/
EXTERN int hu_stamnes_1993(fp_t const liquid_water_path, /**< Liquid water path [g/m^3].*/
                           fp_t const droplet_equivalent_radius, /**< Droplet equivalent radius [microns].*/
                           SpectralGrid_t const grid, /**< Spectral grid object.*/
                           fp_t * const optical_depth, /**< Optical depth (wavenumber).*/
                           fp_t * const single_scatter_albedo, /**< Single-scatter albedo (wavenumber).*/
                           fp_t * const asymmetry_factor /**< Asymmetry factor (wavenumber).*/
                          );


/** @brief Calculate cloud optics using the parameterization described in
           https://doi.org/10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2
    @return RS_SUCCESS or an error code.*/
EXTERN int slingo_1989(fp_t const liquid_water_path, /**< Liquid water path [g/m^3].*/
                       fp_t const droplet_equivalent_radius, /**< Droplet equivalent radius [microns].*/
                       SpectralGrid_t const grid, /**< Spectral grid object.*/
                       fp_t * const optical_depth, /**< Optical depth (wavenumber).*/
                       fp_t * const single_scatter_albedo, /**< Single-scatter albedo (wavenumber).*/
                       fp_t * const asymmetry_factor /**< Asymmetry factor (wavenumber).*/
                      );


/** @brief Calculate liquid cloud optics.
    @return RS_SUCCESS or an error code.*/
EXTERN int cloud_optics(Optics_t * const optics, /**< Optics object.*/
                        fp_t * const liquid_water_path, /**< Liquid water path [g/m^3] (layer).*/
                        fp_t * const droplet_equivalent_radius, /**< Droplet equivalent radius [microns] (layer).*/
                        LiquidCloud_t parameterization /**< Parameterization function pointer.*/
                       );


#endif
