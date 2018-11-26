#ifndef LINE_SHAPE_H_
#define LINE_SHAPE_H_

#include "floating_point_type.h"


typedef struct LineShapeInputs
{
    fp_t w; /*Starting wavenumber where line shape is calculated [1/cm].*/
    uint64_t num_wpoints; /*Number of spectral points were the line shape
                            is calculated.*/
    fp_t wres; /*Spectral resolution [1/cm].*/
    fp_t line_center; /*Line center wavenumber [1/cm].*/
    fp_t lorentz_hwhm; /*Lorentz half-width at half-maximum [1/cm].*/
    fp_t doppler_hwhm; /*Doppler half-width at half-maximum [1/cm].*/
    fp_t eta; /*Mixing parameter for Ida voigt algorithm.*/
} LineShapeInputs_t;


#endif
