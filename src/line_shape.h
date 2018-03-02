#ifndef LINE_SHAPE_H_
#define LINE_SHAPE_H_

#include "floating_point_type.h"


struct LineShapeInputs
{
    fp_t freq;       /*Frequency where line shape is calculated (1/cm).*/
    fp_t lineCenter; /*Line center frequency (1/cm).*/
    fp_t lorHWHM;    /*Lorentz half-width at half-maximum (1/cm).*/
    fp_t gauHWHM;    /*Gaussian half-width at half-maximum (1/cm).*/
    fp_t eta;        /*Mixing parameter for Ida voigt algorithm.*/
};
typedef struct LineShapeInputs LineShapeInputs_t;


#endif
