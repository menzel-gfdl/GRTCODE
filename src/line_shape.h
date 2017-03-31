#ifndef SET_LINE_SHAPE_H_
#define SET_LINE_SHAPE_H_

#include "myreal.h"

struct LineShapeInputs
{
    REAL_t freq;       /*Frequency where line shape is calculated (1/cm).*/
    REAL_t lineCenter; /*Line center frequency (1/cm).*/
    REAL_t lorHWHM;    /*Lorentz half-width at half-maximum (1/cm).*/
    REAL_t gauHWHM;    /*Gaussian half-width at half-maximum (1/cm).*/
    REAL_t eta;        /*Mixing parameter for Ida voigt algorithm.*/
};
typedef struct LineShapeInputs LineShapeInputs_t;

#endif
