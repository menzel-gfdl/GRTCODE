#ifndef RAYLEIGH_H_
#define RAYLEIGH_H_

#include "floating_point_type.h"
#include "optics.h"


/** @brief Calculate the optical properties due to Rayleigh scattering.
    @return RS_SUCCESS or an error code.*/
int rayleigh_scattering(Optics_t * const optics, /**< Optics object.*/
                        fp_t * const pressure /**< Pressure [atm] (levels).*/
                       );


#endif
