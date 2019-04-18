#ifndef LAUNCH_H_
#define LAUNCH_H_

#include "floating_point_type.h"
#include "molecular_lines.h"


/** @brief Driver for optical depth calculation.
    @return RS_SUCCESS or an error code.*/
int launch(MolecularLines_t * const ml, /**< Molecular lines object.*/
           fp_t *p, /**< Pressure [atm] (level).*/
           fp_t *t, /**< Temperature [K] (level).*/
           fp_t * const tau /**< Optical depth (level, wavenumber).*/
          );


#endif
