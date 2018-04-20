#ifndef CONTINUUM_HELPERS_H_
#define CONTINUUM_HELPERS_H_

#include "floating_point_type.h"


int get_coefs(char const * const filepath, /*Path to csv file.*/
              fp_t ** coefs, /*Array of pointers (one pointer for each column
                               in the file, except for the first column).*/
              int const num_coefs, /*The number of columns in the file - 1.*/
              unsigned int const num_grid_points, /*Model grid associated
                                                    with the first column in
                                                    the input file.*/
              int const first_grid_point, /*Value of the model grid at the
                                            first grid point.*/
              double const grid_spacing); /*Spacing between model grid points
                                            (assumed to be uniform).*/


#endif
