#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <stdlib.h>
#include "floating_point_type.h"


typedef int(*sample1d_t)(fp_t const * const, fp_t const * const, fp_t const * const, fp_t * const, size_t);


fp_t angstrom_exponent(fp_t tau1,
                       fp_t tau2,
                       fp_t lambda1,
                       fp_t lambda2
                      );


int angstrom_exponent_sample(fp_t const * const x,
                             fp_t const * const y,
                             fp_t const * const newx,
                             fp_t * const newy,
                             size_t n
                            );


int integrate(fp_t const * const x,
              fp_t const * const y,
              size_t n,
              fp_t * const s
             );


int interpolate2(fp_t const * const x,
                 fp_t const * const y,
                 size_t n,
                 fp_t const * const newx,
                 fp_t * const newy,
                 size_t newn,
                 sample1d_t interp,
                 sample1d_t extrap
                );


int linear_sample(fp_t const * const x,
                  fp_t const * const y,
                  fp_t const * const newx,
                  fp_t * const newy,
                  size_t n
                 );


int monotonically_increasing(fp_t const * const x,
                             size_t n
                            );


#endif
