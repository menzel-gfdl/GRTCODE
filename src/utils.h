#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include "floating_point_type.h"


int malloc_ptr(void ** const p,
               size_t const num_bytes);


#define open_file(f,n,a) \
    {f=fopen(n,a); if (f == NULL) {fatal("failed to open file %s.",n);}}


int to_int(char *s,
           int *i);


int to_double(char *s,
              double *d);


int linear_interp(double *in,
                  int in_size,
                  fp_t *out,
                  int out_size);


int input_bounds_check(int lower,
                       int min,
                       int *upper,
                       int max);


int check_launch_mode(int const launch_type);


#endif
