#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include "floating_point_type.h"


int malloc_ptr(void ** const p,
               size_t const num_bytes);


#define open_file(f,n,a) \
    {f=fopen(n,a); if (f == NULL) {fatal(IO_ERR,"failed to open file %s.",n);}}


int to_int(char const * const s,
           int * const i);


int to_double(char const * const s,
              double * const d);


int to_fp_t(double const d,
            fp_t * const f);


int linear_interp(double *in,
                  int in_size,
                  fp_t *out,
                  int out_size);


int reimann_sum(fp_t const * const data,
                int const data_size,
                fp_t const dx,
                fp_t * const out);


int get_sorted_bounds(fp_t const val,
                      fp_t const * const array,
                      int const array_size,
                      int * const left,
                      int * const right);


int linear_interpolation(fp_t const * const x,
                         fp_t const * const y,
                         int const xy_size,
                         fp_t const val,
                         fp_t * const out);


#endif
