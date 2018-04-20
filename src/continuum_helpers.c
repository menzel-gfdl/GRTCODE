#include <string.h>
#include "continuum_helpers.h"
#include "debug.h"
#include "floating_point_type.h"
#include "parse_csv.h"
#include "utils.h"


int get_coefs(char const * const filepath, /*Path to csv file.*/
              fp_t ** coefs, /*Array of pointers (one pointer for each column
                               in the file, except for the first column).*/
              int const num_coefs, /*The number of columns in the file - 1.*/
              unsigned int const num_grid_points, /*Model grid associated
                                                    with the first column in
                                                    the input file.*/
              int const first_grid_point, /*Value of the model grid at the
                                            first grid point.*/
              double const grid_spacing) /*Spacing between model grid points
                                           (assumed to be uniform).*/
{
    not_null(coefs);

    /*Parse the input file.*/
    int num_lines;
    int num_cols;
    char ** buf;
    check(parse_csv(filepath,
                    &num_lines,
                    &num_cols,
                    &buf));
    if (num_cols != num_coefs+1)
    {
        fatal("input csv file %s must have %d columns.",
              filepath,
              num_coefs+1);
    }
    if (num_lines <= 1)
    {
        fatal("input csv file %s does not contain any data.",
              filepath);
    }

    /*Convert the read in data from strings to floating point.*/
    int num_non_header_lines = num_lines - 1;
    fp_t fbuf[num_cols*num_non_header_lines];
    int i;
    for (i=0;i<num_cols;++i)
    {
        int j;
        for (j=0;j<num_non_header_lines;++j)
        {
            int offset1 = i*num_lines + j + 1;
            int offset2 = i*num_non_header_lines + j;
            double d;
            check(to_double(buf[offset1],
                            &d));
            check(to_fp_t(d,
                          &(fbuf[offset2])));
        }
    }
    for (i=0;i<(num_cols*num_lines);++i)
    {
        free(buf[i]);
    }
    free(buf);

    fp_t *x = &(fbuf[0]);
    for (i=0;i<num_coefs;++i)
    {
        /*Allocate necessary memory.*/
        check(malloc_ptr((void **)(&(coefs[i])),
                         sizeof(*(coefs[i]))*num_grid_points));
        memset(coefs[i],
               0,
               sizeof(*(coefs[i]))*num_grid_points);

        /*Calculate the coefficients at each grid point.*/
        fp_t *c = coefs[i];
        fp_t *y = &(fbuf[(i+1)*num_non_header_lines]);
        unsigned int j;
        for (j=0;j<num_grid_points;++j)
        {
            fp_t w = first_grid_point + j*grid_spacing;
            int left;
            int right;
            check(get_sorted_bounds(w,
                                    x,
                                    num_non_header_lines,
                                    &left,
                                    &right));
            if (left == right)
            {
                c[j] = y[left];
            }
            else if (left >= 0 && right >= 0)
            {
                /*Do a linear interpolation.*/
                fp_t m = (y[right] - y[left])/
                         (x[right] - x[left]);
                fp_t b = y[right] - m*x[right];
                c[j] = w*m + b;
            }
        }
    }
    return SUCCESS;
}
