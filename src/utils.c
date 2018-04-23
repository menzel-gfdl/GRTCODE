#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include "constants.h"
#include "debug.h"
#include "floating_point_type.h"
#include "utils.h"


/*Helper function to perform a malloc, with error checks.*/
int malloc_ptr(void ** const p,
               size_t const num_bytes)
{
    not_null(p);
    *p = malloc(num_bytes);
    if (*p == NULL)
    {
        fatal("malloc of %zu bytes failed.",
              num_bytes);
    }
    return SUCCESS;
}


/*Helper function for converting a string to an integer.*/
int to_int(char *s,
           int *i)
{
    not_null(s);
    not_null(i);
    char *end;
    errno = 0;
    long n = strtol(s,&end,10);
    if (n == 0 && errno != 0)
    {
        if (errno == ERANGE)
        {
            fatal("the input string %s is out of range.",
                  s);
        }
        else if (end == s || errno == EINVAL)
        {
            fatal("invalid input string %s format, expecting the string to"
                      " contain an integer.",
                  s);
        }
        else
        {
            fatal("unknown errno code (%d) returned from strtol.",
                  errno);
        }
    }
    if (n >= INT_MIN && n <= INT_MAX)
    {
        *i = n;
    }
    else
    {
        fatal("input string %s cannot be represented as an int.",
              s);
    }
    return SUCCESS;
}


/*Helper function for converting a string to a double.*/
int to_double(char *s,
              double *d)
{
    not_null(s);
    not_null(d);
    char *end;
    errno = 0;
    *d = strtod(s,&end);
    if (*d == 0.0 && errno != 0)
    {
        if (errno == ERANGE)
        {
            fatal("the input string %s is out of range.",
                  s);
        }
        else if (end == s)
        {
            fatal("invalid input string %s format, expecting the string to"
                      " contain a floating point number.",
                  s);
        }
        else
        {
            fatal("unknown errno code (%d) returned from strtod.",
                  errno);
        }
    }
    return SUCCESS;
}


/*Helper function for converting a double to fp_t.*/
int to_fp_t(double const d,
            fp_t *f)
{
    not_null(f);
    if (sizeof(fp_t) == sizeof(double))
    {
        *f = d;
    }
    else if (sizeof(fp_t) == sizeof(float))
    {
        if (d >= -1.*FLT_MAX && d <= FLT_MAX)
        {
            *f = (fp_t)d;
        }
        else
        {
            fatal("input double value %le cannot be represented as a float.",
                  d);
        }
    }
    else
    {
        fatal("fp_t (size=%lu) must be represent either float or double.",
              sizeof(fp_t));
    }
    return SUCCESS;
}


/*Interpolate values to get new values.*/
int linear_interp(double *in,
                  int in_size,
                  fp_t *out,
                  int out_size)
{
    not_null(in);
    not_null(out);

    /*Make sure sizes are compatible.*/
    if (out_size != in_size + 1)
    {
        fatal("currently the output array size (%d) must be equal to the"
                  " input array size (%d) + 1.",
              out_size,
              in_size);
    }

    /*Handle edge cases.*/
    out[0] = in[0];
    out[out_size-1] = in[in_size-1];

    /*Do linear average.*/
    double const HALF = 0.5;
    int i;
    for (i=1;i<out_size-1;++i)
    {
        out[i] = (fp_t)(HALF*(in[i] + in[i-1]));
    }
    return SUCCESS;
}


/*Check that input bounds are vaild.*/
int input_bounds_check(int lower,
                       int min,
                       int *upper,
                       int max)
{
    not_null(upper);
    if (*upper < min)
    {
        *upper = max;
    }
    else if (*upper > max)
    {
        fatal("upper bound %d is greater than the maximum allowed (%d).",
              *upper,
              max);
    }
    if (lower < min)
    {
        fatal("lower bound %d is less than the minimum allowed (%d).",
              lower,
              min);
    }
    else if (lower > *upper)
    {
        fatal("lower bound %d is greater than upper bound %d.",
              lower,
              *upper);
    }
    return SUCCESS;
}


/*Make sure that a launch mode is valid.*/
int check_launch_mode(int const launch_type)
{
    if (launch_type != HOST_LAUNCH && launch_type != DEVICE_LAUNCH)
    {
        fatal("unsupported launch type %d.  Must be either HOST_LAUNCH (%d)"
                  " of DEVICE_LAUNCH (%d).",
              launch_type,
              HOST_LAUNCH,
              DEVICE_LAUNCH);
    }
    return SUCCESS;
}


/*Perform an integral using a reimann sum.*/
int reimann_sum(fp_t const * const data,
                int const data_size,
                fp_t const dx,
                fp_t * const out)
{
    not_null(data);
    not_null(out);
    *out = 0.;
    int i;
    for (i=0;i<(data_size-1);++i)
    {
        *out += dx*0.5*(data[i] + data[i+1]);
    }
    return SUCCESS;
}


int get_sorted_bounds(fp_t const val,
                      fp_t const * const array,
                      int const array_size,
                      int * const left,
                      int * const right)
{
    not_null(array);
    not_null(left);
    not_null(right);
    if (array_size <= 0)
    {
        fatal("input array size (%d) must be >= 1.",
              array_size);
    }

    if (val < array[0])
    {
        *left = -1;
        *right = 0;
    }
    else if (val > array[array_size-1])
    {
        *left = array_size - 1;
        *right = -1;
    }
    else
    {
        /*Since the input array is sorted, use a binary search.*/
        *left = 0;
        *right = array_size - 1;
        if (val == array[*left])
        {
            *right = *left;
        }
        else if (val == array[*right])
        {
            *left = *right;
        }
        else
        {
            while (1)
            {
                int mid = (*right + *left)/2;
                if (val == array[mid])
                {
                    *left = mid;
                    *right = mid;
                    break;
                }
                else if (val < array[mid])
                {
                    *right = mid;
                }
                else
                {
                    *left = mid;
                }
                if (*right - *left == 1)
                {
                    break;
                }
                else if (*right - *left == 0)
                {
                    fatal("This branch should never occur (val=%e).",
                          val);
                }
            }
        }
    }
    return SUCCESS;
}
