#include <errno.h>
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
