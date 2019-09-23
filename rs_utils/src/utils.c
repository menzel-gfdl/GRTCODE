#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "extern.h"
#include "floating_point_type.h"
#include "utils.h"


/*Malloc memory, with error checks.*/
EXTERN int malloc_ptr(void ** const p, size_t const num_bytes)
{
    not_null(p);
    *p = malloc(num_bytes);
    not_null(*p);
    return RS_SUCCESS;
}


/*Free malloced memory, with error checks.*/
EXTERN int free_ptr(void ** const p)
{
    not_null(p);
    not_null(*p);
    free(*p);
    *p = NULL;
    return RS_SUCCESS;
}


/*Copy a string into a buffer, checking its length.*/
int copy_str(char * const dest, char const * const src, size_t const len)
{
    if (strlen(src) > len)
    {
        char const *mesg = "input string (%s) is larger than the input buffer"
                           " and would be truncated.";
        raise(RS_VALUE_ERR, mesg, src);
    }
    snprintf(dest, len, "%s", src);
    return RS_SUCCESS;
}


/*Open a file, with error checks.*/
int open_file(FILE **file, char const * const name, char const * const mode)
{
    not_null(file);
    not_null(name);
    not_null(mode);
    *file = fopen(name, mode);
    not_null(*file);
    return RS_SUCCESS;
}


/*Convert a string to an integer.*/
int to_int(char const * const s, int * const i)
{
    not_null(s);
    not_null(i);
    char *end;
    errno = 0;
    long n = strtol(s, &end, 10);
    if ((errno == ERANGE && (n == LONG_MAX || n == LONG_MIN)) ||
        (n == 0 && errno != 0))
    {
        char const *mesg = "the input string %s is out of range.";
        raise(RS_RANGE_ERR, mesg, s);
    }
    if (end == s || errno == EINVAL)
    {
        char const *mesg = "invalid input string %s, expecting the string to"
                           " contain an integer.";
        raise(RS_VALUE_ERR, mesg, s);
    }
    if (n >= INT_MIN && n <= INT_MAX)
    {
        *i = n;
    }
    else
    {
        char const *mesg = "input string %s cannot be represented as an int.";
        raise(RS_VALUE_ERR, mesg, s);
    }
    return RS_SUCCESS;
}


/*Convert a string to a double.*/
int to_double(char const * const s, double * const d)
{
    not_null(s);
    not_null(d);
    char *end;
    errno = 0;
    *d = strtod(s, &end);
    if ((errno == ERANGE && (*d == HUGE_VAL || *d == -1.*HUGE_VAL)) ||
        (*d == 0. && errno != 0))
    {
        char const *mesg = "the input string %s is out of range.";
        raise(RS_RANGE_ERR, mesg, s);
    }
    if (end == s)
    {
        char const *mesg = "invalid input string %s, expecting the string to"
                           " contain a floating point number.";
        raise(RS_VALUE_ERR, mesg, s);
    }
    return RS_SUCCESS;
}


/*Convert a double to fp_t.*/
int to_fp_t(double const d, fp_t * const f)
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
            char const *mesg = "input double value %le cannot be represented as a float.";
            raise(RS_RANGE_ERR, mesg, d);
        }
    }
    else
    {
        char const *mesg = "fp_t (size=%lu) must be represent either float or double.";
        raise(RS_VALUE_ERR, mesg, sizeof(fp_t));
    }
    return RS_SUCCESS;
}


/*Find the array indices that bracket the input value.*/
int get_sorted_bounds(fp_t const val, fp_t const * const array, int const array_size,
                      int * const left, int * const right)
{
    not_null(array);
    not_null(left);
    not_null(right);
    min_check(array_size, 1)
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
                    sentinel();
                }
            }
        }
    }
    return RS_SUCCESS;
}


/*Perform a linear interpolation.*/
int linear_interpolation(fp_t const * const x, fp_t const * const y, int const xy_size,
                         fp_t const val, fp_t * const out)
{
    not_null(x);
    not_null(y);
    not_null(out);
    int left;
    int right;
    catch(get_sorted_bounds(val, x, xy_size, &left, &right));
    if (left == right)
    {
        *out = y[left];
    }
    else if (left >= 0 && right >= 0)
    {
        /*Do a linear interpolation.*/
        fp_t m = (y[right] - y[left])/(x[right] - x[left]);
        fp_t b = y[right] - m*x[right];
        *out = val*m + b;
    }
    return RS_SUCCESS;
}


/*Perform an integral using a reimann sum.*/
int reimann_sum(fp_t const * const data, int const data_size, fp_t const dx,
                fp_t * const out)
{
    not_null(data);
    not_null(out);
    *out = 0.;
    int i;
    for (i=0; i<(data_size-1); ++i)
    {
        *out += dx*0.5*(data[i] + data[i+1]);
    }
    return RS_SUCCESS;
}


/*Turn on bit in bit field.*/
int activate(uint64_t * const bit_field, int const index)
{
    not_null(bit_field);
    in_range(index, 0, 63);
    uint64_t const one = 1;
    *bit_field = (*bit_field) | (one << index);
    return RS_SUCCESS;
}


/*Check if bit is turned on in bit field.*/
int is_active(uint64_t const bit_field, int const index)
{
    in_range(index, 0, 63);
    uint64_t const one = 1;
    return bit_field & (one << index);
}
