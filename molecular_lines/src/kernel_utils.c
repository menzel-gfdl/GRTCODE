#include <stdint.h>
#include "debug.h"
#include "floating_point_type.h"
#include "kernel_utils.h"


/** @brief Get the indices of an array that bracket a value.  The input
           array must be sorted.
    @return SUCCESS or an error code.*/
HOST DEVICE int bracket(uint64_t const array_size, /**< Size of input array.*/
                        fp_t const * const array, /**< Input array.*/
                        fp_t const val, /**< Value to bracket.*/
                        uint64_t * const left, /**< Index of maximum array
                                                    value < val.*/
                        uint64_t * const right /**< Index of minimum array
                                                    value > val.*/
                       )
{
    not_null(array);
    not_null(left);
    not_null(right);
    min_check(array_size,1);
    uint64_t l = 0;
    uint64_t r = array_size - 1;
    if (val < array[l] || val > array[r])
    {
        *left = l;
        *right = r;
        return RANGE_ERR;
    }
    if (array[l] == val)
    {
        r = l;
    }
    else if (array[r] == val)
    {
        l = r;
    }
    else
    {
        while (r-l > 1)
        {
            uint64_t mid = l + (r-l)/2;
            if (array[mid] == val)
            {
                l = mid;
                r = mid;
                break;
            }
            else if (val > array[mid])
            {
                l = mid;
            }
            else
            {
                r = mid;
            }
            if (l > r || l == r)
            {
                raise(VALUE_ERR,
                      "Something went wrong (l=%zu,r=%zu).",
                      l,
                      r);
            }
        }
    }
    *left = l;
    *right = r;
    return SUCCESS;
}


/** @brief Perform a quadratic interpolation, but set all values less than
           zero equal to zero.
    @return SUCCESS or an error code.*/
HOST DEVICE int bin_quad_interp(fp_t const * const x, /**<*/
                                fp_t const * const y, /**<*/
                                uint64_t const left, /**<*/
                                uint64_t const right, /**<*/
                                fp_t const w0, /**<*/
                                fp_t const wres, /**<*/
                                fp_t * const tau /**<*/
                               )
{
    not_null(x);
    not_null(y);
    not_null(tau);
    uint64_t j;
    for (j=left;j<=right;++j)
    {
        fp_t w = w0 + j*wres;
        fp_t t = (w-x[1])*(w-x[2])*y[0]/((x[0]-x[1])*(x[0]-x[2])) +
                 (w-x[0])*(w-x[2])*y[1]/((x[1]-x[0])*(x[1]-x[2])) +
                 (w-x[0])*(w-x[1])*y[2]/((x[2]-x[0])*(x[2]-x[1]));
        if (t < 0.f)
        {
            t = 0.f;
        }
        tau[j] += t;
    }
    return SUCCESS;
}


/** @brief Copy optical depths from the coarse mesh to the fine mesh in a
           bin.
    @return SUCCESS or an error code.*/
HOST DEVICE int bin_no_interp(uint64_t const left,
                              uint64_t const right,
                              fp_t const * const taub,
                              fp_t * const tau
                             )
{
    not_null(taub);
    not_null(tau);
    uint64_t j;
    for (j=left;j<=right;++j)
    {
        tau[j] += taub[j-left];
    }
    return SUCCESS;
}
