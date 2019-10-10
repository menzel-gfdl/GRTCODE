#include <math.h>
#include <stdlib.h>
#include "debug.h"
#include "floating_point_type.h"
#include "utilities.h"



fp_t angstrom_exponent(fp_t tau1, fp_t tau2, fp_t lambda1, fp_t lambda2)
{
    fp_t const c = -1.;
    return c*log(tau1/tau2)/log(lambda1/lambda2);
}


int angstrom_exponent_sample(fp_t const * const x, fp_t const * const y,
                             fp_t const * const newx, fp_t * const newy, size_t n)
{
    size_t i;
    for (i=0; i<2; ++i)
    {
        if (y[i] <= 0.)
        {
            char const * mesg = "Cannot calculate the angstrom exponent because"
                                " y[%zu] <= 0 (%e)";
            raise(RS_VALUE_ERR, mesg, i, y[i]);
        }
    }
    fp_t const alpha = -1.*angstrom_exponent(y[1], y[0], x[0], x[1]);
    for (i=0; i<n; ++i)
    {
        newy[i] = y[0]*pow((x[0]/newx[i]), alpha);
    }
    return RS_SUCCESS;
}


int integrate(fp_t const * const x, fp_t const * const y, size_t n, fp_t * const s)
{
    if (!monotonically_increasing(x, n))
    {
        char const * mesg = "x (%p) must be monotonically increasing.";
        raise(RS_VALUE_ERR, mesg, x);
    }
    fp_t const half = 0.5;
    *s = 0.;
    size_t i;
    for (i=0; i<(n-1); ++i)
    {
        *s += half*(y[i] + y[i+1])*(x[i+1] - x[i]);
    }
    return RS_SUCCESS;
}


int interpolate2(fp_t const * const x, fp_t const * const y, size_t n,
                 fp_t const * const newx, fp_t * const newy, size_t newn,
                 sample1d_t interp, sample1d_t extrap)
{
    if (n < 2)
    {
        char const * mesg = "at least two x points required (%zu given).";
        raise(RS_VALUE_ERR, mesg, n);
    }
    if (!monotonically_increasing(x, n))
    {
        char const * mesg = "x (%p) must be monotonically increasing.";
        raise(RS_VALUE_ERR, mesg, x);
    }
    if (!monotonically_increasing(newx, newn))
    {
        char const * mesg = "newx (%p) must be monotonically increasing.";
        raise(RS_VALUE_ERR, mesg, newx);
    }
    size_t i;
    for (i=0; i<newn; ++i)
    {
        if (newx[i] > x[0])
        {
            break;
        }
    }
    if (i > 0)
    {
        /*Handle all points where newx[:i-1] <= x[0].*/
        if (extrap != NULL)
        {
            catch(extrap(x, y, newx, newy, i));
        }
        else
        {
            size_t j;
            for (j=0; j<i; ++j)
            {
                newy[j] = y[0];
            }
        }
        if (i == newn)
        {
            /*All newx[:] <= x[0].*/
            return RS_SUCCESS;
        }
    }
    size_t j;
    for (j=0; j<n-1; ++j)
    {
        size_t k;
        for (k=i; k<newn; ++k)
        {
            if (newx[k] > x[j+1])
            {
                break;
            }
        }
        if (k > i)
        {
            /*Handle all points where x[j] < newx[i:k-1] <= x[j+1].*/
            catch(interp(&(x[j]), &(y[j]), &(newx[i]), &(newy[i]), k-i));
            i = k;
            if (i == newn)
            {
                /*All newx points have been handled.*/
                return RS_SUCCESS;
            }
        }
    }
    /*Handle all points where newx[i:] > x[n-1].*/
    if (extrap != NULL)
    {
        catch(extrap(&(x[n-2]), &(y[n-2]), &(newx[i]), &(newy[i]), newn-i));
    }
    else
    {
        for (j=i; j<newn; ++j)
        {
            newy[j] = y[n-1];
        }
    }
    return RS_SUCCESS;
}


int linear_sample(fp_t const * const x, fp_t const * const y,
                  fp_t const * const newx, fp_t * const newy, size_t n)
{
    fp_t const m = (y[1] - y[0])/(x[1] - x[0]);
    fp_t const b = y[0] - m*x[0];
    size_t i;
    for (i=0; i<n; ++i)
    {
        newy[i] = m*newx[i] + b;
    }
    return RS_SUCCESS;
}


int monotonically_increasing(fp_t const * const x, size_t n)
{
  size_t i;
  for (i=0; i<(n-1); ++i)
  {
      if (x[i+1] <= x[i])
      {
          return 0;
      }
  }
  return 1;
}
