#ifndef FLOATING_POINT_TYPE_H_
#define FLOATING_POINT_TYPE_H_

#ifdef SINGLE_PRECISION
#define TYPE float
/*ln(max(float)) = 88.72284.  Let's use 80 so we have some runway.*/
#define MAX_EXP_ARG 80.f
#define EXP expf
#define POW powf
#define SQRT sqrtf
#define ABS fabsf
#else
#define TYPE double
/*ln(max(double)) = 709.782712893384.  Let's use 700 so we have some runway.*/
#define MAX_EXP_ARG 700.
#define EXP exp
#define POW pow
#define SQRT sqrt
#define ABS fabs
#endif


typedef TYPE fp_t;


#endif
