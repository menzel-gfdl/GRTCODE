#ifndef FLOATING_POINT_TYPE_H_
#define FLOATING_POINT_TYPE_H_

#ifdef SINGLE_PRECISION
#define TYPE float
#define EXP expf
#define POW powf
#define SQRT sqrtf
#define ABS fabsf
#else
#define TYPE double
#define EXP exp
#define POW pow
#define SQRT sqrt
#define ABS fabs
#endif

typedef TYPE fp_t;

#define cat(a,b) a##b
#define funcname(func,type) cat(func,type)


#endif
