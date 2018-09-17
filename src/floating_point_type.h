#ifndef FLOATING_POINT_TYPE_H_
#define FLOATING_POINT_TYPE_H_

#ifdef SINGLE_PRECISION
#define TYPE float
#else
#define TYPE double
#endif

typedef TYPE fp_t;

#define cat(a,b) a##b
#define funcname(func,type) cat(func,type)


#endif
