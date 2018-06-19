#ifndef FLOATING_POINT_TYPE_H_
#define FLOATING_POINT_TYPE_H_

#ifdef DOUBLE_PRECISION
#define TYPE double
#else
#define TYPE float
#endif

typedef TYPE fp_t;

#define cat(a,b) a##b
#define funcname(func,type) cat(func,type)


#endif
