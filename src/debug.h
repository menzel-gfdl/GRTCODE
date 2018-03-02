#ifndef DEBUG_H_
#define DEBUG_H_

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "netcdf.h"


enum return_codes
{
    SUCCESS = 0,
    ERR
};


#define backtrace() {\
    fprintf(stderr, \
            "\t%s: %d\n", \
            __FILE__, \
            __LINE__);}


#define log_err(mesg,...) {\
    fprintf(stderr, \
            "[Error] %s: " mesg "\nBacktrace:\n", \
            __func__, \
            __VA_ARGS__); \
    backtrace();}


#define log_mesg(mesg,...) {\
    fprintf(stderr, \
            "[%s:%d]: " mesg "\n", \
            __FILE__, \
            __LINE__, \
            __VA_ARGS__);}


#define fatal(mesg,...) {\
    log_err(mesg, \
            __VA_ARGS__); \
    return ERR;}


#ifdef __NVCC__
#define using_gpu() {}
#else
#define using_gpu() {\
    fatal("%s","you cannot make CUDA calls unless you compile with nvcc." \
               "  To run on only a host CPU, include the -h option.");}
#endif


#define check(val) {\
    int e_ = val; \
    if (e_ != SUCCESS) \
    { \
        backtrace(); \
        return e_; \
    }}


#define not_null(p) {\
    if (p == NULL) \
    { \
        fatal("null pointer at address %p.",(void *)(&p)); \
    }}


#define is_null(p) {\
    if (p != NULL) \
    { \
        fatal("pointer at address %p is not null.",(void *)(&p)); \
    }}


#define netcdf_check(val) {\
    int e_ = val; \
    if (e_ != NC_NOERR) \
    { \
        fatal("netcdf returned error code %d. %s.", \
              e_, \
              nc_strerror(e_)); \
    }}


#ifdef __NVCC__
#define kernel_err(mesg,...) {}
#else
#define kernel_err(mesg,...) {log_err(mesg,__VA_ARGS__);exit(ERR);}
#endif


#endif
