#ifndef DEBUG_H_
#define DEBUG_H_

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>


enum return_codes
{
    SUCCESS,
    INVALID_ERR,
    DIVBYZERO_ERR,
    OVERFLOW_ERR,
    UNDERFLOW_ERR,
    SENTINEL_ERR,
    NULL_ERR,
    NON_NULL_ERR,
    RANGE_ERR,
    VALUE_ERR,
    COMPILER_ERR,
    IO_ERR
/*
    SUCCESS = GRT_SUCCESS,
    INVALID_ERR = GRT_INVALID,
    DIVBYZERO_ERR = GRT_DIVBYZERO,
    OVERFLOW_ERR = GRT_OVERFLOW,
    UNDERFLOW_ERR = GRT_UNDERFLOW,
    SENTINEL_ERR = GRT_SENTINEL,
    NULL_ERR = GRT_NULL,
    NON_NULL_ERR = GRT_NON_NULL,
    RANGE_ERR = GRT_OUT_OF_RANGE,
    VALUE_ERR = GRT_VALUE_ERR,
    COMPILER_ERR = GRT_COMPILER_ERR,
    IO_ERR = GRT_IO_ERR
*/
};


#ifdef __CUDA_ARCH__
#undef DEBUG
#undef VERBOSE
#undef FAST
#else
#ifdef DEBUG
#define VERBOSE
#undef FAST
#endif
#endif


/*Macros that provide debugging information.*/
#ifdef VERBOSE


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


#define log_warn(mesg,...) {\
    fprintf(stderr, \
            "[Warn] (%s:%d)" mesg "\n", \
            __FILE__, \
            __LINE__, \
            __VA_ARGS__);}


#else


#define backtrace() {}
#define log_err(mesg,...) {}
#define log_warn(mesg,...) {}


#endif


#define log_mesg(mesg,...) {\
    fprintf(stderr, \
            "[%s:%d]: " mesg "\n", \
            __FILE__, \
            __LINE__, \
            __VA_ARGS__);}


/*Macros that return error codes.*/
#define fatal(err,mesg,...) {\
    log_err(mesg, \
            __VA_ARGS__); \
    return err;}


#define check(val) {\
    int e_ = val; \
    if (e_ != SUCCESS) \
    { \
        backtrace(); \
        return e_; \
    }}


#define sentinel() {\
    fatal(SENTINEL_ERR, \
          "This branch should never be reached (%s,%d).", \
          __FILE__, \
          __LINE__)};


/*Safety checks.*/
#if defined(__CUDA_ARCH__) || defined(FAST)


#define not_null(p) {}
#define is_null(p) {}
#define not_nan(v) {}
#define min_check(v,min) {}
#define max_check(v,max) {}
#define in_range(v,min,max) {}


#else


#define not_null(p) {\
    if (p == NULL) \
    { \
        fatal(NULL_ERR, \
              "null pointer at address %p.",(void *)(&p)); \
    }}


#define is_null(p) {\
    if (p != NULL) \
    { \
        fatal(NON_NULL_ERR, \
              "pointer at address %p is not null.",(void *)(&p)); \
    }}


#define not_nan(v) {\
    if (isnan((double)v)) \
    { \
        fatal(INVALID_ERR, \
              "input value (%e) is Nan.", \
              (double)v); \
    }}


#define min_check(v,min) {\
    not_nan(v) \
    not_nan(min) \
    if (v < min) \
    { \
        fatal(RANGE_ERR, \
              "value (%e) less than minimum allowed (%e).", \
              (double)v, \
              (double)min); \
    }}


#define max_check(v,max) {\
    not_nan(v) \
    not_nan(max) \
    if (v > max) \
    { \
        fatal(RANGE_ERR, \
              "value (%e) greater than maximum allowed (%e).", \
              (double)v, \
              (double)max); \
    }}


#define in_range(v,min,max) {\
    if (min > max) \
    { \
        fatal(RANGE_ERR, \
              "min value (%e) greater tha max value (%e).", \
              (double)min, \
              (double)max); \
    } \
    min_check(v,min); \
    max_check(v,max);}


#endif


#ifdef __NVCC__
#define kernel_err(mesg,...) {}
#else
#define kernel_err(mesg,...) {log_err(mesg,__VA_ARGS__);exit(1);}
#endif


#endif
