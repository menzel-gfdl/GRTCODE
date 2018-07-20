#ifndef DEBUG_H_
#define DEBUG_H_

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "verbosity.h"


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
};


#define backtrace() { \
    if (get_verbosity() <= ERROR) { \
        fprintf(stderr, \
                "\t%s: %d\n", \
                __FILE__, \
                __LINE__); \
    }}


#define log_err(mesg,...) { \
    if (get_verbosity() <= ERROR) { \
        fprintf(stderr, \
                "[%s] error: " mesg "\nBacktrace:\n", \
                __func__, \
                __VA_ARGS__); \
        backtrace(); \
    }}


#define log_warn(mesg,...) { \
    if (get_verbosity() <= WARN) { \
        fprintf(stderr, \
                "[%s:%d] warning: " mesg "\n", \
                __FILE__, \
                __LINE__, \
                __VA_ARGS__); \
    }}


#define log_info(mesg,...) {\
    if (get_verbosity() <= INFO) { \
        fprintf(stderr, \
                "[%s:%d] info: " mesg "\n", \
                __FILE__, \
                __LINE__, \
                __VA_ARGS__); \
    }}


#define log_mesg(mesg,...) {\
        fprintf(stdout, \
                mesg "\n", \
                __VA_ARGS__); \
    }


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
