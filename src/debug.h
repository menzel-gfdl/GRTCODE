#ifndef DEBUG_H_
#define DEBUG_H_

#include <stdarg.h>
#include <stdio.h>
#include "netcdf.h"

enum return_codes
{
    SUCCESS = 0,
    ERR,
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

#define fatal(mesg,...) {\
    log_err(mesg, \
            __VA_ARGS__); \
    return ERR;}

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
        fatal("null pointer at address %p.",&p); \
    }}

#define is_null(p) {\
    if (p != NULL) \
    { \
        fatal("pointer at address %p is not null.",&p); \
    }}

#define netcdf_check(val) {\
    int e_ = val; \
    if (e_ != NC_NOERR) \
    { \
        fatal("netcdf returned error code %d. %s.", \
              e_, \
              nc_strerror(e_)); \
    }}

#endif
