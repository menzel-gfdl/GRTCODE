#ifndef DEBUG_H_
#define DEBUG_H_

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "utils.h"
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
    IO_ERR,
    GPU_ERR
};


#define backtrace() { \
    if (get_verbosity() >= ERROR) { \
        fprintf(stderr, \
                "\r\33[2K\t%s: %d\n", \
                __FILE__, \
                __LINE__); \
    }}


#define log_err(mesg,...) { \
    if (get_verbosity() >= ERROR) { \
        fprintf(stderr, \
                "\r\33[2K[%s] error: " mesg "\nBacktrace:\n", \
                __func__, \
                __VA_ARGS__); \
        backtrace(); \
    }}


#define log_warn(mesg,...) { \
    if (get_verbosity() >= WARN) { \
        fprintf(stderr, \
                "\r\33[2K[%s:%d] warning: " mesg "\n", \
                __FILE__, \
                __LINE__, \
                __VA_ARGS__); \
    }}


#define log_info(mesg,...) { \
    if (get_verbosity() >= INFO) { \
        fprintf(stderr, \
                "\r\33[2K[%s:%d] info: " mesg "\n", \
                __FILE__, \
                __LINE__, \
                __VA_ARGS__); \
    }}


#define log_mesg(mesg,...) { \
    if (get_verbosity() >= NONE) { \
        fprintf(stdout, \
                "\r\33[2K" mesg "\n", \
                __VA_ARGS__); \
    }}


/*Macros that return error codes.*/
#define raise(err,mesg,...) { \
    log_err(mesg, \
            __VA_ARGS__); \
    return err;}


#define throw(val) { \
    int e_ = val; \
    if (e_ != SUCCESS) \
    { \
        backtrace(); \
        return e_; \
    }}


#define sentinel() { \
    raise(SENTINEL_ERR, \
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
#define not_null(p) { \
    if (p == NULL) \
    { \
        raise(NULL_ERR, \
              "null pointer at address %p.",(void *)(&p)); \
    }}


#define is_null(p) { \
    if (p != NULL) \
    { \
        raise(NON_NULL_ERR, \
              "pointer at address %p is not null.",(void *)(&p)); \
    }}


#define not_nan(v) { \
    if (isnan((double)v)) \
    { \
        raise(INVALID_ERR, \
              "input value (%e) is Nan.", \
              (double)v); \
    }}


#define min_check(v,min) { \
    not_nan(v) \
    not_nan(min) \
    if (v < min) \
    { \
        raise(RANGE_ERR, \
              "value (%e) less than minimum allowed (%e).", \
              (double)v, \
              (double)min); \
    }}


#define max_check(v,max) { \
    not_nan(v) \
    not_nan(max) \
    if (v > max) \
    { \
        raise(RANGE_ERR, \
              "value (%e) greater than maximum allowed (%e).", \
              (double)v, \
              (double)max); \
    }}


#define in_range(v,min,max) { \
    if (min > max) \
    { \
        raise(RANGE_ERR, \
              "min value (%e) greater tha max value (%e).", \
              (double)min, \
              (double)max); \
    } \
    min_check(v,min); \
    max_check(v,max);}
#endif


#define HOST_ONLY -1
#define cat(a,b) a##b


#ifdef __NVCC__
#define gpu_throw(val) {\
    cudaError_t e_ = val; \
    if (e_ != cudaSuccess) \
    { \
        raise(GPU_ERR, \
              "cuda: %s", \
              cudaGetErrorString(e_)); \
    }}
#define _glaunch(func,threads,loc,...) { \
    gpu_throw(cudaSetDevice(loc)); \
    int min_grid_size; \
    int dim_block; \
    gpu_throw(cudaOccupancyMaxPotentialBlockSize(&min_grid_size, \
                                                 &dim_block, \
                                                 cat(func,_d), \
                                                 0, \
                                                 (int)threads)); \
    int dim_grid = (((int)threads) + dim_block - 1)/dim_block; \
    cat(func,_d)<<<dim_grid,dim_block,0,0>>>(__VA_ARGS__);}
#define FROM_HOST cudaMemcpyHostToDevice
#define FROM_DEVICE cudaMemcpyDeviceToHost
#define HOST __host__
#define DEVICE __device__
#else
#define gpu_throw(val) {}
#define _glaunch(func,threads,loc,...) {}
#define FROM_HOST
#define FROM_DEVICE
#define HOST
#define DEVICE
#endif


#ifndef _OPENMP
#define omp_set_num_threads(n) {}
#define omp_get_max_threads() 1
#endif


#define gmalloc(ptr,size,loc) { \
    if (loc == HOST_ONLY) \
    { \
        throw(malloc_ptr((void **)&ptr, \
                         sizeof(*ptr)*size)); \
    } \
    else \
    { \
        gpu_throw(cudaSetDevice(loc)); \
        gpu_throw(cudaMalloc(&ptr, \
                             sizeof(*ptr)*size)); \
    }}


#define gfree(ptr,loc) { \
    if (loc == HOST_ONLY) \
    { \
        throw(free_ptr((void **)&ptr)); \
    } \
    else \
    { \
        gpu_throw(cudaSetDevice(loc)); \
        gpu_throw(cudaFree(ptr)); \
    }}


#define gmemset(ptr,val,size,loc) { \
    if (loc == HOST_ONLY) \
    { \
        memset(ptr, \
               val, \
               sizeof(*ptr)*size); \
    } \
    else \
    { \
        gpu_throw(cudaSetDevice(loc)); \
        gpu_throw(cudaMemset(ptr, \
                             val, \
                             sizeof(*ptr)*size)); \
    }}


#define gmemcpy(dst,src,size,loc,dir) { \
    if (loc == HOST_ONLY) \
    { \
        memcpy(dst, \
               src, \
               sizeof(*dst)*size); \
    } \
    else \
    { \
        gpu_throw(cudaSetDevice(loc)); \
        gpu_throw(cudaMemcpy(dst, \
                             src, \
                             sizeof(*dst)*size, \
                             dir)); \
    }}


#define glaunch(func,threads,loc,...) { \
    if (loc == HOST_ONLY) \
    { \
        int t_ = omp_get_max_threads(); \
        if ((uint64_t)threads < (uint64_t)t_) \
        { \
            omp_set_num_threads(threads); \
        } \
        throw(func(__VA_ARGS__)); \
        if ((uint64_t)threads < (uint64_t)t_) \
        { \
            omp_set_num_threads(t_); \
        } \
    } \
    else \
    { \
        _glaunch(func,threads,loc,__VA_ARGS__); \
    }}


#endif
