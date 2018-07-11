#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __NVCC__
#include "cuda_helpers.cuh"
#endif
#include "debug.h"
#include "device_launch.h"
#include "host_launch.h"
#include "new.h"
#include "ozone_continuum.h"
#include "parse_HITRAN_file.h"
#ifdef __NVCC__
#include "query_gpu.cuh"
#endif
#include "TIPS_2011.h"
#include "utils.h"
#include "water_vapor_continuum.h"


/*Bounds for input parameters.*/
int const MIN_NUM_LEVELS = 2;
int const MAX_NUM_LEVELS = 61;
int const MIN_NUM_MOLECULES = 1;
int const MAX_NUM_MOLECULES = 10;
double const MIN_WAVENUMBER = 1.;
double const MAX_WAVENUMBER = 50000.;
double const MIN_RESOLUTION = 0.01;
double const MAX_RESOLUTION = 10.;
double const MIN_CUTOFF = 1.;
double const MAX_CUTOFF = 50.;
int const MAX_NUM_LINES = 524288; /*2^19*/


/*Optional parameter default values.*/
double const DEFAULT_CUTOFF = 25.;
int const DEFAULT_GPU = 0;
int const HOST_ONLY = -1;
int const DEFAULT_H2O_CONTINUUM = 1;
int const DEFAULT_O3_CONTINUUM = 1;


/*Initialize the library.*/
#ifdef __NVCC__
extern "C"
#endif
int grt_context_init(GrtContext_t **context,
                     int const num_levels,
                     double const w0,
                     double const wn,
                     double const wres,
                     uint64_t * const num_wpoints,
                     double const * const wcutoff,
                     int const * const gpu_id,
                     int const * const num_threads,
                     char const * const h2o_ctm_dir,
                     char const * const o3_ctm_dir)
{
    /*Guard against bad constants.*/
    assert(MIN_NUM_LEVELS >= 2);
    assert(MAX_NUM_LEVELS >= MIN_NUM_LEVELS);
    assert(MIN_NUM_MOLECULES >= 1);
    assert(MAX_NUM_MOLECULES >= MIN_NUM_MOLECULES);
    assert(MIN_WAVENUMBER >= 1.);
    assert(MAX_WAVENUMBER >= MIN_WAVENUMBER);
    assert(MIN_RESOLUTION > 0.);
    assert(MAX_RESOLUTION >= MIN_RESOLUTION);
    assert((MAX_WAVENUMBER-MIN_WAVENUMBER)/MIN_RESOLUTION + 1. <=
           (double)UINT64_MAX);
    assert(MIN_CUTOFF >= 1.);
    assert(MAX_CUTOFF >= MIN_CUTOFF);

    /*Define the problem size.*/
    GrtContext_t c;
    in_range(num_levels,MIN_NUM_LEVELS,MAX_NUM_LEVELS);
    c.num_levels = num_levels;
    c.num_layers = num_levels - 1;
    in_range(w0,MIN_WAVENUMBER,MAX_WAVENUMBER);
    c.w0 = w0;
    in_range(wn,MIN_WAVENUMBER,MAX_WAVENUMBER);
    c.wn = wn;
    in_range(wres,MIN_RESOLUTION,MAX_RESOLUTION);
    c.wres = wres;
    if (wn <= w0)
    {
        fatal(VALUE_ERR,
              "Spectral grid upper bound (%e) must be > than the spectral"
                  " grid lower bound (%e).",
              wn,
              w0);
    }
    c.num_wpoints = (wn-w0)/wres + 1.;
    not_null(num_wpoints);
    *num_wpoints = c.num_wpoints;

    /*Set optional parameters.*/
    if (wcutoff != NULL)
    {
        in_range(*wcutoff,MIN_CUTOFF,MAX_CUTOFF);
        c.wcutoff = *wcutoff;
    }
    else
    {
        c.wcutoff = DEFAULT_CUTOFF;
    }

#ifdef __NVCC__
    int num_devices;
    check(get_num_gpus(&num_devices,
                       1));
#else
    int num_devices = 0;
#endif

    if (gpu_id != NULL)
    {
        if (*gpu_id != HOST_ONLY)
        {
            in_range(*gpu_id,0,num_devices);
        }
        c.gpu_id = *gpu_id;
    }
    else
    {
        if (num_devices > 0)
        {
            c.gpu_id = DEFAULT_GPU;
        }
        else
        {
            c.gpu_id = HOST_ONLY;
        }
    }
    if (c.gpu_id != HOST_ONLY)
    {
#ifndef __NVCC__
        fatal(COMPILER_ERR,
              "you must build with nvcc in order to use GPUs (gpu_id=%d.)",
              c.gpu_id);
#else
        HANDLE_ERROR(cudaSetDevice(c.gpu_id));
        log_mesg("Molecular lines will be calculated using GPU device %d.",
                 c.gpu_id);
#endif
    }

    int const min_num_threads = 1;
#ifdef _OPENMP
    int const max_num_threads = omp_get_max_threads();
#else
    int const max_num_threads = min_num_threads;
#endif

    if (num_threads != NULL)
    {
        in_range(*num_threads,min_num_threads,max_num_threads);
        c.num_threads = *num_threads;
    }
    else
    {
        c.num_threads = max_num_threads;
    }

    /*Pepare water vapor continuum.*/
    if (h2o_ctm_dir != NULL)
    {
        c.use_h2o_ctm = 1;
    }
    else
    {
        c.use_h2o_ctm = 0;
    }
    if (c.use_h2o_ctm)
    {
        /*Read in the water vapor continuum coefficients.*/
        WaterVaporContinuumCoefs_t h2o_cc;
        check(get_water_vapor_continuum_coefs(&h2o_cc,
                                              h2o_ctm_dir,
                                              c.num_wpoints,
                                              c.w0,
                                              c.wres));
        check(malloc_ptr((void **) (&c.h2o_cc),
                         sizeof(*(c.h2o_cc))));
        if (c.gpu_id != HOST_ONLY)
        {
            check(put_water_vapor_coefs_on_device(&h2o_cc,
                                                  c.h2o_cc));
        }
        else
        {
            memcpy(c.h2o_cc,
                   &h2o_cc,
                   sizeof(*(c.h2o_cc)));
        }
    }

    /*Prepare ozone continuum.*/
    if (o3_ctm_dir != NULL)
    {
        c.use_o3_ctm = 1;
    }
    else
    {
        c.use_o3_ctm = 0;
    }
    if (c.use_o3_ctm)
    {
        /*Read in the ozone continuum coefficients.*/
        OzoneContinuumCoefs_t o3_cc;
        check(get_ozone_continuum_coefs(&o3_cc,
                                        o3_ctm_dir,
                                        c.num_wpoints,
                                        c.w0,
                                        c.wres));
        check(malloc_ptr((void **) (&c.o3_cc),
                         sizeof(*(c.o3_cc))));
        if (c.gpu_id != HOST_ONLY)
        {
            check(put_ozone_coefs_on_device(&o3_cc,
                                            c.o3_cc));
        }
        else
        {
            memcpy(c.o3_cc,
                   &o3_cc,
                   sizeof(*(c.o3_cc)));
        }
    }

    /*Reserve memory.*/
    c.num_molecules = 0;
    check(malloc_ptr((void **)(&(c.line_params)),
                     sizeof(*(c.line_params))*MAX_NUM_MOLECULES));
    not_null(c.line_params);
    if (c.gpu_id != HOST_ONLY)
    {
#ifdef __NVCC__
        size_t num_elements = c.num_levels;
        HANDLE_ERROR(cudaMalloc(&(c.P),
                                sizeof(*(c.P))*num_elements));
        HANDLE_ERROR(cudaMalloc(&(c.T),
                                sizeof(*(c.T))*num_elements));
        num_elements *= sizeof(*(c.x))*MAX_NUM_MOLECULES;
        HANDLE_ERROR(cudaMalloc(&(c.x),
                                num_elements));
        HANDLE_ERROR(cudaMemset(c.x,
                                0,
                                num_elements));
        num_elements = c.num_layers;
        HANDLE_ERROR(cudaMalloc(&(c.Pavg),
                                sizeof(*(c.Pavg))*num_elements));
        HANDLE_ERROR(cudaMalloc(&(c.Tavg),
                                sizeof(*(c.Tavg))*num_elements));
        HANDLE_ERROR(cudaMalloc(&(c.N),
                                sizeof(*(c.N))*num_elements));
        HANDLE_ERROR(cudaMalloc(&(c.Ns),
                                sizeof(*(c.Ns))*num_elements));
        HANDLE_ERROR(cudaMalloc(&(c.Psavg),
                                sizeof(*(c.Psavg))*num_elements));
        num_elements *= MAX_NUM_LINES;
        HANDLE_ERROR(cudaMalloc(&(c.gamma),
                                sizeof(*(c.gamma))*num_elements));
        HANDLE_ERROR(cudaMalloc(&(c.Pshift),
                                sizeof(*(c.Pshift))*num_elements));
        HANDLE_ERROR(cudaMalloc(&(c.s),
                                sizeof(*(c.s))*num_elements));
        num_elements = c.num_layers*c.num_wpoints;
        HANDLE_ERROR(cudaMalloc(&(c.tau),
                                sizeof(*(c.tau))*num_elements));
#endif
        check(alloc_line_params_device(&(c.lines),
                                       MAX_NUM_LINES));
        check(initTIPS_d());
    }
    else
    {
        check(malloc_ptr((void **)(&c.x),
                         sizeof(*(c.x))*c.num_levels*MAX_NUM_MOLECULES));
        check(malloc_ptr((void **)(&c.Pavg),
                         sizeof(*(c.Pavg))*c.num_layers));
        check(malloc_ptr((void **)(&c.Tavg),
                         sizeof(*(c.Tavg))*c.num_layers));
        check(malloc_ptr((void **)(&c.N),
                         sizeof(*(c.N))*c.num_layers));
        check(malloc_ptr((void **)(&c.Ns),
                         sizeof(*(c.Ns))*c.num_layers));
        check(malloc_ptr((void **)(&c.Psavg),
                         sizeof(*(c.Psavg))*c.num_layers));
        check(malloc_ptr((void **)(&c.snn_ref),
                         sizeof(*(c.snn_ref))*c.num_layers*MAX_NUM_LINES));
        check(malloc_ptr((void **)(&c.gamma),
                         sizeof(*(c.gamma))*c.num_layers*MAX_NUM_LINES));
        check(malloc_ptr((void **)(&c.Pshift),
                         sizeof(*(c.Pshift))*c.num_layers*MAX_NUM_LINES));
        check(malloc_ptr((void **)(&c.s),
                         sizeof(*(c.s))*c.num_layers*MAX_NUM_LINES));
    }

    /*Copy data into the input context.*/
    not_null(context);
    check(malloc_ptr((void **)context,
                     sizeof(**context)));
    not_null(*context);
    memcpy(*context,
           &c,
           sizeof(**context));
    return SUCCESS;
}


/*Finalize the library.*/
#ifdef __NVCC__
extern "C"
#endif
int grt_context_free(GrtContext_t **context)
{
    not_null(context);
    GrtContext_t *c = *context;
    not_null(c);
    LineFlags_t flags = {((unsigned int) -1),1,0};
    int i;
    for (i=0;i<(c->num_molecules);++i)
    {
        check(free_line_params_host(&(c->line_params[i]),
                                    flags));
    }
    free(c->line_params);
    if (c->gpu_id != HOST_ONLY)
    {
#ifdef __NVCC__
        HANDLE_ERROR(cudaSetDevice(c->gpu_id));
        HANDLE_ERROR(cudaFree(c->P));
        HANDLE_ERROR(cudaFree(c->T));
        HANDLE_ERROR(cudaFree(c->x));
        HANDLE_ERROR(cudaFree(c->Pavg));
        HANDLE_ERROR(cudaFree(c->Tavg));
        HANDLE_ERROR(cudaFree(c->N));
        HANDLE_ERROR(cudaFree(c->Ns));
        HANDLE_ERROR(cudaFree(c->Psavg));
        HANDLE_ERROR(cudaFree(c->gamma));
        HANDLE_ERROR(cudaFree(c->Pshift));
        HANDLE_ERROR(cudaFree(c->s));
        HANDLE_ERROR(cudaFree(c->tau));
        check(free_line_params_device(&(c->lines)));
#endif
    }
    else
    {
        free(c->x);
        free(c->Pavg);
        free(c->Tavg);
        free(c->N);
        free(c->Ns);
        free(c->Psavg);
        free(c->snn_ref);
        free(c->gamma);
        free(c->Pshift);
        free(c->s);
    }
    if (c->use_h2o_ctm)
    {
        if (c->gpu_id != HOST_ONLY)
        {
            check(remove_water_vapor_coefs_from_device(c->h2o_cc));
        }
        else
        {
            check(free_water_vapor_continuum_coefs(c->h2o_cc));
        }
        free(c->h2o_cc);
    }
    if (c->use_o3_ctm)
    {
        if (c->gpu_id != HOST_ONLY)
        {
            check(remove_ozone_coefs_from_device(c->o3_cc));
        }
        else
        {
            check(free_ozone_continuum_coefs(c->o3_cc));
        }
        free(c->o3_cc);
    }
    free(c);
    *context = NULL;
    return SUCCESS;
}


/*Add a molecule.*/
#ifdef __NVCC__
extern "C"
#endif
int add_molecule(GrtContext_t *context,
                 char const * const hitran_filepath,
                 int * const molecule_id,
                 double const * const min_line_center_wavenumber,
                 double const * const max_line_center_wavenumber)
{
    not_null(context);
    not_null(hitran_filepath);
    int index = (context->num_molecules)++;
    in_range(context->num_molecules,1,MAX_NUM_MOLECULES);
    LineFlags_t flags = {((unsigned int) -1),1,0};
    double w0;
    if (min_line_center_wavenumber != NULL)
    {
        in_range(*min_line_center_wavenumber,MIN_WAVENUMBER,MAX_WAVENUMBER);
        w0 = *min_line_center_wavenumber;
    }
    else
    {
        w0 = context->w0;
    }
    double wn;
    if (max_line_center_wavenumber != NULL)
    {
        in_range(*max_line_center_wavenumber,MIN_WAVENUMBER,MAX_WAVENUMBER);
        wn = *max_line_center_wavenumber;
    }
    else
    {
        wn = context->wn;
    }
    min_check(wn,w0);
    check(parse_hitran_file(&(context->line_params[index]),
                            hitran_filepath,
                            flags,
                            w0,
                            wn));
    not_null(molecule_id);
    *molecule_id = index;
    return SUCCESS;
}


/*Update a molecule's ppmv.*/
#ifdef __NVCC__
extern "C"
#endif
int set_molecule_ppmv(GrtContext_t *context,
                      int const molecule_id,
                      fp_t const * const ppmv)
{
    not_null(context);
    not_null(ppmv);
    in_range(molecule_id,0,context->num_molecules-1);
    int offset = molecule_id*context->num_levels;
    size_t num_bytes = sizeof(*ppmv)*context->num_levels;
    if (context->gpu_id != HOST_ONLY)
    {
#ifdef __NVCC__
        HANDLE_ERROR(cudaSetDevice(context->gpu_id));
        HANDLE_ERROR(cudaMemcpy(&(context->x[offset]),
                                ppmv,
                                num_bytes,
                                cudaMemcpyHostToDevice));
#endif
    }
    else
    {
        memcpy(&(context->x[offset]),
               ppmv,
               num_bytes);
    }
    return SUCCESS;
}


/*Calcluate the total optical depth in each layer at each spectral grid
  point.*/
#ifdef __NVCC__
extern "C"
#endif
int calculate_optical_depth(GrtContext_t *context,
                            fp_t const * const pressure,
                            fp_t const * const temperature,
                            fp_t *optical_depth)
{
    not_null(context);
    not_null(pressure);
    not_null(temperature);
    not_null(optical_depth);
    if (context->gpu_id != HOST_ONLY)
    {
#ifdef __NVCC__
        HANDLE_ERROR(cudaSetDevice(context->gpu_id));
        size_t num_elements = context->num_levels;
        HANDLE_ERROR(cudaMemcpy(context->P,
                                pressure,
                                sizeof(*pressure)*num_elements,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(context->T,
                                temperature,
                                sizeof(*temperature)*num_elements,
                                cudaMemcpyHostToDevice));
        check(launch(context->num_levels,
                     context->P,
                     context->T,
                     context->x,
                     context->Pavg,
                     context->Tavg,
                     context->N,
                     context->Ns,
                     context->Psavg,
                     context->gamma,
                     context->Pshift,
                     context->s,
                     context->lines,
                     context->num_molecules,
                     context->line_params,
                     context->w0,
                     context->wres,
                     context->num_wpoints,
                     context->wcutoff,
                     context->use_h2o_ctm,
                     context->h2o_cc,
                     context->use_o3_ctm,
                     context->o3_cc,
                     context->tau));
        num_elements = context->num_layers*context->num_wpoints;
        HANDLE_ERROR(cudaMemcpy(optical_depth,
                                context->tau,
                                sizeof(*optical_depth)*num_elements,
                                cudaMemcpyDeviceToHost));
#endif
    }
    else
    {
        check(launch_h(context->num_levels,
                       pressure,
                       temperature,
                       context->x,
                       context->Pavg,
                       context->Tavg,
                       context->N,
                       context->Ns,
                       context->Psavg,
                       context->snn_ref,
                       context->gamma,
                       context->Pshift,
                       context->s,
                       context->lines,
                       context->num_molecules,
                       context->line_params,
                       context->w0,
                       context->wres,
                       context->num_wpoints,
                       context->wcutoff,
                       context->use_h2o_ctm,
                       context->h2o_cc,
                       context->use_o3_ctm,
                       context->o3_cc,
                       optical_depth));
    }
    return SUCCESS;
}


int grt_errstr(int const code,
               char * const buf,
               int const buf_size)
{
    not_null(buf);
    min_check(buf_size,1);
    switch (code)
    {
        case SUCCESS:
            break;
        case INVALID_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: detected a floating point invalid.");
            break;
        case DIVBYZERO_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: detected a floating point divide-by-zero.");
            break;
        case OVERFLOW_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: detected a floating point overflow.");
            break;
        case UNDERFLOW_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: detected a floating point underflow.");
            break;
        case SENTINEL_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: entered unexpected code branch.");
            break;
        case NULL_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: attempt to dereference a null pointer.");
            break;
        case NON_NULL_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: expected a null pointer, but pointer already"
                         " has a value assigned to it.");
            break;
        case RANGE_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: detected a value out of its expected range.");
            break;
        case VALUE_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: detected a bad value.");
            break;
        case COMPILER_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: build was done with an incorrect compiler.");
            break;
        case IO_ERR:
            snprintf(buf,
                     buf_size,
                     "GRT: error while performing I/O.");
            break;
        default:
            snprintf(buf,
                     buf_size,
                     "Unknown code %d.",
                     code);
    }
    return SUCCESS;
}
