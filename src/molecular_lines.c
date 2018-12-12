/** @file*/
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "debug.h"
#include "floating_point_type.h"
#include "host_launch.h"
#include "molecular_lines.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "spectral_bin.h"
#include "utils.h"
#include "verbosity.h"
#include "water_vapor_continuum.h"


static int const MIN_NUM_LEVELS = 2; /**< Minimum number of atmospheric
                                          levels allowed.*/
static int const MAX_NUM_LEVELS = 101; /**< Maximum number of atmospheric
                                            levels allowed.*/
static double const MIN_WAVENUMBER = 1.; /**< Smallest wavenumber [1/cm]
                                              allowed.*/
static double const MAX_WAVENUMBER = 3250.; /**< Larget wavenumber [1/cm]
                                                 allowed.*/
static double const MIN_RESOLUTION = 0.001; /**< Finest spectral resolution
                                                 [1/cm] allowed.*/
static double const MAX_RESOLUTION = 10.; /**< Coarsest spectral resolution
                                               [1/cm] allowed.*/
static double const MIN_CUTOFF = 1.; /**< Smallest cut-off [1/cm]
                                          from a line center allowed.*/
static double const MAX_CUTOFF = 50.; /**< Larget cut-off [1/cm]
                                           from a line center allowed.*/
static int const MAX_NUM_LINES = 524288; /**< Largest number of spectral
                                              lines per molecule allowed.*/
static double const DEFAULT_CUTOFF = 25.; /**< Default cut-off [1/cm] from
                                               a line center.*/
static int const DEFAULT_GPU = 0; /**< Default GPU device to use.*/


/** @brief Initialize a library context.
    @return SUCCESS or an error code.*/
EXTERN int grt_context_init(GrtContext_t **context,
                            int const num_levels,
                            double const w0,
                            double const wn,
                            double const wres,
                            char const * const hitran_path,
                            char const * const h2o_ctm_dir,
                            char const * const o3_ctm_dir,
                            double const * const wcutoff,
                            int const * const gpu_id,
                            int const * const num_threads,
                            int const * const optical_depth_method)
{
    /*Determine whether this context will be associated with a specific
      GPU or the host CPU.*/
    GrtContext_t c;
    int num_devices = 0;
    gpu_throw(get_num_gpus(&num_devices));
    if (gpu_id != NULL)
    {
        c.gpu_id = *gpu_id;
    }
    else if (num_devices > 0)
    {
        c.gpu_id = DEFAULT_GPU;
    }
    else
    {
        c.gpu_id = HOST_ONLY;
    }
    if (c.gpu_id == HOST_ONLY)
    {
        int const min_num_threads = 1;
        int const max_num_threads = omp_get_max_threads();
        if (num_threads != NULL)
        {
            in_range(*num_threads,min_num_threads,max_num_threads);
            c.num_threads = *num_threads;
        }
        else
        {
            c.num_threads = max_num_threads;
        }
        omp_set_num_threads(c.num_threads);
#ifdef _OPENMP
        log_mesg("Using %d OpenMP threads.",
                 c.num_threads);
#endif
    }
    else
    {
#ifdef __NVCC__
        in_range(c.gpu_id,0,num_devices);
        log_mesg("Using GPU device %d.",
                 c.gpu_id);
#else
        raise(COMPILER_ERR,
              "you must build with nvcc in order to use GPUs (gpu_id=%d.)",
              c.gpu_id);
#endif
    }

    /*Set the size of the atmospheric column.*/
    in_range(num_levels,MIN_NUM_LEVELS,MAX_NUM_LEVELS);
    c.num_levels = num_levels;
    c.num_layers = num_levels - 1;
    log_mesg("Atmospheric column properties:\n\tnumber of levels: %d\n\t"
                 "number of layers: %d",
             c.num_levels,
             c.num_layers);

    /*Set the spectral grid properties.*/
    in_range(w0,MIN_WAVENUMBER,MAX_WAVENUMBER);
    c.w0 = w0;
    in_range(wn,MIN_WAVENUMBER,MAX_WAVENUMBER);
    c.wn = wn;
    in_range(wres,MIN_RESOLUTION,MAX_RESOLUTION);
    c.wres = wres;
    if (wn <= w0)
    {
        raise(VALUE_ERR,
              "Spectral grid upper bound (%e) must be > than the spectral"
                  " grid lower bound (%e).",
              wn,
              w0);
    }
    c.num_wpoints = ceil((wn-w0)/wres) + 1.;
    log_mesg("Spectral grid properties:\n\tlower bound: %e [1/cm]\n\t"
                 "upper bound: %e [1/cm]\n\tresolution: %e [1/cm]\n\t"
                 "total size: %zu grid points",
             c.w0,
             c.wn,
             c.wres,
             c.num_wpoints);

    /*Create the spectral bins.*/
    double bin_width = 1.;
    throw(create_spectral_bins(&(c.bins),
                               c.num_layers,
                               c.w0,
                               c.num_wpoints,
                               c.wres,
                               bin_width,
                               c.gpu_id));
    log_mesg("Spectral bin properties:\n\tnumber of bins: %zu\n\t"
                 "bin width: %e\n\tspectral grid points per bin: %d\n\t"
                 "interpolation: %d\n\tspectral gid points in last bin:"
                 " %d\n\tinterpolation in last bin: %d",
             c.bins.n,
             bin_width,
             c.bins.ppb,
             c.bins.do_interp,
             c.bins.last_ppb,
             c.bins.do_last_interp);

    /*Store the path to the hitran database file.*/
    not_null(hitran_path);
    snprintf(c.hitran_path,
             DIR_PATH_LEN,
             "%s",
             hitran_path);
    log_mesg("Using HITRAN database file %s.",
             c.hitran_path);

    /*Set the molecular line cutoff.*/
    if (wcutoff != NULL)
    {
        in_range(*wcutoff,MIN_CUTOFF,MAX_CUTOFF);
        c.wcutoff = *wcutoff;
    }
    else
    {
        c.wcutoff = DEFAULT_CUTOFF;
    }
    log_mesg("Using spectral line cut-off of %e [1/cm].",
             c.wcutoff);


    /*Set the method that will be used to calculate the optical depths.*/
    if (optical_depth_method != NULL)
    {
        in_range(*optical_depth_method,wavenumber_sweep,line_sample);
        c.optical_depth_method = *optical_depth_method;
    }
    else
    {
        c.optical_depth_method = wavenumber_sweep;
    }

    /*Prepare to add molecules.*/
    c.num_molecules = 0;
    c.molecule_bit_field = 0;

    /*Pepare water vapor continuum.*/
    c.use_h2o_ctm = 0;
    if (h2o_ctm_dir != NULL)
    {
        if (strcmp(h2o_ctm_dir,"none") != 0)
        {
            c.use_h2o_ctm = 1;
            throw(copy_str(c.h2o_ctm_dir,
                           h2o_ctm_dir,
                           DIR_PATH_LEN));
        }
    }

    /*Prepare ozone continuum.*/
    c.use_o3_ctm = 0;
    if (o3_ctm_dir != NULL)
    {
        if (strcmp(o3_ctm_dir,"none") != 0)
        {
            c.use_o3_ctm = 1;
            throw(copy_str(c.o3_ctm_dir,
                           o3_ctm_dir,
                           DIR_PATH_LEN));
        }
    }

    /*Reserve memory.*/
    gmalloc(c.x,c.num_levels*NUM_MOLS,c.gpu_id);
    gmalloc(c.n,c.num_layers,c.gpu_id);
    gmalloc(c.pavg,c.num_layers,c.gpu_id);
    gmalloc(c.tavg,c.num_layers,c.gpu_id);
    gmalloc(c.psavg,c.num_layers,c.gpu_id);
    gmalloc(c.ns,c.num_layers,c.gpu_id);
    gmalloc(c.linecenter,c.num_layers*MAX_NUM_LINES,c.gpu_id);
    gmalloc(c.snn,c.num_layers*MAX_NUM_LINES,c.gpu_id);
    gmalloc(c.gamma,c.num_layers*MAX_NUM_LINES,c.gpu_id);
    gmalloc(c.alpha,c.num_layers*MAX_NUM_LINES,c.gpu_id);
    if (c.gpu_id != HOST_ONLY)
    {
        gmalloc(c.p,c.num_levels,c.gpu_id);
        gmalloc(c.t,c.num_levels,c.gpu_id);
        gmalloc(c.tau,c.num_layers*c.num_wpoints,c.gpu_id);
    }

    /*Copy data into the input context.*/
    not_null(context);
    throw(malloc_ptr((void **)context,
                     sizeof(**context)));
    not_null(*context);
    memcpy(*context,
           &c,
           sizeof(**context));
    return SUCCESS;
}


/** @brief Finalize a library context
    @return SUCCESS or an error code.*/
EXTERN int grt_context_free(GrtContext_t **context /**< Library context.*/
                           )
{
    not_null(context);
    GrtContext_t *c = *context;
    not_null(c);
    int i;
    for (i=0;i<c->num_molecules;++i)
    {
        throw(free_molecule(&(c->mols[i])));
    }
    throw(destroy_spectral_bins(&(c->bins)));
    gfree(c->x,c->gpu_id);
    gfree(c->n,c->gpu_id);
    gfree(c->pavg,c->gpu_id);
    gfree(c->tavg,c->gpu_id);
    gfree(c->psavg,c->gpu_id);
    gfree(c->ns,c->gpu_id);
    gfree(c->linecenter,c->gpu_id);
    gfree(c->snn,c->gpu_id);
    gfree(c->gamma,c->gpu_id);
    gfree(c->alpha,c->gpu_id);
    if (c->gpu_id != HOST_ONLY)
    {
        gfree(c->p,c->gpu_id);
        gfree(c->t,c->gpu_id);
        gfree(c->tau,c->gpu_id);
    }
    if (c->use_h2o_ctm && is_molecule_active(c->molecule_bit_field,H2O))
    {
        throw(free_water_vapor_continuum_coefs(&(c->h2o_cc)));
    }
    if (c->use_o3_ctm && is_molecule_active(c->molecule_bit_field,O3))
    {
        throw(free_ozone_continuum_coefs(&(c->o3_cc)));
    }
    free(c);
    *context = NULL;
    return SUCCESS;
}


/** @brief Add a molecule to a context
    @return SUCCESS or an error code.*/
EXTERN int grt_add_molecule(GrtContext_t *context, /**< Library context.*/
                            int const molecule_id, /**< Molecule id.*/
                            double const * const min_line_center, /**< Lower bound [1/cm]
                                                                       for spectral line
                                                                       centers.*/
                            double const * const max_line_center /**< Upper bound [1/cm]
                                                                      for spectral line
                                                                      centers.*/
                           )
{
    not_null(context);
    if (is_molecule_active(context->molecule_bit_field,molecule_id))
    {
        raise(VALUE_ERR,
              "molecule %d has already been added.",
              molecule_id);
    }
    int index = context->num_molecules;
    (context->num_molecules)++;
    in_range(context->num_molecules,1,NUM_MOLS);
    throw(activate_molecule(&(context->molecule_bit_field),
                            molecule_id));
    double w0;
    if (min_line_center != NULL)
    {
        in_range(*min_line_center,MIN_WAVENUMBER,MAX_WAVENUMBER);
        w0 = *min_line_center;
    }
    else
    {
        w0 = context->w0;
    }
    double wn;
    if (max_line_center != NULL)
    {
        in_range(*max_line_center,MIN_WAVENUMBER,MAX_WAVENUMBER);
        wn = *max_line_center;
    }
    else
    {
        wn = context->wn;
    }
    min_check(wn,w0);
    throw(molecule(&(context->mols[index]),
                   molecule_id,
                   context->hitran_path,
                   w0,
                   wn));
    log_mesg("Using %s (%zu lines in range %e - %e [1/cm]).",
             context->mols[index].name,
             context->mols[index].line_params.num_lines,
             w0,
             wn);

    if (molecule_id == H2O && context->use_h2o_ctm)
    {
        /*Read in the water vapor continuum coefficients.*/
        log_mesg("Using the %s continuum.",
                 context->mols[index].name);
        throw(get_water_vapor_continuum_coefs(&(context->h2o_cc),
                                              context->h2o_ctm_dir,
                                              context->num_wpoints,
                                              context->w0,
                                              context->wres,
                                              context->gpu_id));
    }

    if (molecule_id == O3 && context->use_o3_ctm)
    {
        /*Read in the ozone continuum coefficients.*/
        log_mesg("Using the %s continuum.",
                 context->mols[index].name);
        throw(get_ozone_continuum_coefs(&(context->o3_cc),
                                        context->o3_ctm_dir,
                                        context->num_wpoints,
                                        context->w0,
                                        context->wres,
                                        context->gpu_id));
    }
    return SUCCESS;
}


/*Update a molecule's ppmv.*/
EXTERN int grt_set_molecule_ppmv(GrtContext_t *context,
                                 int const molecule_id,
                                 fp_t const * const ppmv)
{
    not_null(context);
    not_null(ppmv);
    if (!is_molecule_active(context->molecule_bit_field,molecule_id))
    {
        log_warn("molecule %d is not being used.",
                 molecule_id);
        return SUCCESS;
    }
    int index;
    throw(molecule_hash(molecule_id,
                        &index));
    fp_t a[context->num_levels];
    int i;
    for (i=0;i<context->num_levels;++i)
    {
        a[i] = ppmv[i]*1.e-6;
    }
    int offset = index*context->num_levels;
    gmemcpy(&(context->x[offset]),
            a,
            context->num_levels,
            context->gpu_id,
            FROM_HOST);
    return SUCCESS;
}


/*Calcluate the total optical depth in each layer at each spectral grid
  point.*/
EXTERN int grt_calculate_optical_depth(GrtContext_t *context,
                                       fp_t const * const pressure,
                                       fp_t const * const temperature,
                                       fp_t *optical_depth)
{
    not_null(context);
    not_null(pressure);
    not_null(temperature);
    not_null(optical_depth);
    fp_t const mbtoatm = 0.000986923f;
    fp_t p[context->num_levels];
    int i;
    for (i=0;i<context->num_levels;++i)
    {
        p[i] = pressure[i]*mbtoatm;
    }
    if (context->gpu_id == HOST_ONLY)
    {
        log_mesg("Calculating optical depths for %d molecules in %d"
                     " atmospheric layers on the host CPU.",
                 context->num_molecules,
                 context->num_levels-1);
    }
    else
    {
        log_mesg("Calculating optical depths for %d molecules in %d"
                     " atmospheric layers on GPU %d.",
                 context->num_molecules,
                 context->num_levels-1,
                 context->gpu_id);
    }
    throw(launch_h(context,
                   p,
                   temperature,
                   optical_depth));
/*
#ifdef __NVCC__
        size_t num_elements = context->num_levels;
        HANDLE_ERROR(cudaMemcpy(context->P,
                                p,
                                sizeof(*p)*num_elements,
                                cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(context->T,
                                temperature,
                                sizeof(*temperature)*num_elements,
                                cudaMemcpyHostToDevice));
        throw(launch(context->num_levels,
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
                     context->molecule_bit_field,
                     context->line_params,
                     context->w0,
                     context->wres,
                     context->num_wpoints,
                     context->wcutoff,
                     context->use_h2o_ctm,
                     context->h2o_cc,
                     context->use_o3_ctm,
                     context->o3_cc,
                     context->tau,
                     context->fine_factor,
                     context->bins));
        num_elements = context->num_layers*context->num_wpoints;
        HANDLE_ERROR(cudaMemcpy(optical_depth,
                                context->tau,
                                sizeof(*optical_depth)*num_elements,
                                cudaMemcpyDeviceToHost));
#endif
*/
    return SUCCESS;
}


/*Get the number of molecules that have been added to the context.*/
EXTERN int grt_get_num_molecules(GrtContext_t const * const context,
                                 int * const n)
{
    not_null(context);
    not_null(n);
    *n = context->num_molecules;
    return SUCCESS;
}


/*Get the number of spectral grid points for the input context.*/
EXTERN int grt_get_spectral_grid_size(GrtContext_t const * const context,
                                      uint64_t * const n)
{
    not_null(context);
    not_null(n);
    *n = context->num_wpoints;
    return SUCCESS;
}


/*Return a message for an input return code.*/
EXTERN int grt_errstr(int const code,
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


/*Set the verbosity level for the library.*/
EXTERN void grt_set_verbosity(int const level)
{
    set_verbosity(level);
}


/*Get the verbosity level for the library.*/
EXTERN int grt_get_verbosity(void)
{
    return get_verbosity();
}
