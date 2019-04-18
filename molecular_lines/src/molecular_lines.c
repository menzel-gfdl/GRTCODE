/** @file*/
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "cfcs.h"
#include "debug.h"
#include "floating_point_type.h"
#include "launch.h"
#include "molecular_lines.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "rs_config.h"
#include "spectral_bin.h"
#include "tips2017.h"
#include "utils.h"
#include "verbosity.h"
#include "water_vapor_continuum.h"


static double const MIN_CUTOFF = 1.; /**< Smallest cut-off [1/cm]
                                          from a line center allowed.*/
static double const MAX_CUTOFF = 50.; /**< Larget cut-off [1/cm]
                                           from a line center allowed.*/
static int const MAX_NUM_LINES = 524288; /**< Largest number of spectral
                                              lines per molecule allowed.*/
static double const DEFAULT_CUTOFF = 25.; /**< Default cut-off [1/cm] from
                                               a line center.*/


/** @brief Initialize a library context.
    @return RS_SUCCESS or an error code.*/
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
    catch(get_num_gpus(&num_devices,
                       1));
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
        char *mesg = "Using %d OpenMP threads.";
        log_mesg(mesg, c.num_threads);
#endif
    }
    else
    {
#ifdef __NVCC__
        in_range(c.gpu_id, 0, num_devices);
        char *mesg = "Using GPU device %d.";
        log_mesg(mesg, c.gpu_id);
#endif
    }

    /*Set the size of the atmospheric column.*/
    in_range(num_levels,MIN_NUM_LEVELS,MAX_NUM_LEVELS);
    c.num_levels = num_levels;
    c.num_layers = num_levels - 1;
    char *mesg = "Atmospheric column properties:\n\tnumber of levels: %d\n\t"
                 "number of layers: %d";
    log_mesg(mesg, c.num_levels, c.num_layers);

    /*Set the spectral grid properties.*/
    in_range(w0, MIN_WAVENUMBER, MAX_WAVENUMBER);
    c.w0 = w0;
    in_range(wn, MIN_WAVENUMBER, MAX_WAVENUMBER);
    c.wn = wn;
    in_range(wres, MIN_RESOLUTION, MAX_RESOLUTION);
    c.wres = wres;
    if (wn <= w0)
    {
        mesg = "Spectral grid upper bound (%e) must be > than the spectral"
               " grid lower bound (%e).";
        raise(RS_VALUE_ERR, mesg, wn, w0);
    }
    c.num_wpoints = ceil((wn-w0)/wres) + 1.;
    mesg = "Spectral grid properties:\n\tlower bound: %e [1/cm]\n\t"
           "upper bound: %e [1/cm]\n\tresolution: %e [1/cm]\n\t"
           "total size: %zu grid points";
    log_mesg(mesg, c.w0, c.wn, c.wres, c.num_wpoints);

    /*Create the spectral bins.*/
    double bin_width = 1.;
    catch(create_spectral_bins(&(c.bins), c.num_layers, c.w0, c.num_wpoints, c.wres,
                               bin_width, c.gpu_id));
    mesg = "Spectral bin properties:\n\tnumber of bins: %zu\n\t"
           "bin width: %e\n\tspectral grid points per bin: %d\n\t"
           "interpolation: %d\n\tspectral gid points in last bin:"
           " %d\n\tinterpolation in last bin: %d";
    log_mesg(mesg, c.bins.n, bin_width, c.bins.ppb, c.bins.do_interp, c.bins.last_ppb,
             c.bins.do_last_interp);

    /*Store the path to the hitran database file.*/
    not_null(hitran_path);
    snprintf(c.hitran_path, DIR_PATH_LEN, "%s", hitran_path);
    mesg = "Using HITRAN database file %s.";
    log_mesg(mesg, c.hitran_path);

    /*Set the molecular line cutoff.*/
    if (wcutoff != NULL)
    {
        in_range(*wcutoff, MIN_CUTOFF, MAX_CUTOFF);
        c.wcutoff = *wcutoff;
    }
    else
    {
        c.wcutoff = DEFAULT_CUTOFF;
    }
    mesg = "Using spectral line cut-off of %e [1/cm].";
    log_mesg(mesg, c.wcutoff);


    /*Set the method that will be used to calculate the optical depths.*/
    if (optical_depth_method != NULL)
    {
        in_range(*optical_depth_method, wavenumber_sweep, line_sample);
        c.optical_depth_method = *optical_depth_method;
    }
    else
    {
        c.optical_depth_method = wavenumber_sweep;
    }

    /*Prepare to add molecules/cfcs.*/
    c.num_molecules = 0;
    c.molecule_bit_field = 0;
    c.num_cfcs = 0;
    c.cfc_bit_field = 0;

    /*Pepare water vapor continuum.*/
    c.use_h2o_ctm = 0;
    if (h2o_ctm_dir != NULL)
    {
        if (strcmp(h2o_ctm_dir, "none") != 0)
        {
            c.use_h2o_ctm = 1;
            catch(copy_str(c.h2o_ctm_dir, h2o_ctm_dir, DIR_PATH_LEN));
        }
    }

    /*Prepare ozone continuum.*/
    c.use_o3_ctm = 0;
    if (o3_ctm_dir != NULL)
    {
        if (strcmp(o3_ctm_dir, "none") != 0)
        {
            c.use_o3_ctm = 1;
            catch(copy_str(c.o3_ctm_dir, o3_ctm_dir, DIR_PATH_LEN));
        }
    }

    /*Reserve memory.*/
    gmalloc(c.x, c.num_levels*NUM_MOLS, c.gpu_id);
    gmalloc(c.x_cfc, c.num_levels*NUM_CFCS, c.gpu_id);
    gmalloc(c.n, c.num_layers, c.gpu_id);
    gmalloc(c.pavg, c.num_layers, c.gpu_id);
    gmalloc(c.tavg, c.num_layers, c.gpu_id);
    gmalloc(c.psavg, c.num_layers, c.gpu_id);
    gmalloc(c.ns, c.num_layers, c.gpu_id);
    gmalloc(c.linecenter, c.num_layers*MAX_NUM_LINES, c.gpu_id);
    gmalloc(c.snn, c.num_layers*MAX_NUM_LINES, c.gpu_id);
    gmalloc(c.gamma, c.num_layers*MAX_NUM_LINES, c.gpu_id);
    gmalloc(c.alpha, c.num_layers*MAX_NUM_LINES, c.gpu_id);
    if (c.gpu_id != HOST_ONLY)
    {
        gmalloc(c.p, c.num_levels, c.gpu_id);
        gmalloc(c.t, c.num_levels, c.gpu_id);
        gmalloc(c.tau, c.num_layers*c.num_wpoints, c.gpu_id);
    }

    /*Initialize TIPS.*/
    if (c.gpu_id != HOST_ONLY)
    {
        catch(inittips_d());
    }

    /*Copy data into the input context.*/
    not_null(context);
    catch(malloc_ptr((void **)context, sizeof(**context)));
    not_null(*context);
    memcpy(*context, &c, sizeof(**context));
    return RS_SUCCESS;
}


/** @brief Finalize a library context
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_context_free(GrtContext_t **context /**< Library context.*/
                           )
{
    not_null(context);
    GrtContext_t *c = *context;
    not_null(c);
    int i;
    for (i=0; i<c->num_molecules; ++i)
    {
        catch(free_molecule(&(c->mols[i])));
    }
    for (i=0; i<c->num_cfcs; ++i)
    {
        catch(free_cfc_cross_sections(&(c->cfcs[i])));
    }
    catch(destroy_spectral_bins(&(c->bins)));
    gfree(c->x, c->gpu_id);
    gfree(c->x_cfc, c->gpu_id);
    gfree(c->n, c->gpu_id);
    gfree(c->pavg, c->gpu_id);
    gfree(c->tavg, c->gpu_id);
    gfree(c->psavg, c->gpu_id);
    gfree(c->ns, c->gpu_id);
    gfree(c->linecenter, c->gpu_id);
    gfree(c->snn, c->gpu_id);
    gfree(c->gamma, c->gpu_id);
    gfree(c->alpha, c->gpu_id);
    if (c->gpu_id != HOST_ONLY)
    {
        gfree(c->p, c->gpu_id);
        gfree(c->t, c->gpu_id);
        gfree(c->tau, c->gpu_id);
    }
    if (c->use_h2o_ctm && is_molecule_active(c->molecule_bit_field, H2O))
    {
        catch(free_water_vapor_continuum_coefs(&(c->h2o_cc)));
    }
    if (c->use_o3_ctm && is_molecule_active(c->molecule_bit_field, O3))
    {
        catch(free_ozone_continuum_coefs(&(c->o3_cc)));
    }
    free(c);
    *context = NULL;
    return RS_SUCCESS;
}


/** @brief Add a molecule to a context
    @return RS_SUCCESS or an error code.*/
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
    if (is_molecule_active(context->molecule_bit_field, molecule_id))
    {
        char *mesg = "molecule %d has already been added.";
        raise(RS_VALUE_ERR, mesg, molecule_id);
    }
    int index = context->num_molecules;
    (context->num_molecules)++;
    in_range(context->num_molecules, 1, NUM_MOLS);
    catch(activate_molecule(&(context->molecule_bit_field), molecule_id));
    double w0;
    if (min_line_center != NULL)
    {
        in_range(*min_line_center, MIN_WAVENUMBER, MAX_WAVENUMBER);
        w0 = *min_line_center;
    }
    else
    {
        w0 = context->w0;
    }
    double wn;
    if (max_line_center != NULL)
    {
        in_range(*max_line_center, MIN_WAVENUMBER, MAX_WAVENUMBER);
        wn = *max_line_center;
    }
    else
    {
        wn = context->wn;
    }
    min_check(wn, w0);
    catch(molecule(&(context->mols[index]), molecule_id, context->hitran_path,
                   w0, wn, context->num_layers, context->gpu_id));
    char *mesg = "Using %s (%zu lines in range %e - %e [1/cm]).";
    log_mesg(mesg, context->mols[index].name,
             context->mols[index].line_params.num_lines, w0, wn);

    if (molecule_id == H2O && context->use_h2o_ctm)
    {
        /*Read in the water vapor continuum coefficients.*/
        mesg ="Using the %s continuum.";
        log_mesg(mesg, context->mols[index].name);
        catch(get_water_vapor_continuum_coefs(&(context->h2o_cc), context->h2o_ctm_dir,
                                              context->num_wpoints, context->w0,
                                              context->wres, context->gpu_id));
    }

    if (molecule_id == O3 && context->use_o3_ctm)
    {
        /*Read in the ozone continuum coefficients.*/
        mesg = "Using the %s continuum.";
        log_mesg(mesg, context->mols[index].name);
        catch(get_ozone_continuum_coefs(&(context->o3_cc), context->o3_ctm_dir,
                                        context->num_wpoints, context->w0,
                                        context->wres, context->gpu_id));
    }
    return RS_SUCCESS;
}


/*Update a molecule's ppmv.*/
EXTERN int grt_set_molecule_ppmv(GrtContext_t *context,
                                 int const molecule_id,
                                 fp_t const * const ppmv)
{
    not_null(context);
    not_null(ppmv);
    if (!is_molecule_active(context->molecule_bit_field, molecule_id))
    {
        char *mesg = "molecule %d is not being used.";
        log_warn(mesg, molecule_id);
        return RS_SUCCESS;
    }
    int index;
    catch(molecule_hash(molecule_id, &index));
    fp_t a[context->num_levels];
    int i;
    for (i=0; i<context->num_levels; ++i)
    {
        a[i] = ppmv[i]*1.e-6;
    }
    int offset = index*context->num_levels;
    gmemcpy(&(context->x[offset]), a, context->num_levels, context->gpu_id, FROM_HOST);
    return RS_SUCCESS;
}


/** @brief Add a CFC to a context
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_add_cfc(GrtContext_t *context, int const cfc_id,
                       char const * const filepath)
{
    not_null(context);
    if (is_cfc_active(context->cfc_bit_field, cfc_id))
    {
        char *mesg = "cfc %d has already been added.";
        raise(RS_VALUE_ERR, mesg, cfc_id);
    }
    int index = context->num_cfcs;
    (context->num_cfcs)++;
    in_range(context->num_cfcs, 1, NUM_CFCS);
    catch(activate_cfc(&(context->cfc_bit_field), cfc_id));

    /*Read in the CFC cross section values.*/
    catch(get_cfc_cross_sections(&(context->cfcs[index]), cfc_id, filepath,
                                 context->num_wpoints, context->w0, context->wres,
                                 context->gpu_id));
    char *mesg = "Using CFC %s.";
    log_mesg(mesg, context->cfcs[index].name);
    return RS_SUCCESS;
}


/*Update a CFC's ppmv.*/
EXTERN int grt_set_cfc_ppmv(GrtContext_t *context, int const cfc_id,
                            fp_t const * const ppmv)
{
    not_null(context);
    not_null(ppmv);
    if (!is_cfc_active(context->cfc_bit_field, cfc_id))
    {
        char *mesg = "CFC %d is not being used.";
        log_warn(mesg, cfc_id);
        return RS_SUCCESS;
    }
    fp_t a[context->num_levels];
    int i;
    for (i=0; i<context->num_levels; ++i)
    {
        a[i] = ppmv[i]*1.e-6;
    }
    int offset = cfc_id*context->num_levels;
    gmemcpy(&(context->x_cfc[offset]), a, context->num_levels, context->gpu_id, FROM_HOST);
    return RS_SUCCESS;
}


/*Calcluate the total optical depth in each layer at each spectral grid
  point.*/
EXTERN int grt_calculate_optical_depth(GrtContext_t *context, fp_t *pressure,
                                       fp_t *temperature, fp_t *optical_depth)
{
    not_null(context);
    not_null(pressure);
    not_null(temperature);
    not_null(optical_depth);
    fp_t const mbtoatm = 0.000986923f;
    fp_t p[MAX_NUM_LEVELS];
    int i;
    for (i=0; i<context->num_levels; ++i)
    {
        p[i] = pressure[i]*mbtoatm;
    }
    if (context->gpu_id == HOST_ONLY)
    {
        char *mesg = "Calculating optical depths for %d molecules in %d"
                     " atmospheric layers on the host CPU.";
        log_mesg(mesg, context->num_molecules, context->num_levels-1);
    }
    else
    {
        char *mesg = "Calculating optical depths for %d molecules in %d"
                     " atmospheric layers on GPU %d.";
        log_mesg(mesg, context->num_molecules, context->num_levels-1, context->gpu_id);
    }
    catch(launch(context, p, temperature, optical_depth));
    return RS_SUCCESS;
}


/*Get the number of molecules that have been added to the context.*/
EXTERN int grt_get_num_molecules(GrtContext_t const * const context,
                                 int * const n)
{
    not_null(context);
    not_null(n);
    *n = context->num_molecules;
    return RS_SUCCESS;
}


/*Get the number of spectral grid points for the input context.*/
EXTERN int grt_get_spectral_grid_size(GrtContext_t const * const context,
                                      uint64_t * const n)
{
    not_null(context);
    not_null(n);
    *n = context->num_wpoints;
    return RS_SUCCESS;
}


/*Return a message for an input return code.*/
EXTERN int grt_errstr(int const code, char * const buf, int const buf_size)
{
    not_null(buf);
    min_check(buf_size, 1);
    switch (code)
    {
        case RS_SUCCESS:
            break;
        case RS_INVALID_ERR:
            snprintf(buf, buf_size, "GRT: detected a floating point invalid.");
            break;
        case RS_DIVBYZERO_ERR:
            snprintf(buf, buf_size, "GRT: detected a floating point divide-by-zero.");
            break;
        case RS_OVERFLOW_ERR:
            snprintf(buf, buf_size, "GRT: detected a floating point overflow.");
            break;
        case RS_UNDERFLOW_ERR:
            snprintf(buf, buf_size, "GRT: detected a floating point underflow.");
            break;
        case RS_SENTINEL_ERR:
            snprintf(buf, buf_size, "GRT: entered unexpected code branch.");
            break;
        case RS_NULL_ERR:
            snprintf(buf, buf_size, "GRT: attempt to dereference a null pointer.");
            break;
        case RS_NON_NULL_ERR:
            snprintf(buf, buf_size,
                     "GRT: expected a null pointer, but pointer already has a value"
                         " assigned to it.");
            break;
        case RS_RANGE_ERR:
            snprintf(buf, buf_size, "GRT: detected a value out of its expected range.");
            break;
        case RS_VALUE_ERR:
            snprintf(buf, buf_size, "GRT: detected a bad value.");
            break;
        case RS_COMPILER_ERR:
            snprintf(buf, buf_size, "GRT: build was done with an incorrect compiler.");
            break;
        case RS_IO_ERR:
            snprintf(buf, buf_size, "GRT: error while performing I/O.");
            break;
        case RS_GPU_ERR:
            snprintf(buf, buf_size, "GRT: error while running on GPU.");
            break;
        default:
            snprintf(buf, buf_size, "Unknown code %d.", code);
    }
    return RS_SUCCESS;
}
