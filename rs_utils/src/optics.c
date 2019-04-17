#include <stdint.h>
#include "config.h"
#include "debug.h"
#include "floating_point_type.h"
#include "optics.h"
#include "spectral_grid.h"


/*Reserve memory for the optics.*/
int create_optics(Optics_t * const optics, int const num_layers, 
                  SpectralGrid_t const * const grid)
{
    not_null(optics);
    not_null(grid);
    in_range(num_layers, MIN_NUM_LAYERS, MAX_NUM_LAYERS);
    optics->num_layers = num_layers;
    char *mesg = "Initializing shortwave context:\nAtmospheric column properties:\n\t"
                     "number of levels: %d\n\tnumber of layers: %d";
    log_info(mesg, num_layers+1, num_layers);
    catch(create_spectral_grid(&(optics->grid), grid->w0, grid->wn, grid->dw,
                               &(grid->gpu_id)));
    uint64_t n = num_layers*(optics->grid.n);
    gmalloc(optics->g, n, optics->grid.gpu_id);
    gmemset(optics->g, 0, n, optics->grid.gpu_id);
    gmalloc(optics->omega, n, optics->grid.gpu_id);
    gmemset(optics->omega, 0, n, optics->grid.gpu_id);
    gmalloc(optics->tau, n, optics->grid.gpu_id);
    gmemset(optics->tau, 0, n, optics->grid.gpu_id);
    return RS_SUCCESS;
}


/*Free memory for the optics.*/
int destroy_optics(Optics_t * const optics)
{
    not_null(optics);
    gfree(optics->g, optics->grid.gpu_id);
    gfree(optics->omega, optics->grid.gpu_id);
    gfree(optics->tau, optics->grid.gpu_id);
    catch(destroy_spectral_grid(&(optics->grid)));
    return RS_SUCCESS;
}


/*Determine if two optics objects are compatible.*/
int optics_compatible(Optics_t const * const one, Optics_t const * const two,
                      int * const result)
{
    not_null(one);
    not_null(two);
    not_null(result);
    int same_grids;
    catch(compare_spectral_grids(&(one->grid), &(two->grid), &same_grids));
    if ((one->num_layers == two->num_layers) && (same_grids == 1))
    {
        *result = 1;
    }
    else
    {
        *result = 0;
    }
    return RS_SUCCESS;
}


/** @brief Add optics objects together.
    @return RS_SUCCESS or an error code.*/
static int add_optics_objects(uint64_t const n, /**< Size of arrays.*/
                              int const num_optics,
                              fp_t const * const g_in,
                              fp_t const * const omega_in,
                              fp_t const * const tau_in,
                              fp_t * const g_out,
                              fp_t * const omega_out,
                              fp_t * const tau_out
                             )
{
    uint64_t i;
#pragma omp parallel for default(none) private(i)
    for (i=0; i<n; ++i)
    {
        int j;
        for (j=0; j<num_optics; ++j)
        {
            uint64_t offset = j*n + i;
            g_out[i] += g_in[offset]*omega_in[offset]*tau_in[offset];
            omega_out[i] += omega_in[offset]*tau_in[offset];
            tau_out[i] += tau_in[offset];
        }
        g_out[i] /= omega_out[i];
        omega_out[i] /= tau_out[i];
    }
    return RS_SUCCESS;
}


#ifdef __NVCC__
/** @brief Add optics objects together.*/
__global__ static void add_optics_objects(uint64_t const n, /**< Size of arrays.*/
                                          int const num_optics,
                                          fp_t const * const g_in,
                                          fp_t const * const omega_in,
                                          fp_t const * const tau_in,
                                          fp_t * const g_out,
                                          fp_t * const omega_out,
                                          fp_t * const tau_out
                                         )
{
    uint64_t i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n)
    {
        int j;
        for (j=0; j<num_optics; ++j)
        {
            uint64_t offset = j*n + i;
            g_out[i] += g_in[offset]*omega_in[offset]*tau_in[offset];
            omega_out[i] += omega_in[offset]*tau_in[offset];
            tau_out[i] += tau_in[offset];
        }
        g_out[i] /= omega_out[i];
        omega_out[i] /= tau_out[i];
    }
    return;
}
#endif


/*Add optical properties together.*/
int add_optics(Optics_t const * const * const optics, int const num_optics,
               Optics_t * const result)
{
    not_null(optics);
    not_null(result);
    Optics_t const *o;
    int j;
    for (j=0; j<num_optics; ++j)
    {
        not_null(optics[j]);
        if (j == 0)
        {
            o = optics[j];
        }
        else
        {
            int ok;
            catch(optics_compatible(optics[j], o, &ok));
            if (ok == 0)
            {
                char *mesg = "input optics objects (%p, %p) are incompatible.";
                raise(RS_VALUE_ERR, mesg, o, optics[j]);
            }
        }
    }

    catch(create_optics(result, o->num_layers, &(o->grid)));
    fp_t *g = NULL;
    fp_t *omega = NULL;
    fp_t *tau = NULL;
    uint64_t n = o->num_layers*o->grid.n;
    gmalloc(g, n*num_optics, o->grid.gpu_id);
    gmalloc(omega, n*num_optics, o->grid.gpu_id);
    gmalloc(tau, n*num_optics, o->grid.gpu_id);
    for (j=0; j<num_optics; ++j)
    {
        gmemcpy(&(g[j*n]), optics[j]->g, n, o->grid.gpu_id, FROM_HOST);
        gmemcpy(&(omega[j*n]), optics[j]->omega, n, o->grid.gpu_id, FROM_HOST);
        gmemcpy(&(tau[j*n]), optics[j]->tau, n, o->grid.gpu_id, FROM_HOST);
    }
    glaunch(add_optics_objects, n, o->grid.gpu_id, n, num_optics, g, omega, tau, result->g,
            result->omega, result->tau);
    gfree(g, o->grid.gpu_id);
    gfree(omega, o->grid.gpu_id);
    gfree(tau, o->grid.gpu_id);
    return RS_SUCCESS;
}
