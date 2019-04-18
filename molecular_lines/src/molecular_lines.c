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


/** @brief Reserve memory for molecular lines.
    @return RS_SUCCESS or an error code.*/
EXTERN int create_molecular_lines(MolecularLines_t * const ml, /**< Molecular lines object.*/
                                  int const num_levels, /**< Number of atmospheric levels.*/
                                  SpectralGrid_t const * const grid, /**< Spectral grid.*/
                                  Device_t const * const device, /**< Device.*/
                                  char const * const hitran_path, /**< Path to HITRAN database file.*/
                                  char const * const h2o_ctm_dir, /**< Path to water vapor continuum directory.*/
                                  char const * const o3_ctm_dir, /**< Path to ozone continuum directory.*/
                                  double const * const wcutoff, /**< Cutoff from line center [1/cm].*/
                                  int const * const optical_depth_method /**< Method to use to calculate optical depths.*/
                                 )
{
    not_null(ml);
    in_range(num_levels, MIN_NUM_LEVELS, MAX_NUM_LEVELS);
    ml->num_levels = num_levels;
    ml->num_layers = num_levels - 1;
    char *mesg = "Molecular lines properties:\n\tnumber of levels: %d\n\t"
                 "number of layers: %d";
    log_info(mesg, ml->num_levels, ml->num_layers);
    not_null(grid);
    ml->grid = *grid;
    not_null(device);
    ml->device = *device;

    /*Create the spectral bins.*/
    double bin_width = 1.;
    catch(create_spectral_bins(&(ml->bins), ml->num_layers, ml->grid.w0, ml->grid.n,
                               ml->grid.dw, bin_width, ml->device));
    mesg = "Spectral bin properties:\n\tnumber of bins: %zu\n\t"
           "bin width: %e\n\tspectral grid points per bin: %d\n\t"
           "interpolation: %d\n\tspectral gid points in last bin:"
           " %d\n\tinterpolation in last bin: %d";
    log_info(mesg, ml->bins.n, bin_width, ml->bins.ppb, ml->bins.do_interp,
             ml->bins.last_ppb, ml->bins.do_last_interp);

    /*Store the path to the hitran database file.*/
    not_null(hitran_path);
    snprintf(ml->hitran_path, DIR_PATH_LEN, "%s", hitran_path);
    mesg = "Using HITRAN database file %s.";
    log_info(mesg, ml->hitran_path);

    /*Set the molecular line cutoff.*/
    if (wcutoff != NULL)
    {
        in_range(*wcutoff, MIN_CUTOFF, MAX_CUTOFF);
        ml->wcutoff = *wcutoff;
    }
    else
    {
        ml->wcutoff = DEFAULT_CUTOFF;
    }
    mesg = "Using spectral line cut-off of %e [1/cm].";
    log_info(mesg, ml->wcutoff);

    /*Set the method that will be used to calculate the optical depths.*/
    if (optical_depth_method != NULL)
    {
        in_range(*optical_depth_method, wavenumber_sweep, line_sample);
        ml->optical_depth_method = *optical_depth_method;
    }
    else
    {
        ml->optical_depth_method = wavenumber_sweep;
    }

    /*Prepare to add molecules/cfcs.*/
    ml->num_molecules = 0;
    ml->molecule_bit_field = 0;
    ml->num_cfcs = 0;
    ml->cfc_bit_field = 0;

    /*Pepare water vapor continuum.*/
    ml->use_h2o_ctm = 0;
    if (h2o_ctm_dir != NULL)
    {
        if (strcmp(h2o_ctm_dir, "none") != 0)
        {
            ml->use_h2o_ctm = 1;
            catch(copy_str(ml->h2o_ctm_dir, h2o_ctm_dir, DIR_PATH_LEN));
        }
    }

    /*Prepare ozone continuum.*/
    ml->use_o3_ctm = 0;
    if (o3_ctm_dir != NULL)
    {
        if (strcmp(o3_ctm_dir, "none") != 0)
        {
            ml->use_o3_ctm = 1;
            catch(copy_str(ml->o3_ctm_dir, o3_ctm_dir, DIR_PATH_LEN));
        }
    }

    /*Reserve memory.*/
    gmalloc(ml->x, ml->num_levels*NUM_MOLS, ml->device);
    gmalloc(ml->x_cfc, ml->num_levels*NUM_CFCS, ml->device);
    gmalloc(ml->n, ml->num_layers, ml->device);
    gmalloc(ml->pavg, ml->num_layers, ml->device);
    gmalloc(ml->tavg, ml->num_layers, ml->device);
    gmalloc(ml->psavg, ml->num_layers, ml->device);
    gmalloc(ml->ns, ml->num_layers, ml->device);
    gmalloc(ml->linecenter, ml->num_layers*MAX_NUM_LINES, ml->device);
    gmalloc(ml->snn, ml->num_layers*MAX_NUM_LINES, ml->device);
    gmalloc(ml->gamma, ml->num_layers*MAX_NUM_LINES, ml->device);
    gmalloc(ml->alpha, ml->num_layers*MAX_NUM_LINES, ml->device);
    if (ml->device != HOST_ONLY)
    {
        gmalloc(ml->p, ml->num_levels, ml->device);
        gmalloc(ml->t, ml->num_levels, ml->device);
        gmalloc(ml->tau, ml->num_layers*ml->grid.n, ml->device);
    }

    /*Initialize TIPS.*/
    if (ml->device != HOST_ONLY)
    {
        catch(inittips_d());
    }
    return RS_SUCCESS;
}


/** @brief Free memory for the molecular lines.
    @return RS_SUCCESS or an error code.*/
EXTERN int destroy_molecular_lines(MolecularLines_t * const ml /**< Molecular lines object.*/
                                  )
{
    not_null(ml);
    int i;
    for (i=0; i<ml->num_molecules; ++i)
    {
        catch(free_molecule(&(ml->mols[i])));
    }
    for (i=0; i<ml->num_cfcs; ++i)
    {
        catch(free_cfc_cross_sections(&(ml->cfcs[i])));
    }
    catch(destroy_spectral_bins(&(ml->bins)));
    gfree(ml->x, ml->device);
    gfree(ml->x_cfc, ml->device);
    gfree(ml->n, ml->device);
    gfree(ml->pavg, ml->device);
    gfree(ml->tavg, ml->device);
    gfree(ml->psavg, ml->device);
    gfree(ml->ns, ml->device);
    gfree(ml->linecenter, ml->device);
    gfree(ml->snn, ml->device);
    gfree(ml->gamma, ml->device);
    gfree(ml->alpha, ml->device);
    if (ml->device != HOST_ONLY)
    {
        gfree(ml->p, ml->device);
        gfree(ml->t, ml->device);
        gfree(ml->tau, ml->device);
    }
    if (ml->use_h2o_ctm && is_molecule_active(ml->molecule_bit_field, H2O))
    {
        catch(free_water_vapor_continuum_coefs(&(ml->h2o_cc)));
    }
    if (ml->use_o3_ctm && is_molecule_active(ml->molecule_bit_field, O3))
    {
        catch(free_ozone_continuum_coefs(&(ml->o3_cc)));
    }
    return RS_SUCCESS;
}


/** @brief Add a molecule.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_add_molecule(MolecularLines_t * const ml, /**< Molecular lines object.*/
                            int const molecule_id, /**< Molecule id.*/
                            double const * const min_line_center, /**< Lower bound [1/cm]
                                                                       for spectral line
                                                                       centers.*/
                            double const * const max_line_center /**< Upper bound [1/cm]
                                                                      for spectral line
                                                                      centers.*/
                           )
{
    not_null(ml);
    if (is_molecule_active(ml->molecule_bit_field, molecule_id))
    {
        char *mesg = "molecule %d has already been added.";
        raise(RS_VALUE_ERR, mesg, molecule_id);
    }
    int index = ml->num_molecules;
    (ml->num_molecules)++;
    in_range(ml->num_molecules, 1, NUM_MOLS);
    catch(activate_molecule(&(ml->molecule_bit_field), molecule_id));
    double w0;
    if (min_line_center != NULL)
    {
        in_range(*min_line_center, MIN_WAVENUMBER, MAX_WAVENUMBER);
        w0 = *min_line_center;
    }
    else
    {
        w0 = ml->grid.w0;
    }
    double wn;
    if (max_line_center != NULL)
    {
        in_range(*max_line_center, MIN_WAVENUMBER, MAX_WAVENUMBER);
        wn = *max_line_center;
    }
    else
    {
        wn = ml->grid.wn;
    }
    min_check(wn, w0);
    catch(molecule(&(ml->mols[index]), molecule_id, ml->hitran_path,
                   w0, wn, ml->num_layers, ml->device));
    char *mesg = "Using %s (%zu lines in range %e - %e [1/cm]).";
    log_mesg(mesg, ml->mols[index].name, ml->mols[index].line_params.num_lines, w0, wn);

    if (molecule_id == H2O && ml->use_h2o_ctm)
    {
        /*Read in the water vapor continuum coefficients.*/
        mesg ="Using the %s continuum.";
        log_mesg(mesg, ml->mols[index].name);
        catch(get_water_vapor_continuum_coefs(&(ml->h2o_cc), ml->h2o_ctm_dir,
                                              ml->grid.n, ml->grid.w0,
                                              ml->grid.dw, ml->device));
    }

    if (molecule_id == O3 && ml->use_o3_ctm)
    {
        /*Read in the ozone continuum coefficients.*/
        mesg = "Using the %s continuum.";
        log_mesg(mesg, ml->mols[index].name);
        catch(get_ozone_continuum_coefs(&(ml->o3_cc), ml->o3_ctm_dir,
                                        ml->grid.n, ml->grid.w0,
                                        ml->grid.dw, ml->device));
    }
    return RS_SUCCESS;
}


/** @brief Update a molecule's ppmv.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_set_molecule_ppmv(MolecularLines_t * const ml, /**< Molecular lines object.*/
                                 int const molecule_id, /**< Molecule id.*/
                                 fp_t const * const ppmv /**< Abundance [ppmv] (level).*/
                                )
{
    not_null(ml);
    not_null(ppmv);
    if (!is_molecule_active(ml->molecule_bit_field, molecule_id))
    {
        char *mesg = "molecule %d is not being used.";
        log_warn(mesg, molecule_id);
        return RS_SUCCESS;
    }
    int index;
    catch(molecule_hash(molecule_id, &index));
    fp_t a[ml->num_levels];
    int i;
    for (i=0; i<ml->num_levels; ++i)
    {
        a[i] = ppmv[i]*1.e-6;
    }
    int offset = index*ml->num_levels;
    gmemcpy(&(ml->x[offset]), a, ml->num_levels, ml->device, FROM_HOST);
    return RS_SUCCESS;
}


/** @brief Add a CFC.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_add_cfc(MolecularLines_t * const ml, /**< Molecular lines object.*/
                       int const cfc_id, /**< CFC id.*/
                       char const * const filepath /**< Path to CFC cross section csv file.*/
                      )
{
    not_null(ml);
    if (is_cfc_active(ml->cfc_bit_field, cfc_id))
    {
        char *mesg = "cfc %d has already been added.";
        raise(RS_VALUE_ERR, mesg, cfc_id);
    }
    int index = ml->num_cfcs;
    (ml->num_cfcs)++;
    in_range(ml->num_cfcs, 1, NUM_CFCS);
    catch(activate_cfc(&(ml->cfc_bit_field), cfc_id));

    /*Read in the CFC cross section values.*/
    catch(get_cfc_cross_sections(&(ml->cfcs[index]), cfc_id, filepath,
                                 ml->grid.n, ml->grid.w0, ml->grid.dw,
                                 ml->device));
    char *mesg = "Using CFC %s.";
    log_mesg(mesg, ml->cfcs[index].name);
    return RS_SUCCESS;
}


/** @brief Update a CFC's ppmv.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_set_cfc_ppmv(MolecularLines_t * const ml, /**< Molecular lines object.*/
                            int const cfc_id, /**< CFC id.*/
                            fp_t const * const ppmv /**< Abundance [ppmv] (level).*/
                           )
{
    not_null(ml);
    not_null(ppmv);
    if (!is_cfc_active(ml->cfc_bit_field, cfc_id))
    {
        char *mesg = "CFC %d is not being used.";
        log_warn(mesg, cfc_id);
        return RS_SUCCESS;
    }
    fp_t a[ml->num_levels];
    int i;
    for (i=0; i<ml->num_levels; ++i)
    {
        a[i] = ppmv[i]*1.e-6;
    }
    int offset = cfc_id*ml->num_levels;
    gmemcpy(&(ml->x_cfc[offset]), a, ml->num_levels, ml->device, FROM_HOST);
    return RS_SUCCESS;
}


/** @brief Calcluate the total optical depth in each layer at each spectral grid point.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_calculate_optical_depth(MolecularLines_t * const ml, /**< Molecular lines object.*/
                                       fp_t * const pressure, /**< Pressure [mb] (level).*/
                                       fp_t * const temperature, /**< Temperature [K] (level).*/
                                       Optics_t * const optics /**< Optics object.*/
                                      )
{
    not_null(ml);
    not_null(pressure);
    not_null(temperature);
    not_null(optics);
    assert(ml->device, optics->device);
    assert(ml->num_layers, optics->num_layers);
    int same_grids;
    catch(compare_spectral_grids(&(ml->grid), &(optics->grid), &same_grids));
    assert(same_grids, 1);
    fp_t const mbtoatm = 0.000986923f;
    fp_t p[MAX_NUM_LEVELS];
    int i;
    for (i=0; i<ml->num_levels; ++i)
    {
        p[i] = pressure[i]*mbtoatm;
    }
    catch(launch(ml, p, temperature, optics->tau));
    return RS_SUCCESS;
}


/** @brief Get the number of molecules.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_get_num_molecules(MolecularLines_t const * const ml, /**< Molecular lines object.*/
                                 int * const n /**< Number of molecules.*/
                                )
{
    not_null(ml);
    not_null(n);
    *n = ml->num_molecules;
    return RS_SUCCESS;
}


/** @brief Get the number of spectral grid points.
    @return RS_SUCCESS or an error code.*/
EXTERN int grt_get_spectral_grid_size(MolecularLines_t const * const ml, /**< Molecular lines object.*/
                                      uint64_t * const n /**< Spectral grid size.*/
                                     )
{
    not_null(ml);
    not_null(n);
    *n = ml->grid.n;
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
