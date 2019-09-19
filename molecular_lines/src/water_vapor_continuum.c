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
#include "parse_csv.h"
#include "utils.h"
#include "water_vapor_continuum.h"


/*Read in the water vapor continuum coefficients.*/
int get_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc,
                                    char const * const h2o_ctm_dir,
                                    uint64_t const num_wpoints,
                                    double const w0,
                                    double const res,
                                    int const gpu_id)
{
    not_null(cc);
    not_null(h2o_ctm_dir);

    /*Set file names.*/
    char *filepath[NUM_COEFS];
    int num_vals[NUM_COEFS];
    size_t s = strlen(h2o_ctm_dir) + 64;
    int i;
    for (i=0; i<NUM_COEFS; ++i)
    {
        gmalloc(filepath[i], s, HOST_ONLY);
        switch (i)
        {
            case MTCKD25_F296:
                snprintf(filepath[i], s, "%s/296MTCKD25_F.csv", h2o_ctm_dir);
                num_vals[i] = 1;
                break;
            case MTCKD25_S296:
                snprintf(filepath[i], s, "%s/296MTCKD25_S.csv", h2o_ctm_dir);
                num_vals[i] = 1;
                break;
            case CKDF:
                snprintf(filepath[i], s, "%s/CKDF.csv", h2o_ctm_dir);
                num_vals[i] = 3;
                break;
            case CKDS:
                snprintf(filepath[i], s, "%s/CKDS.csv", h2o_ctm_dir);
                num_vals[i] = 3;
                break;
            default:
                sentinel();
        }
    }

    /*Allocate memory for each of the coefficient pointer.*/
    gmalloc(cc->coefs,NUM_COEFS,HOST_ONLY);
    for (i=0; i<NUM_COEFS; ++i)
    {
        /*Read in the data.*/
        char const *mesg = "Reading in water vapor continuum coefficients from file %s.";
        log_info(mesg, filepath[i]);
        int num_lines;
        int num_cols;
        char **buf;
        catch(parse_csv(filepath[i], &num_lines, &num_cols, 1, &buf));
        if ((num_vals[i]+1) != num_cols)
        {
            mesg = "The number of columns (%d) in file %s does not match"
                   " the expected number (%d).";
            raise(RS_VALUE_ERR, mesg, num_cols, filepath[i], num_vals[i]+1);
        }

        /*Convert the data from strings to floating point.*/
        fp_t *fbuf;
        int data_size = num_lines*num_cols;
        gmalloc(fbuf, data_size, HOST_ONLY);
        int j;
        for (j=0; j<data_size; ++j)
        {
            double d;
            catch(to_double(buf[j], &d));
            catch(to_fp_t(d, &(fbuf[j])));
            gfree(buf[j], HOST_ONLY);
        }
        gfree(buf, HOST_ONLY);

        /*Allocate space for the coefficient values at each wavenumber.*/
        fp_t *c;
        gmalloc(c, num_wpoints, HOST_ONLY);
        gmemset(c, 0, num_wpoints, HOST_ONLY);

        /*Interpolate to wavenumber grid.*/
        fp_t *x = &(fbuf[0]);
        fp_t *y = &(fbuf[num_lines]);
        for (j=0; (unsigned int)j<num_wpoints; ++j)
        {
            catch(linear_interpolation(x, y, num_lines, (fp_t)(w0 + j*res), &(c[j])));
        }
        gfree(fbuf, HOST_ONLY);
        if (gpu_id == HOST_ONLY)
        {
            cc->coefs[i] = c;
        }
        else
        {
            gmalloc(cc->coefs[i], num_wpoints, gpu_id);
            gmemcpy(cc->coefs[i], c, num_wpoints, gpu_id, FROM_HOST);
            gfree(c, HOST_ONLY);
        }
    }
    cc->num_wpoints = num_wpoints;
    cc->gpu_id = gpu_id;
    for (i=0; i<NUM_COEFS; ++i)
    {
        gfree(filepath[i], HOST_ONLY);
    }
    return RS_SUCCESS;
}


int free_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc)
{
    not_null(cc);
    int i;
    for (i=0; i<NUM_COEFS; ++i)
    {
        gfree(cc->coefs[i], cc->gpu_id);
    }
    gfree(cc->coefs, HOST_ONLY);
    return RS_SUCCESS;
}
