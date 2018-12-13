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
#include "floating_point_type.h"
#include "ozone_continuum.h"
#include "parse_csv.h"
#include "utils.h"


/*Read in the ozone continuum coefficients.*/
int get_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc,
                              char const * const o3_ctm_dir,
                              uint64_t const num_wpoints,
                              double const w0,
                              double const res,
                              int const gpu_id)
{
    not_null(cc);
    not_null(o3_ctm_dir);

    /*Set file name.*/
    char *filepath;
    size_t s = strlen(o3_ctm_dir) + 64;
    gmalloc(filepath,s,HOST_ONLY);
    snprintf(filepath,
             s,
             "%s/ozone_continuum.csv",
             o3_ctm_dir);
    int num_vals = 1;

    /*Read in the data.*/
    log_info("Reading in ozone continuum coefficients from file %s.",
             filepath);
    int num_lines;
    int num_cols;
    char **buf;
    throw(parse_csv(filepath,
                    &num_lines,
                    &num_cols,
                    1,
                    &buf));
    if ((num_vals + 1) != num_cols)
    {
        raise(VALUE_ERR,
              "The number of columns (%d) in file %s does not match"
                 " the expected number (%d).",
              num_cols,
              filepath,
              num_vals + 1);
    }

    /*Convert the data from strings to floating point.*/
    fp_t *fbuf;
    int data_size = num_lines*num_cols;
    gmalloc(fbuf,data_size,HOST_ONLY);
    int j;
    for (j=0;j<data_size;++j)
    {
        double d;
        throw(to_double(buf[j],
                        &d));
        throw(to_fp_t(d,
                      &(fbuf[j])));
        gfree(buf[j],HOST_ONLY);
    }
    gfree(buf,HOST_ONLY);

    /*Allocate space for the coefficient values at each wavenumber.*/
    fp_t *c;
    gmalloc(c,num_wpoints,HOST_ONLY);
    gmemset(c,0,num_wpoints,HOST_ONLY);

    /*Interpolate to wavenumber grid.*/
    fp_t *x = &(fbuf[0]);
    fp_t *y = &(fbuf[num_lines]);
    for (j=0;(unsigned int)j<num_wpoints;++j)
    {
        throw(linear_interpolation(x,
                                   y,
                                   num_lines,
                                   (fp_t)(w0 + j*res),
                                   &(c[j])));
    }
    gfree(fbuf,HOST_ONLY);
    if (gpu_id == HOST_ONLY)
    {
        cc->cross_section = c;
    }
    else
    {
        gmalloc(cc->cross_section,num_wpoints,gpu_id);
        gmemcpy(cc->cross_section,c,num_wpoints,gpu_id,FROM_HOST);
        gfree(c,HOST_ONLY);
    }
    cc->num_wpoints = num_wpoints;
    cc->gpu_id = gpu_id;
    gfree(filepath,HOST_ONLY);
    return SUCCESS;
}


int free_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc)
{
    not_null(cc);
    gfree(cc->cross_section,cc->gpu_id);
    return SUCCESS;
}
