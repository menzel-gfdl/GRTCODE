#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "floating_point_type.h"
#include "parse_csv.h"
#include "solar_flux.h"
#include "utils.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif


/*Read in the incoming solar flux values.*/
int get_solar_flux(SolarFlux_t *sf,
                   unsigned int const nws,
                   int const w0,
                   double const res)
{
    not_null(sf);

    /*Set the file name.*/
    char *filepath = "INPUT/solar_flux.csv";
    int num_vals = 1;

    /*Read in the data.*/
    log_mesg("Reading in solar flux values from file %s.",
             filepath);
    int num_lines;
    int num_cols;
    char **buf;
    check(parse_csv(filepath,
                    &num_lines,
                    &num_cols,
                    1,
                    &buf));
    if ((num_vals + 1) != num_cols)
    {
        fatal("The number of columns (%d) in file %s does not match"
                  " the expected number (%d).",
              num_cols,
              filepath,
              num_vals + 1);
    }

    /*Convert the data from strings to floating point.*/
    fp_t *fbuf = NULL;
    int data_size = num_lines*num_cols;
    check(malloc_ptr((void **)(&fbuf),
                     sizeof(*fbuf)*data_size));
    int j;
    for (j=0;j<data_size;++j)
    {
        double d;
        check(to_double(buf[j],
                        &d));
        check(to_fp_t(d,
                      &(fbuf[j])));
        free(buf[j]);
    }
    free(buf);

    /*Allocate space for the solar fluxes.*/
    fp_t *c = NULL;
    int num_bytes = sizeof(*c)*nws;
    check(malloc_ptr((void **)(&c),
                     num_bytes));
    memset(c,
           0,
           num_bytes);

    /*Interpolate to wavenumber grid.*/
    fp_t *x = &(fbuf[0]);
    fp_t *y = &(fbuf[num_lines]);
    for (j=0;j<nws;++j)
    {
        check(linear_interpolation(x,
                                   y,
                                   num_lines,
                                   (fp_t)(w0 + j*res),
                                   &(c[j])));
    }
    free(fbuf);

    /*Integrate fluxes over wavenumber grid.*/
    check(reimann_sum(c,
                      nws,
                      res,
                      &(sf->total_sw_flux)));
    sf->incident_sw_flux = c;
    sf->nws = nws;
    return SUCCESS;
}


int free_solar_flux(SolarFlux_t *sf)
{
    not_null(sf);
    free(sf->incident_sw_flux);
    sf->incident_sw_flux = NULL;
    return SUCCESS;
}


int put_solar_flux_on_device(SolarFlux_t const * const in,
                             SolarFlux_t * const out)
{
    not_null(in);
    not_null(out);
    using_gpu();
    int num_bytes = sizeof(*(in->incident_sw_flux))*(in->nws);
    HANDLE_ERROR(cudaMalloc(&(out->incident_sw_flux),
                            num_bytes));
    HANDLE_ERROR(cudaMemcpy(out->incident_sw_flux,
                            in->incident_sw_flux,
                            num_bytes,
                            cudaMemcpyHostToDevice));
    out->total_sw_flux = in->total_sw_flux;
    return SUCCESS;
}


int remove_solar_flux_from_device(SolarFlux_t * const in)
{
    not_null(in);
    using_gpu();
    HANDLE_ERROR(cudaFree(in->incident_sw_flux));
    return SUCCESS;
}
