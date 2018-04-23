#include <stdlib.h>
#include "continuum_helpers.h"
#include "debug.h"
#include "floating_point_type.h"
#include "solar_flux.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif


/*Read in the incoming solar flux values.*/
int get_solar_flux(char const * const filepath,
                   SolarFlux_t *sf,
                   unsigned int const nws,
                   int const w0,
                   double const res,
                   int put_on_device)
{
    not_null(sf);
    enum solar_flux_csv_coefs
    {
        FLUX = 0,
        CSV_NUM_COEFS
    };
    fp_t *coefs[CSV_NUM_COEFS];
    log_mesg("Reading in solar flux values from file %s.",
             filepath);
    fp_t *data = NULL;
    int data_size;
    check(get_coefs(filepath,
                    coefs,
                    CSV_NUM_COEFS,
                    nws,
                    w0,
                    res,
                    &data,
                    &data_size));
    sf->incident_sw_flux = coefs[FLUX];

    int s = data_size/(CSV_NUM_COEFS + 1);
    double total = 0.;
    int i;
    for (i=0;i<(s-1);++i)
    {
        double dx = data[i+1] - data[i];
        total += 0.5*dx*(data[(FLUX+1)*s + i] + data[(FLUX+1)*s + i + 1]);
    }
    sf->total_sw_flux = (fp_t)total;
    if (put_on_device)
    {
        using_gpu();
#ifdef __NVCC__
        fatal("implement this on the device.");
#endif
    }
    return SUCCESS;
}


int free_solar_flux(SolarFlux_t *sf,
                    int const on_device)
{
    not_null(sf);
    if (on_device)
    {
        using_gpu();
#ifdef __NVCC__
        fatal("implement this on the device.");
#endif
    }
    else
    {
        free(sf->incident_sw_flux);
    }
    return SUCCESS;
}
