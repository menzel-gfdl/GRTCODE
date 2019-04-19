#include <stdint.h>
#include "debug.h"
#include "extern.h"
#include "floating_point_type.h"
#include "parse_csv.h"
#include "solar_flux.h"
#include "utils.h"


/*Read in data for the solar flux.*/
EXTERN int create_solar_flux(SolarFlux_t * const solar_flux, SpectralGrid_t const * const grid,
                             char const * const filepath)
{
    not_null(solar_flux);
    not_null(grid);
    not_null(filepath);

    /*Read in the data.*/
    char *mesg = "Reading in solar flux values from file %s.";
    log_info(mesg, filepath);
    int num_vals = 1;
    int num_lines;
    int num_cols;
    char **buf;
    catch(parse_csv(filepath, &num_lines, &num_cols, 1, &buf));
    if ((num_vals + 1) != num_cols)
    {
        mesg = "The number of columns (%d) in file %s does not match"
               " the expected number (%d).";
        raise(RS_VALUE_ERR, mesg, num_cols, filepath, num_vals+1);
    }

    /*Convert the data from strings to floating point.*/
    fp_t *fbuf = NULL;
    uint64_t data_size = num_lines*num_cols;
    gmalloc(fbuf, data_size, HOST_ONLY);
    uint64_t j;
    for (j=0; j<data_size; ++j)
    {
        double d;
        catch(to_double(buf[j], &d));
        catch(to_fp_t(d, &(fbuf[j])));
        gfree(buf[j], HOST_ONLY);
    }
    gfree(buf, HOST_ONLY);

    /*Allocate space for the solar fluxes.*/
    solar_flux->grid = *grid;
    fp_t *c = NULL;
    gmalloc(c, grid->n, HOST_ONLY);
    gmemset(c, 0, grid->n, HOST_ONLY);

    /*Interpolate to spectral grid.*/
    fp_t *x = &(fbuf[0]);
    fp_t *y = &(fbuf[num_lines]);
    for (j=0; j<grid->n; ++j)
    {
        catch(linear_interpolation(x, y, num_lines, (fp_t)(grid->w0 + j*grid->dw),
                                   &(c[j])));
    }
    gfree(fbuf, HOST_ONLY);

    /*Integrate fluxes over spectral grid.*/
    fp_t total_flux;
    catch(reimann_sum(c, grid->n, grid->dw, &total_flux));

    /*Adjust the spectral.*/
    for (j=0; j<grid->n; ++j)
    {
        c[j] /= total_flux;
    }
    solar_flux->incident_flux = c;
    solar_flux->n = grid->n;
    return RS_SUCCESS;
}


/*Free memory for the solar flux.*/
int destroy_solar_flux(SolarFlux_t * const solar_flux)
{
    not_null(solar_flux);
    gfree(solar_flux->incident_flux, HOST_ONLY);
    return RS_SUCCESS;
}
