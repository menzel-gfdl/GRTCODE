#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "cfcs.h"
#include "debug.h"
#include "floating_point_type.h"
#include "parse_csv.h"
#include "utils.h"


/*Read in the cfc cross sections.*/
int get_cfc_cross_sections(CfcCrossSection_t *xsc, int const id, char const * const filepath,
                           uint64_t const num_wpoints, double const w0,
                           double const res, int const gpu_id)
{
    not_null(xsc);
    not_null(filepath);

    /*Set CFC metadata.*/
    xsc->id = id;
    switch(id)
    {
        case F11:
            snprintf(xsc->name, CFC_NAME_LEN, "F11");
            break;
        case F12:
            snprintf(xsc->name, CFC_NAME_LEN, "F12");
            break;
        default:
            {char *mesg = "unrecognized CFC id %d.";
            raise(RS_VALUE_ERR, mesg, id);}
    }

    /*Read in the data.*/
    int num_lines;
    int num_cols;
    char **buf;
    catch(parse_csv(filepath, &num_lines, &num_cols, 1, &buf));
    int const ncols_req = 2;
    if (num_cols != ncols_req)
    {
        char *mesg = "The number of columns (%d) in file %s does not match"
                     " the expected number (%d).";
        raise(RS_VALUE_ERR, mesg, num_cols, filepath, ncols_req);
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

    /*Allocate space for the cross section values at each wavenumber.*/
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
        xsc->cross_section = c;
    }
    else
    {
        gmalloc(xsc->cross_section, num_wpoints, gpu_id);
        gmemcpy(xsc->cross_section, c, num_wpoints, gpu_id, FROM_HOST);
        gfree(c, HOST_ONLY);
    }
    xsc->num_wpoints = num_wpoints;
    xsc->gpu_id = gpu_id;
    return RS_SUCCESS;
}


/*Free memory.*/
int free_cfc_cross_sections(CfcCrossSection_t *xsc)
{
    not_null(xsc);
    gfree(xsc->cross_section, xsc->gpu_id);
    return RS_SUCCESS;
}


int activate_cfc(uint32_t * const cfc_bit_field, int const id)
{
    not_null(cfc_bit_field);
    in_range(id, 0, NUM_CFCS);
    uint32_t const one = 1;
    *cfc_bit_field = (*cfc_bit_field) | (one << id);
    return RS_SUCCESS;
}


int is_cfc_active(uint32_t const cfc_bit_field, int const id)
{
    in_range(id, 0, NUM_CFCS);
    uint32_t const one = 1;
    return cfc_bit_field & (one << id);
}
