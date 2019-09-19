#include <stdint.h>
#include <stdio.h>
#include "collision_induced_absorption.h"
#include "debug.h"
#include "floating_point_type.h"
#include "parse_csv.h"
#include "utils.h"


/*Read in the collision-induced absorption cross sections.*/
int get_collision_induced_cross_sections(CollisionInducedAbsorption_t * const cia,
                                         int const id[2], char const * const filepath,
                                         uint64_t const num_wpoints, double const w0,
                                         double const res, int const gpu_id)
{
    not_null(cia);
    not_null(filepath);

    /*Set CIA metadata.*/
    int j;
    for (j=0; j<2; ++j)
    {
        int offset = j*CIA_NAME_LEN;
        switch(id[j])
        {
            case CIA_N2:
                snprintf(&(cia->name_buf[offset]), CIA_NAME_LEN, "N2");
                break;
            case CIA_O2:
                snprintf(&(cia->name_buf[offset]), CIA_NAME_LEN, "O2");
                break;
            default:
                {char const *mesg = "unrecognized CIA id %d.";
                raise(RS_VALUE_ERR, mesg, id[j]);}
        }
        cia->id[j] = id[j];
        cia->name[j] = &(cia->name_buf[offset]);
    }

    /*Read in the data.*/
    int num_vals = 1;
    char const *mesg = "Reading in collision-induced absorption cross sections from file %s.";
    log_info(mesg, filepath);
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
    fp_t *fbuf;
    int data_size = num_lines*num_cols;
    gmalloc(fbuf, data_size, HOST_ONLY);
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
        cia->cross_section = c;
    }
    else
    {
        gmalloc(cia->cross_section, num_wpoints, gpu_id);
        gmemcpy(cia->cross_section, c, num_wpoints, gpu_id, FROM_HOST);
        gfree(c, HOST_ONLY);
    }
    cia->num_wpoints = num_wpoints;
    cia->gpu_id = gpu_id;
    return RS_SUCCESS;
}


/*Free memory reserved for collision-induced absorption cross sections.*/
int free_collision_induced_cross_sections(CollisionInducedAbsorption_t * const cia)
{
    not_null(cia);
    gfree(cia->cross_section, cia->gpu_id);
    return RS_SUCCESS;
}
