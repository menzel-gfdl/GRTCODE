#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "molecular_lines.h"


/*Utility macro to catch errors from library functions.*/
#define check_rc(rc) { \
    if (rc != 0) { \
        char buf[256]; \
        int e_ = grt_errstr(rc,buf,256); \
        fprintf(stderr,"%s\n",buf); \
        fprintf(stderr,"%s: %d Error\n",__FILE__,__LINE__); \
        return EXIT_FAILURE; \
    }}


/*Utility macro for switching precision.  If you want to run in double
  precision, you must include the -DDOUBLE_PRECISION when building the
  library.*/
#ifdef DOUBLE_PRECISION
#define FP_t double
#else
#define FP_t float
#endif


int main(void)
{
    GrtContext_t *context; /*Declare a library context pointer.*/
    int num_levels = 25; /*Number of atmospheric levels.*/
    double w0 = 2.; /*Lower bound [1/cm] for the spectral grid.*/
    double wn = 100.; /*Upper bound [1/cm] for the spectral grid.*/
    double wres = 0.1; /*Resoultion [1/cm] for the spectral grid.*/
    char *h2o_ctm_dir = "water_vapor_continuum"; /*Path of the directory
                                                   that contains the
                                                   water vapor continuum
                                                   input files.*/
    char *o3_ctm_dir = "ozone_continuum"; /*Path of the directory
                                            that contains the
                                            ozone continuum input file.*/

    /*Initalize the library context pointer.*/
    check_rc(grt_context_init(&context,
                              num_levels,
                              w0,
                              wn,
                              wres,
                              NULL,
                              NULL,
                              NULL,
                              h2o_ctm_dir,
                              o3_ctm_dir));

    char hitran_path[128];

    /*Set path to the water vapor HITRAN database input file.*/
    snprintf(hitran_path,128,"HITRAN_files/water_vapor.hitran12.par");

    /*Only water vapor lines with line centers in the range 1 - 1000 [1/cm]
      will be included in the calculation.*/
    double min_line_center_wavenumber = 1.;
    double max_line_center_wavenumber = 1000.;

    /*Add water vapor to the library context.*/
    int h2o; /*Water vapor molecule id.  Set by the library.*/
    check_rc(grt_add_molecule(context, 
                              hitran_path, 
                              &h2o, 
                              &min_line_center_wavenumber,
                              &max_line_center_wavenumber));

    /*Set path to the ozone HITRAN database input file.*/
    snprintf(hitran_path,128,"HITRAN_files/ozone.hitran12.par");

    /*Add ozone to the library context.*/
    int o3;
    check_rc(grt_add_molecule(context, 
                              hitran_path, 
                              &o3, 
                              NULL,
                              NULL));

    uint64_t num_wpoints; /*Size of the spectral grid.*/
    check_rc(grt_get_spectral_grid_size(context,
                                        &num_wpoints));
    int num_layers = num_levels - 1; /*Number of atmospheric layers.  This
                                       value is always equal to the number of
                                       atmospheric levels - 1.*/

    /*Allocate necessary arrays.*/
    FP_t *pressure = (FP_t *)malloc(sizeof(*pressure)*num_levels);
    FP_t *temperature = (FP_t *)malloc(sizeof(*temperature)*num_levels);
    FP_t *ppmv = (FP_t *)malloc(sizeof(*ppmv)*num_levels);
    FP_t *optical_depth = (FP_t *)malloc(sizeof(*optical_depth)*
                                         num_wpoints*num_layers);

    /*Loop over some columns.*/
    int num_columns = 4;
    int i;
    for (i=0;i<num_columns;++i)
    {
        /*Make up some data for the column.*/
        int j;
        for (j=0;j<num_levels;++j)
        {
            pressure[j] = 0.1 + 150.*j;
            temperature[j] = 230. + 2.3*j;
            ppmv[j] = 300. + 0.2*j;
        }

        /*Set the water vapor abundance for the library context.*/
        check_rc(grt_set_molecule_ppmv(context, 
                                       h2o, 
                                       ppmv));

        /*Make up some more data for the column.*/
        for (j=0;j<num_levels;++j)
        {
            ppmv[j] = 325. - 3.3*j;
        }

        /*Set the water vapor abundance for the library context.*/
        check_rc(grt_set_molecule_ppmv(context, 
                                       o3, 
                                       ppmv));

        /*Calculate the optical depths.*/
        check_rc(grt_calculate_optical_depth(context, 
                                             pressure, 
                                             temperature, 
                                             optical_depth));
    }

    /*Clean up.*/
    free(pressure);
    free(temperature);
    free(ppmv);
    free(optical_depth);

    /*Free memory allocated by the library context.*/
    check_rc(grt_context_free(&context));
}
