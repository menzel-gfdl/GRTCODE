#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#else
#error "You must build this code with OpenMP."
#endif
#include "molecular_lines.h"


/*Set the number GPUs you want to run on.*/
#define NUM_GPUS 1


/*Utility macro to catch errors from library functions.*/
#define check_rc(rc) { \
    if (rc != 0) { \
        char buf[256]; \
        int e_ = grt_errstr(rc,buf,256); \
        fprintf(stderr,"%s\n",buf); \
        fprintf(stderr,"%s: %d Error\n",__FILE__,__LINE__); \
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
    GrtContext_t *context[NUM_GPUS]; /*Declare library context pointers.*/
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
    int h2o[NUM_GPUS]; /*Water vapor molecule id.  Set by the library.*/
    int o3[NUM_GPUS]; /*Ozone molecule id.  Set by the library.*/
    int i;

#pragma omp parallel for num_threads(NUM_GPUS) \
                         default(none) \
                         shared(context,num_levels,w0,wn,wres, \
                                h2o_ctm_dir,o3_ctm_dir,h2o,o3,stderr) \
                         private(i) /*is watching you.*/
    for (i=0;i<NUM_GPUS;++i)
    {
        int gpu_id = omp_get_thread_num();
        GrtContext_t *c = context[gpu_id];

        /*Initalize library context pointers.*/
        check_rc(grt_context_init(&c,
                                  num_levels,
                                  w0,
                                  wn,
                                  wres,
                                  NULL,
                                  &gpu_id,
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
        check_rc(grt_add_molecule(c, 
                                  hitran_path, 
                                  &(h2o[gpu_id]), 
                                  &min_line_center_wavenumber,
                                  &max_line_center_wavenumber));

        /*Set path to the ozone HITRAN database input file.*/
        snprintf(hitran_path,128,"HITRAN_files/ozone.hitran12.par");

        /*Add ozone to the library context.*/
        check_rc(grt_add_molecule(c, 
                                  hitran_path, 
                                  &(o3[gpu_id]), 
                                  NULL,
                                  NULL));
    }

    uint64_t num_wpoints; /*Size of the spectral grid.  Here both context
                            pointers have the same spectral grid, so we'll
                            just use one of them to get the size.*/
    check_rc(grt_get_spectral_grid_size(context[0],
                                        &num_wpoints));
    int num_layers = num_levels - 1; /*Number of atmospheric layers.  This
                                       value is always equal to the number of
                                       atmospheric levels - 1.*/

    /*Allocate necessary arrays.*/
    int num_columns = 2*NUM_GPUS;
    FP_t *pressure = (FP_t *)malloc(sizeof(*pressure)*num_levels*num_columns);
    FP_t *temperature = (FP_t *)malloc(sizeof(*temperature)*num_levels*
                                       num_columns);
    FP_t *ppmv = (FP_t *)malloc(sizeof(*ppmv)*num_levels*num_columns);
    FP_t *optical_depth = (FP_t *)malloc(sizeof(*optical_depth)*
                                         num_wpoints*num_layers*num_columns);

    /*Loop over some columns.*/
#pragma omp parallel for num_threads(NUM_GPUS) \
                         default(none) \
                         shared(context,num_levels,num_columns,num_wpoints, \
                                pressure,temperature,ppmv,h2o,o3,stderr,\
                                optical_depth) \
                         private(i)
    for (i=0;i<num_columns;++i)
    {
        int gpu_id = omp_get_thread_num();
        GrtContext_t *c = context[gpu_id];
        int offset = i*num_columns;

        /*Make up some data for the column.*/
        int j;
        for (j=0;j<num_levels;++j)
        {
            int o = i*num_levels + j;
            pressure[o] = 0.1 + 150.*j;
            temperature[o] = 230. + 2.3*j;
            ppmv[o] = 300. + 0.2*j;
        }

        /*Set the water vapor abundance for the library context.*/
        check_rc(grt_set_molecule_ppmv(c, 
                                       h2o[gpu_id], 
                                       &(ppmv[offset])));

        /*Make up some more data for the column.*/
        for (j=0;j<num_levels;++j)
        {
            ppmv[i*num_levels+j] = 325. - 3.3*j;
        }

        /*Set the water vapor abundance for the library context.*/
        check_rc(grt_set_molecule_ppmv(c, 
                                       o3[gpu_id], 
                                       &(ppmv[offset])));

        /*Calculate the optical depths.*/
        check_rc(grt_calculate_optical_depth(c, 
                                             &(pressure[offset]), 
                                             &(temperature[offset]), 
                                             &(optical_depth[offset*num_wpoints])));
    }

    /*Clean up.*/
    free(pressure);
    free(temperature);
    free(ppmv);
    free(optical_depth);

    /*Free memory allocated by the library context.*/
#pragma omp parallel for num_threads(NUM_GPUS) \
                         default(none) \
                         shared(context,stderr) \
                         private(i)
    for (i=0;i<NUM_GPUS;++i)
    {
        int gpu_id = omp_get_thread_num();
        check_rc(grt_context_free(&(context[gpu_id])));
    }
}
