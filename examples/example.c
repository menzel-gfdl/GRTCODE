#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#else
#error "You must build this code with OpenMP."
#endif
#include "molecular_lines.h"


/*Utility macros to catch errors from library functions.*/
#define check_rc(rc) { \
    if (rc != 0) { \
        char buf[256]; \
        grt_errstr(rc,buf,256); \
        fprintf(stderr,"%s\n",buf); \
        fprintf(stderr,"%s: %d Error\n",__FILE__,__LINE__); \
        return rc; \
    } \
}


#define omp_check_rc(rc,e) { \
    if (rc != 0) { \
        char buf[256]; \
        grt_errstr(rc,buf,256); \
        fprintf(stderr,"%s\n",buf); \
        fprintf(stderr,"%s: %d Error\n",__FILE__,__LINE__); \
        e |= rc; \
    } \
}


/*Utility macro for switching precision.  If you want to run in double
  precision, you must include the -DDOUBLE_PRECISION when building the
  library.*/
#ifdef SINGLE_PRECISION
#define FP_t float
#else
#define FP_t double
#endif


int main(int argc,char **argv)
{
    /*Command line argument controls how many GPUs will be used.*/
    int num_contexts = 1;
    int host_only = 0;
    int const host_id = -1;
    if (argc == 2)
    {
        if (strcmp("--host",argv[1]) == 0)
        {
            host_only = 1;
        }
        else
        {
            num_contexts = atoi(argv[1]);
        }
    }
    else if (argc > 2)
    {
        fprintf(stderr,"Usage: %s [--host|num_gpus]\n",argv[0]);
        return EXIT_FAILURE;
    }

    GrtContext_t *context[num_contexts]; /*Declare library context pointers.*/
    int num_levels = 25; /*Number of atmospheric levels.*/
    double w0 = 1.; /*Lower bound [1/cm] for the spectral grid.*/
    double wn = 3000.; /*Upper bound [1/cm] for the spectral grid.*/
    double wres = 0.1; /*Resoultion [1/cm] for the spectral grid.*/
    char *h2o_ctm_dir = "water_vapor_continuum"; /*Path of the directory
                                                   that contains the
                                                   water vapor continuum
                                                   input files.*/
    char *o3_ctm_dir = "ozone_continuum"; /*Path of the directory
                                            that contains the
                                            ozone continuum input file.*/
    int i;

    /*Increase the verbosity of the library output.*/
    grt_set_verbosity(3);

    for (i=0;i<num_contexts;++i)
    {
        /*Determine the GPU id for the context.*/
        int const *g;
        if (argc == 1)
        {
            g = NULL;
        }
        else if (host_only)
        {
            g = &host_id;
        }
        else
        {
            g = &i;
        }

        /*Initalize library context pointers.*/
        GrtContext_t *c;
        check_rc(grt_context_init(&c,
                                  num_levels,
                                  w0,
                                  wn,
                                  wres,
                                  "HITRAN_files/hitran2012.par",
                                  h2o_ctm_dir,
                                  o3_ctm_dir,
                                  NULL,
                                  g,
                                  NULL));
        context[i] = c;

        /*Only water vapor lines with line centers in the range 1 - 1000 [1/cm]
          will be included in the calculation.*/
        double min_line_center_wavenumber = 1.;
        double max_line_center_wavenumber = 1000.;

        /*Add water vapor to the library context.*/
        check_rc(grt_add_molecule(c,
                                  H2O,
                                  &min_line_center_wavenumber,
                                  &max_line_center_wavenumber));

        /*Add ozone to the library context.*/
        check_rc(grt_add_molecule(c,
                                  O3,
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
    int num_columns = 2*num_contexts;
    FP_t *pressure = (FP_t *)malloc(sizeof(*pressure)*num_levels*num_columns);
    FP_t *temperature = (FP_t *)malloc(sizeof(*temperature)*num_levels*
                                       num_columns);
    FP_t *ppmv = (FP_t *)malloc(sizeof(*ppmv)*num_levels*num_columns);
    FP_t *optical_depth = (FP_t *)malloc(sizeof(*optical_depth)*
                                         num_wpoints*num_layers*num_columns);
    int rc[num_contexts];
    memset(rc,
           0,
           sizeof(rc)*num_contexts);

    /*Loop over some columns.*/
#pragma omp parallel for num_threads(num_contexts) \
                         default(none) \
                         shared(context,num_levels,num_columns,num_wpoints, \
                                pressure,temperature,ppmv,stderr, \
                                optical_depth,rc) \
                         private(i) /*is watching you, seeing your every move.*/
    for (i=0;i<num_columns;++i)
    {
        int g = omp_get_thread_num();
        GrtContext_t *c = context[g];
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
        omp_check_rc(grt_set_molecule_ppmv(c,
                                           H2O,
                                           &(ppmv[offset])),
                     rc[g]);

        /*Make up some more data for the column.*/
        for (j=0;j<num_levels;++j)
        {
            ppmv[i*num_levels+j] = 325. - 3.3*j;
        }

        /*Set the water vapor abundance for the library context.*/
        omp_check_rc(grt_set_molecule_ppmv(c,
                                           O3,
                                           &(ppmv[offset])),
                     rc[g]);

        /*Calculate the optical depths.*/
        omp_check_rc(grt_calculate_optical_depth(c,
                                                 &(pressure[offset]),
                                                 &(temperature[offset]),
                                                 &(optical_depth[offset*num_wpoints])),
                     rc[g]);
    }
    for (i=0;i<num_contexts;++i)
    {
        check_rc(rc[i]);
    }

    /*Clean up.*/
    free(pressure);
    free(temperature);
    free(ppmv);
    free(optical_depth);

    /*Free memory allocated by the library context.*/
    for (i=0;i<num_contexts;++i)
    {
        check_rc(grt_context_free(&(context[i])));
    }
}
