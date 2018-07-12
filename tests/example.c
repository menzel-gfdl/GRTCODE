#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "new.h"


#define check_rc(rc) { \
    if (rc != 0) { \
        char buf[256]; \
        int e_ = grt_errstr(rc,buf,256); \
        fprintf(stderr,"%s\n",buf); \
        fprintf(stderr,"%s: %d Error\n",__FILE__,__LINE__); \
        return EXIT_FAILURE; \
    }}


#ifdef DOUBLE_PRECISION
#define FP_t double
#else
#define FP_t float
#endif


int main(void)
{
    /*Initalize library.*/
    GrtContext_t *context;
    int num_levels = 25;
    double w0 = 2.;
    double wn = 100.;
    double wres = 0.1;
    char *h2o_ctm_dir = "water_vapor_continuum";
    char *o3_ctm_dir = "ozone_continuum";
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

    /*Add water vapor.*/
    char h2o_hitran[64];
    snprintf(h2o_hitran,64,"HITRAN_files/01_hit12.par");
    int h2o;
    double min_line_center_wavenumber = 4.;
    double max_line_center_wavenumber = 8.;
    check_rc(grt_add_molecule(context, 
                              h2o_hitran, 
                              &h2o, 
                              &min_line_center_wavenumber,
                              &max_line_center_wavenumber));

    /*Mimic looping over columns.*/
    uint64_t num_wpoints;
    check_rc(grt_get_spectral_grid_size(context,
                                        &num_wpoints));
    int num_layers = num_levels - 1;
    FP_t *pressure = (FP_t *)malloc(sizeof(*pressure)*num_levels);
    FP_t *temperature = (FP_t *)malloc(sizeof(*temperature)*num_levels);
    FP_t *ppmv = (FP_t *)malloc(sizeof(*ppmv)*num_levels);
    FP_t *optical_depth = (FP_t *)malloc(sizeof(*optical_depth)*
                                         num_wpoints*num_layers);
    int num_columns = 4;
    int i;
    for (i=0;i<num_columns;++i)
    {
        /*Make up some data.*/
        int j;
        for (j=0;j<num_levels;++j)
        {
            pressure[j] = 0.1 + 150.*j;
            temperature[j] = 230. + 2.3*j;
            ppmv[j] = 300. + 0.2*j;
        }

        /*Set water vapor ppmv.*/
        check_rc(grt_set_molecule_ppmv(context, 
                                       h2o, 
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

    /*Finalize library.*/
    check_rc(grt_context_free(&context));
}
