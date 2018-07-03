#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "new.h"


#define check_rc(rc) { \
    if (rc != 0) {fprintf(stderr,"%s: %d Error\n",__FILE__,__LINE__);return 1;}}


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
    uint64_t num_wpoints;
    int use_gpu = 0; /*0 = host only run. use_gpu = 1 will attempt to run
                       on the default GPU, but requires the program to be
                       compiled with nvcc.*/
    int use_h2o_ctm = 1;
    int use_o3_ctm = 1;
    check_rc(initialize_grt(&context,
                            num_levels,
                            w0,
                            wn,
                            wres,
                            &num_wpoints,
                            NULL,
                            &use_gpu,
                            NULL,
                            &use_h2o_ctm,
                            &use_o3_ctm));

    /*Add water vapor.*/
    char h2o_hitran[64];
    snprintf(h2o_hitran,64,"h2o_hit12.par");
    int h2o;
    double min_line_center_wavenumber = 4.;
    double max_line_center_wavenumber = 8.;
    check_rc(add_molecule(context, 
                          h2o_hitran, 
                          &h2o, 
                          &min_line_center_wavenumber,
                          &max_line_center_wavenumber));

    /*Mimic looping over columns.*/
    FP_t *pressure = (FP_t *)malloc(sizeof(*pressure)*num_levels);
    FP_t *temperature = (FP_t *)malloc(sizeof(*temperature)*num_levels);
    FP_t *ppmv = (FP_t *)malloc(sizeof(*ppmv)*num_levels);
    FP_t *optical_depth = (FP_t *)malloc(sizeof(*optical_depth)*
                                         num_wpoints*(num_levels-1));
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
        check_rc(set_molecule_ppmv(context, 
                                   h2o, 
                                   ppmv));

        /*Calculate the optical depths.*/
        check_rc(calculate_optical_depth(context, 
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
    check_rc(finalize_grt(&context));
}
