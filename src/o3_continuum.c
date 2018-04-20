#include <stdlib.h>
#include "omp.h"
#include "continuum_helpers.h"
#include "debug.h"
#include "floating_point_type.h"
#include "o3_continuum.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif


/*Read in the ozone continuum coefficients.*/
int get_ozone_continuum_coefs(char const * const filepath,
                              OzoneContinuumCoefs_t *cc,
                              unsigned int const nws,
                              int const w0,
                              double const res,
                              int put_on_device)
{
    not_null(cc);
    enum ozone_csv_coefs
    {
        CROSS_SECTION = 0,
        CSV_NUM_COEFS
    };
    fp_t *coefs[CSV_NUM_COEFS];
    log_mesg("Reading in ozone continuum coefficients from file %s.",
             filepath);
    check(get_coefs(filepath,
                    coefs,
                    CSV_NUM_COEFS,
                    nws,
                    w0,
                    res));
    cc->cross_section = coefs[CROSS_SECTION];
    if (put_on_device)
    {
        using_gpu();
#ifdef __NVCC__
        fatal("implement this on the device.");
#endif
    }
    return SUCCESS;
}


int free_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc,
                               int const on_device)
{
    not_null(cc);
    if (on_device)
    {
        using_gpu();
#ifdef __NVCC__
        fatal("implement this on the device.");
#endif
    }
    else
    {
        free(cc->cross_section);
    }
    return SUCCESS;
}


#ifdef __NVCC__
__global__
void calc_ozone_ctm_optdepth(unsigned int const nws,
                             int const nlayers,
                             fp_t const * const cross_section,
                             fp_t const * const N,
                             fp_t * const tau)
{
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < nws)
    {
        int lyr;

#pragma unroll
        for (lyr=0;lyr<nlayers;++lyr)
        {
            tau[lyr*nws+tid] += N[lyr]*cross_section[tid];
        }
    }
    return;
}
#endif


void calc_ozone_ctm_optdepth_h(unsigned int const nws,
                               int const nlayers,
                               fp_t const * const cross_section,
                               fp_t const * const N,
                               fp_t * const tau)
{
    int lyr;
    unsigned int tid;

#pragma omp parallel for collapse(2) \
                         schedule(static) \
                         default(none) \
                         private(lyr,tid)
    for (lyr=0;lyr<nlayers;++lyr)
    {
        for (tid=0;tid<nws;++tid)
        {
            tau[lyr*nws+tid] += N[lyr]*cross_section[tid];
        }
    }
    return;
}
