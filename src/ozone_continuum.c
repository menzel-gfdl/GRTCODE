#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "debug.h"
#include "floating_point_type.h"
#include "ozone_continuum.h"
#include "parse_csv.h"
#include "utils.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif


/*Read in the ozone continuum coefficients.*/
int get_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc,
                              unsigned int const nws,
                              int const w0,
                              double const res)
{
    not_null(cc);

    /*Set file name.*/
    char *filepath = "INPUT/ozone_continuum/ozone_continuum.csv";
    int num_vals = 1;

    /*Read in the data.*/
    log_mesg("Reading in ozone continuum coefficients from file %s.",
             filepath);
    int num_lines;
    int num_cols;
    char **buf;
    check(parse_csv(filepath,
                    &num_lines,
                    &num_cols,
                    1,
                    &buf));
    if ((num_vals + 1) != num_cols)
    {
        fatal("The number of columns (%d) in file %s does not match"
                 " the expected number (%d).",
              num_cols,
              filepath,
              num_vals + 1);
    }

    /*Convert the data from strings to floating point.*/
    fp_t *fbuf = NULL;
    int data_size = num_lines*num_cols;
    check(malloc_ptr((void **)(&fbuf),
                     sizeof(*fbuf)*data_size));
    int j;
    for (j=0;j<data_size;++j)
    {
        double d;
        check(to_double(buf[j],
                        &d));
        check(to_fp_t(d,
                      &(fbuf[j])));
        free(buf[j]);
    }
    free(buf);

    /*Allocate space for the coefficient values at each wavenumber.*/
    fp_t *c = NULL;
    int num_bytes = sizeof(*c)*nws;
    check(malloc_ptr((void **)&c,
                     num_bytes));
    memset(c,
           0,
           num_bytes);

    /*Interpolate to wavenumber grid.*/
    fp_t *x = &(fbuf[0]);
    fp_t *y = &(fbuf[num_lines]);
    for (j=0;j<nws;++j)
    {
        check(linear_interpolation(x,
                                   y,
                                   num_lines,
                                   (fp_t)(w0 + j*res),
                                   &(c[j])));
    }
    free(fbuf);
    cc->cross_section = c;
    cc->nws = nws;
    return SUCCESS;
}


int free_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc)
{
    not_null(cc);
    free(cc->cross_section);
    cc->cross_section = NULL;
    return SUCCESS;
}


int put_ozone_coefs_on_device(OzoneContinuumCoefs_t const * const in,
                              OzoneContinuumCoefs_t * const out)
{
    not_null(in);
    not_null(out);
    using_gpu();
    int num_bytes = sizeof(*(in->cross_section))*(in->nws);
    HANDLE_ERROR(cudaMalloc(&(out->cross_section),
                            num_bytes));
    HANDLE_ERROR(cudaMemcpy(out->cross_section,
                            in->cross_section,
                            num_bytes,
                            cudaMemcpyHostToDevice));
    return SUCCESS;
}


int remove_ozone_coefs_from_device(OzoneContinuumCoefs_t * const in)
{
    not_null(in);
    using_gpu();
    HANDLE_ERROR(cudaFree(in->cross_section));
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
