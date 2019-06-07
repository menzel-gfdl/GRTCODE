#include "debug.h"
#include "floating_point_type.h"
#include "curtis_godson.h"


/*Calculate integrated number densities.*/
int calc_number_densities(int const num_layers, fp_t const * const p, fp_t * const n)
{
    fp_t const c = 2.147822334314468e+25; /*[1/(cm^2*atm)]*/
    int i;
#pragma omp parallel for default(none) private(i)
    for (i=0; i<num_layers; ++i)
    {
        fp_t dp = p[i] - p[i+1];
        dp = dp >= 0.f ? dp : dp*-1.f;
        n[i] = c*dp;
    }
    return RS_SUCCESS;
}


#ifdef __NVCC__
/*Calculate integrated number densities.*/
__global__ void calc_number_densities_d(int const num_layers, fp_t const * const p,
                                        fp_t * const n)
{
    int const i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < num_layers)
    {
        fp_t const c = 2.147822334314468e+25; /*[1/(cm^2*atm)]*/
        fp_t dp = p[i] - p[i+1];
        dp = dp >= 0.f ? dp : dp*-1.f;
        n[i] = c*dp;
    }
    return;
}
#endif


/*Calculate layer pressures and temperatures.*/
int calc_pressures_and_temperatures(int const num_layers, fp_t const * const p,
                                    fp_t const * const t, fp_t * const pavg,
                                    fp_t * const tavg)
{
    int i;
#pragma omp parallel for default(none) private(i)
    for (i=0;i<num_layers;++i)
    {
        pavg[i] = 0.5f*(p[i] + p[i+1]);
        tavg[i] = 0.5f*(t[i] + t[i+1]);
    }
    return RS_SUCCESS;
}


#ifdef __NVCC__
/*Calculate layer pressures and temperatures.*/
__global__ void calc_pressures_and_temperatures_d(int const num_layers, fp_t const * const p,
                                                  fp_t const * const t, fp_t * const pavg,
                                                  fp_t * const tavg)
{
    int const i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < num_layers)
    {
        pavg[i] = 0.5f*(p[i] + p[i+1]);
        tavg[i] = 0.5f*(t[i] + t[i+1]);
    }
    return;
}
#endif


/*Calculate partial pressures and number densities.*/
int calc_partial_pressures_and_number_densities(int const num_layers, fp_t const * const p,
                                                fp_t const * const x, fp_t const * const n,
                                                fp_t * const ps, fp_t * const ns)
{
    fp_t const third = 1.f/3.f;
    fp_t const sixth = 1.f/6.f;
    int i;
#pragma omp parallel for default(none) private(i)
    for (i=0;i<num_layers;++i)
    {
        ps[i] = third*(x[i]*p[i] + x[i+1]*p[i+1]) + sixth*(x[i]*p[i+1] + x[i+1]*p[i]);
        ns[i] = n[i]*0.5f*(x[i] + x[i+1]);
    }
    return RS_SUCCESS;
}


#ifdef __NVCC__
/*Calculate partial pressures and number densities.*/
__global__ void calc_partial_pressures_and_number_densities_d(int const num_layers, fp_t const * const p,
                                                              fp_t const * const x, fp_t const * const n,
                                                              fp_t * const ps, fp_t * const ns)
{
    int const i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < num_layers)
    {
        fp_t const third = 1.f/3.f;
        fp_t const sixth = 1.f/6.f;
        ps[i] = third*(x[i]*p[i] + x[i+1]*p[i+1]) + sixth*(x[i]*p[i+1] + x[i+1]*p[i]);
        ns[i] = n[i]*0.5f*(x[i] + x[i+1]);
    }
    return;
}
#endif
