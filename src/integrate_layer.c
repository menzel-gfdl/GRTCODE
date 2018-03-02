#include "floating_point_type.h"
#include "integrate_layer.h"


#ifdef __NVCC__
__global__
void get_avg_TP(int const nLayer, /**<Number of atmospheric layers.*/
                fp_t const * const P, /**<Pressure at atmospheric layer interfaces [atm].*/
                fp_t const * const T, /**<Temperature at atmospheric layer interfaces [K].*/
                fp_t * const Pavg, /**<Average layer pressure [atm].*/
                fp_t * const Tavg) /**<Average layer temperature [K].*/
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < nLayer)
    {
        fp_t const HALF = 0.5;
        Pavg[tid] = HALF*(P[tid] + P[tid+1]);
        Tavg[tid] = HALF*(T[tid] + T[tid+1]);
    }
    return;
}
#endif


void get_avg_TP_h(int const nLayer, /**<Number of atmospheric layers.*/
                  fp_t const * const P, /**<Pressure at atmospheric layer interfaces [atm].*/
                  fp_t const * const T, /**<Temperature at atmospheric layer interfaces [K].*/
                  fp_t * const Pavg, /**<Average layer pressure [atm].*/
                  fp_t * const Tavg) /**<Average layer temperature [K].*/
{
    int tid;
    fp_t const HALF = 0.5;
    for (tid=0;tid<nLayer;++tid)
    {
        Pavg[tid] = HALF*(P[tid] + P[tid+1]);
        Tavg[tid] = HALF*(T[tid] + T[tid+1]);
    }
    return;
}


#ifdef __NVCC__
__global__
void get_avg_NPs(int const nlayer, /**<Number of atmospheric layers.*/
                 fp_t const * const x, /**<Molecular abundance at atmospheric layer interfaces.*/
                 fp_t const * const P, /**<Pressure at atmospheric layer interfaces [atm].*/
                 fp_t * const N, /**<Number density molecules integrated across the layer [cm^-2].*/
                 fp_t * const Psavg) /**<Average layer partial pressure [atm].*/
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < nlayer)
    {
        fp_t const HALF = 0.5;
        fp_t const THIRD = 1./3.;
        fp_t const SIXTH = 1./6.;
        fp_t const ATM_TO_PA = 101325;
        fp_t const MAIR = 0.02897/6.0221409e23; /*[kg]*/
        fp_t const G = 9.80665; /*[m s^-2]*/
        fp_t const M_TO_CM = 100;
        fp_t const c = ATM_TO_PA/(MAIR*G*M_TO_CM*M_TO_CM);
        Psavg[tid] = THIRD*(x[tid]*P[tid] + x[tid+1]*P[tid+1])
                     + SIXTH*(x[tid]*P[tid+1] + x[tid+1]*P[tid]);
        fp_t dp = P[tid] - P[tid+1];
        if (dp < 0.f)
        {
            dp *= -1.f;
        }
        N[tid] = dp*c*HALF*(x[tid]+x[tid+1]);
    }
    return;
}
#endif


void get_avg_NPs_h(int const nlayer, /**<Number of atmospheric layers.*/
                   fp_t const * const x, /**<Molecular abundance at atmospheric layer interfaces.*/
                   fp_t const * const P, /**<Pressure at atmospheric layer interfaces [atm].*/
                   fp_t * const N, /**<Number density molecules integrated across the layer [cm^-2].*/
                   fp_t * const Psavg) /**<Average layer partial pressure [atm].*/
{
    int tid;
    fp_t const HALF = 0.5;
    fp_t const THIRD = 1./3.;
    fp_t const SIXTH = 1./6.;
    fp_t const ATM_TO_PA = 101325;
    fp_t const MAIR = 0.02897/6.0221409e23; /*[kg]*/
    fp_t const G = 9.80665; /*[m s^-2]*/
    fp_t const M_TO_CM = 100;
    fp_t const c = ATM_TO_PA/(MAIR*G*M_TO_CM*M_TO_CM);
    for (tid=0;tid<nlayer;++tid)
    {
        Psavg[tid] = THIRD*(x[tid]*P[tid] + x[tid+1]*P[tid+1])
                     + SIXTH*(x[tid]*P[tid+1] + x[tid+1]*P[tid]);
        fp_t dp = P[tid] - P[tid+1];
        if (dp < 0.f)
        {
            dp *= -1.f;
        }
        N[tid] = dp*c*HALF*(x[tid]+x[tid+1]);
    }
    return;
}
