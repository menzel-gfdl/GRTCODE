#include "debug.h"
#include "floating_point_type.h"
#include "integrate_layer.h"


#define average(val1,val2) 0.5*(val1+val2);


/*Integrate the total number density of the air across each atmospheric layer,
  assuming that:
    - each layer is hydrostatic.
    - accerlation due to gravity is constant across each layer.
*/
#ifdef __NVCC__
__global__
void integrated_N(int const num_layers, /**<Number of atmospheric
                                            layers.*/
                  fp_t const * const P, /**<Pressure [atm] at
                                            atmospheric layer edges.*/
                  fp_t * const N) /**<Integrated layer number density
                                      [cm^-2].*/
{
    int const tid = blockIdx.x*blockDim.x + threadIdx.x;
    fp_t const G = 9.80665; /*Acceleration due to gravity [m s^-2].*/
    fp_t const M_AIR = 0.02897/6.0221409e23; /*Mass of an air molecule [kg].*/
    fp_t const ATM_TO_PA = 101325.; /*[Pa/atm].*/
    fp_t const M_TO_CM = 100.; /*[cm/m].*/
    fp_t const c = ATM_TO_PA/(M_AIR*G*M_TO_CM*M_TO_CM); /*[1/(atm*cm^2)].*/

    if (tid < num_layers)
    {
        fp_t dp = P[tid] - P[tid+1];
        if (dp < (fp_t)0.)
        {
            dp *= (fp_t)(-1.);
        }
        N[tid] = c*dp;
    }
    return;
}
#endif


int integrated_N_h(int const num_layers,
                   fp_t const * const P,
                   fp_t * const N)
{
    not_null(P);
    not_null(N);
    fp_t const G = 9.80665;
    fp_t const M_AIR = 0.02897/6.0221409e23;
    fp_t const ATM_TO_PA = 101325.;
    fp_t const M_TO_CM = 100.;
    fp_t const c = ATM_TO_PA/(M_AIR*G*M_TO_CM*M_TO_CM);
    int i;
    for (i=0;i<num_layers;++i)
    {
        fp_t dp = P[i] - P[i+1];
        if (dp < (fp_t)0.)
        {
            dp *= (fp_t)(-1.);
        }
        N[i] = c*dp;
    }
    return SUCCESS;
}


/*Calculate the Curtis-Godson integrals for pressure and temperature
  across each atmospheric layer, assuming that:
    - each layer is hydrostatic.
    - accerlation due to gravity is constant across each layer.
    - temperature is a linear function of pressure.
*/
#ifdef __NVCC__
__global__
void Curtis_Godson_PT(int const num_layers, /**<Number of atmospheric
                                                layers.*/
                      fp_t const * const P, /**<Pressure [atm] at
                                                atmospheric layer edges.*/
                      fp_t const * const T, /**<Temperature [K] at
                                                atmospheric layer edges.*/
                      fp_t * const Pavg, /**<Average layer pressure [atm].*/
                      fp_t * const Tavg) /**<Average layer temperature [K].*/
{
    int const tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num_layers)
    {
        Pavg[tid] = average(P[tid],P[tid+1]);
        Tavg[tid] = average(T[tid],T[tid+1]);
    }
    return;
}
#endif


int Curtis_Godson_PT_h(int const num_layers,
                       fp_t const * const P,
                       fp_t const * const T,
                       fp_t * const Pavg,
                       fp_t * const Tavg)
{
    not_null(P);
    not_null(T);
    not_null(Pavg);
    not_null(Tavg);
    int i;
    for (i=0;i<num_layers;++i)
    {
        Pavg[i] = average(P[i],P[i+1]);
        Tavg[i] = average(T[i],T[i+1]);
    }
    return SUCCESS;
}


/*Calculate the Curtis-Godson integrals for partial pressure and
  number density of a specific molecular species across each atmospheric
  layer, assuming that:
    - each layer is hydrostatic.
    - accerlation due to gravity is constant across each layer.
    - the molecular abundance is a linear function of pressure.
*/
#ifdef __NVCC__
__global__
void Curtis_Godson_PsNs(int const num_layers, /**<Number of atmospheric
                                                  layers.*/
                       fp_t const * const P, /**<Pressure [atm] at
                                                 atmospheric layer edges.*/
                       fp_t const * const x, /**<Molecular abundance at
                                                 atmospheric layer edges.*/
                       fp_t const * const N, /**<Total integrated layer
                                                 number density [cm^-2].*/
                       fp_t * const Psavg, /**<Average layer partial
                                               pressure [atm].*/
                       fp_t * const Ns) /**<Integrated number density [cm^-2]
                                            for the input molecule.*/
{
    int const tid = blockIdx.x*blockDim.x + threadIdx.x;
    fp_t const THIRD = 1./3.;
    fp_t const SIXTH = 1./6.;
    if (tid < num_layers)
    {
        Psavg[tid] = THIRD*(x[tid]*P[tid] + x[tid+1]*P[tid+1])
                     + SIXTH*(x[tid]*P[tid+1] + x[tid+1]*P[tid]);
        Ns[tid] = N[tid]*average(x[tid],x[tid+1]);
    }
    return;
}
#endif


int Curtis_Godson_PsNs_h(int const num_layers,
                         fp_t const * const P,
                         fp_t const * const x,
                         fp_t const * const N,
                         fp_t * const Psavg,
                         fp_t * const Ns)
{
    not_null(P);
    not_null(x);
    not_null(N);
    not_null(Psavg);
    not_null(Ns);
    fp_t const THIRD = 1./3.;
    fp_t const SIXTH = 1./6.;
    int i;
    for (i=0;i<num_layers;++i)
    {
        Psavg[i] = THIRD*(x[i]*P[i] + x[i+1]*P[i+1])
                   + SIXTH*(x[i]*P[i+1] + x[i+1]*P[i]);
        Ns[i] = N[i]*average(x[i],x[i+1]);
    }
    return SUCCESS;
}
