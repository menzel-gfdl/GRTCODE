#include <math.h>
#include "flux.h"
#include "myreal.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __NVCC__
__global__ void calcFlux(unsigned int const nF,
                         unsigned int const numLayers,
                         REAL_t * const fluxOut,
                         REAL_t const * const T,
                         REAL_t const Tsurf,
                         REAL_t const * const tau,
                         REAL_t const w,
                         REAL_t const res)
{
    /*Local variables*/
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    REAL_t I_down;
    REAL_t I_up;
    REAL_t const wv = w + tid*res;
    unsigned int i;
    REAL_t tc_down;
    REAL_t const coef = -1.66;
    REAL_t p_down;
    REAL_t tc_up;
    REAL_t p_up;
    int up_index;

    if (tid < nF)
    {
        /*No longwave intensity from space.*/
        I_down = 0;

        /*Longwave intensity from the surface.*/
        I_up = planckFunc(Tsurf,wv);

        fluxOut[tid] = M_PI*I_down;
        fluxOut[numLayers*nF+tid] = M_PI*I_up;

#pragma unroll
        for (i=0;i<numLayers;++i)
        {
            tc_down = exp(coef*tau[i*nF+tid]);
            p_down = planckFunc(T[i],wv)*(1-tc_down);
            I_down = p_down + I_down*tc_down;
            fluxOut[(i+1)*nF+tid] += M_PI*I_down;

            up_index = numLayers - 1 - i;
            tc_up = exp(coef*tau[up_index*nF+tid]);
            p_up = planckFunc(T[up_index],wv)*(1-tc_up);
            I_up = p_up + I_up*tc_up;
            fluxOut[up_index*nF+tid] += M_PI*I_up;
        }
    }

    return;
}

#endif

void calcFlux_h(unsigned int const nF,
                unsigned int const numLayers,
                REAL_t * const fluxOut,
                REAL_t const * const T,
                REAL_t const Tsurf,
                REAL_t const * const tau,
                REAL_t const w,
                REAL_t const res)
{
    /*Local variables*/
    unsigned int i;
    REAL_t I_down[nF];
    REAL_t I_up[nF];
    unsigned int j;
    REAL_t tc;
    REAL_t coef = -1.66;
    REAL_t p;

    for (i=0;i<nF;++i)
    {
        I_down[i] = 0;
        I_up[i] = planckFunc(Tsurf,w+i*res);
        fluxOut[i] = M_PI*I_down[i];
        fluxOut[numLayers*nF+i] = M_PI*I_up[i];
    }

    for (i=0;i<numLayers;++i)
    {
        for (j=0;j<nF;++j)
        {
            tc = exp(coef*tau[i*nF+j]);
            p = planckFunc(T[i],w+j*res)*(1-tc);

            I_down[j] = p + I_down[j]*tc;
            I_up[j] = p + I_up[j]*tc;

            fluxOut[(i+1)*nF+j] += M_PI*I_down[j];
            fluxOut[(numLayers-1-i)*nF+j] += M_PI*I_up[j];
        }
    }

    return;
}

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t planckFunc(REAL_t const T,
                  REAL_t const w)
{
    /*Local variables*/
    REAL_t const h = 6.626070040E-34; /*Planck constant (J*s)*/
    REAL_t const c = 299792458; /*Speed of light (m/s).*/
    REAL_t const kB = 1.38064852E-23; /*Boltzmann constant (J/K).*/
    REAL_t const MToCm = 100;
    REAL_t const wm = w*MToCm; /*Wavenumber (1/m).*/

    return ((2*h*c*c*wm*wm*wm)/(exp(h*c*wm/(kB*T))-1));
}

void sum_fluxes(unsigned int const nF,
                unsigned int const numLevels,
                REAL_t const * const fluxes,
                REAL_t * const fluxes_accumulated,
                REAL_t const res)
{
    REAL_t const MToCm = 100;
    REAL_t const resm = res*MToCm; /*Wavenumber resolution (m).*/

    for (unsigned int i=0;i<numLevels;++i)
    {
        fluxes_accumulated[i] = 0;
        for (unsigned int j=0;j<nF-1;++j)
        {
            REAL_t a = fluxes[i*nF+j];
            REAL_t b = fluxes[i*nF+j+1] - a;
            fluxes_accumulated[i] += resm*(a+0.5*b);
        }
    }

    return;
}

#ifndef SKIPMAIN
#include <stdlib.h>

int main(int argc,char **argv)
{


    return EXIT_SUCCESS
}

#endif
