#include <math.h>
#include "flux.h"
#include "myreal.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __NVCC__
__global__ void calcFlux(unsigned int const nF,
                         unsigned int const numLayers,
                         REAL_t * const fluxDown,
                         REAL_t * const fluxUp,
                         REAL_t const * const T,
                         REAL_t const Tsurf,
                         REAL_t const * const tau,
                         REAL_t const w,
                         REAL_t const res,
                         REAL_t const emissivity,
                         REAL_t const * const TLEV)
{
    /*Local variables*/
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    REAL_t const wv = w + tid*res;
    REAL_t c1[4];
    REAL_t c2[4];
    int i;
    unsigned int j;
    REAL_t I_down;
    REAL_t tc;
    REAL_t p;
    REAL_t I_up;

    if (tid < nF)
    {
        /*Set constants.*/
        c1[0] = -1./0.0694318442;
        c1[1] = -1./0.3300094782;
        c1[2] = -1./0.6699905218;
        c1[3] = -1./0.9305681558;

        c2[0] = 2*M_PI*0.0694318442*0.1739274226;
        c2[1] = 2*M_PI*0.3300094782*0.3260725774;
        c2[2] = 2*M_PI*0.6699905218*0.3260725774;
        c2[3] = 2*M_PI*0.9305681558*0.1739274226;

        /*Zero out the flux arrays.*/
        for (i=0;i<numLayers+1;++i)
        {
            fluxDown[i*nF+tid] = 0;
            fluxUp[i*nF+tid] = 0;
        }

        for (j=0;j<4;++j)
        {
            /*Downward pass.*/
            I_down = 0;

            for (i=0;i<numLayers;++i)
            {
                tc = exp(c1[j]*tau[i*nF+tid]);
                p = (1-tc)*effectivePlanck(T[i],
                                           TLEV[i+1],
                                           wv,
                                           tau[i*nF+tid]);
                I_down = p + I_down*tc;
                fluxDown[(i+1)*nF+tid] += c2[j]*I_down;
            }

            /*Upward pass.*/
            I_up = emissivity*planckFunc(Tsurf,wv) + (1-emissivity)*I_down;
            fluxUp[numLayers*nF+tid] += c2[j]*I_up;

            for (i=numLayers-1;i>=0;--i)
            {
                tc = exp(c1[j]*tau[i*nF+tid]);
                p = (1-tc)*effectivePlanck(T[i],
                                           TLEV[i],
                                           wv,
                                           tau[i*nF+tid]);
                I_up = p + I_up*tc;
                fluxUp[i*nF+tid] += c2[j]*I_up;
            }
        }
    }

    return;
}
#endif

void calcFlux_h(unsigned int const nF,
                unsigned int const numLayers,
                REAL_t * const fluxDown,
                REAL_t * const fluxUp,
                REAL_t const * const T,
                REAL_t const Tsurf,
                REAL_t const * const tau,
                REAL_t const w,
                REAL_t const res,
                REAL_t const emissivity,
                REAL_t const * const TLEV)
{
    /*Local variables*/
    REAL_t c1[4];
    REAL_t c2[4];
    int i;
    unsigned int j;
    unsigned int k;
    REAL_t I_down[nF];
    REAL_t tc;
    REAL_t p;
    REAL_t I_up[nF];

    /*Set constants.*/
    c1[0] = -1./0.0694318442;
    c1[1] = -1./0.3300094782;
    c1[2] = -1./0.6699905218;
    c1[3] = -1./0.9305681558;

    c2[0] = 2*M_PI*0.0694318442*0.1739274226;
    c2[1] = 2*M_PI*0.3300094782*0.3260725774;
    c2[2] = 2*M_PI*0.6699905218*0.3260725774;
    c2[3] = 2*M_PI*0.9305681558*0.1739274226;

    /*Zero out the flux arrays.*/
    for (i=0;i<numLayers+1;++i)
    {
        for (j=0;j<nF;++j)
        {
            fluxDown[i*nF+j] = 0;
            fluxUp[i*nF+j] = 0;
        }
    }

    for (j=0;j<4;++j)
    {
        /*Downward pass.*/
        for (k=0;k<nF;++k)
        {
            I_down[k] = 0;
        }

        for (i=0;i<numLayers;++i)
        {
            for (k=0;k<nF;++k)
            {
                tc = exp(c1[j]*tau[i*nF+k]);
                p = (1-tc)*effectivePlanck(T[i],
                                           TLEV[i+1],
                                           w+res*k,
                                           tau[i*nF+k]);
                I_down[k] = p + I_down[k]*tc;
                fluxDown[(i+1)*nF+k] += c2[j]*I_down[k];
            }
        }

        /*Upward pass.*/
        for (k=0;k<nF;++k)
        {
            I_up[k] = emissivity*planckFunc(Tsurf,w+res*k) +
                      (1-emissivity)*I_down[k];
            fluxUp[numLayers*nF+k] += c2[j]*I_up[k];
        }

        for (i=numLayers-1;i>=0;--i)
        {
            for (k=0;k<nF;k++)
            {
                tc = exp(c1[j]*tau[i*nF+k]);
                p = (1-tc)*effectivePlanck(T[i],
                                           TLEV[i],
                                           w+res*k,
                                           tau[i*nF+k]);
                I_up[k] = p + I_up[k]*tc;
                fluxUp[i*nF+k] += c2[j]*I_up[k];
            }
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

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t effectivePlanck(REAL_t const Tcenter,
                       REAL_t const Tedge,
                       REAL_t const w,
                       REAL_t const tau)
{
    /*From equation 16 in Clough et al. (1992).*/
    REAL_t const a = 0.193;
    REAL_t const b = 0.013;

    return ((planckFunc(Tcenter,w) + (a*tau + b*tau*tau)*planckFunc(Tedge,w))/
               (1 + a*tau + b*tau*tau));
}

void sum_fluxes(unsigned int const nF,
                unsigned int const numLevels,
                REAL_t const * const fluxes,
                REAL_t * const fluxes_accumulated,
                REAL_t const res)
{
    REAL_t const MToCm = 100;
    REAL_t const resm = res*MToCm; /*Wavenumber resolution (1/m).*/
    unsigned int i;
    unsigned int j;

    for (i=0;i<numLevels;++i)
    {
        fluxes_accumulated[i] = 0;
        for (j=0;j<nF-1;++j)
        {
            REAL_t a = fluxes[i*nF+j];
            REAL_t b = fluxes[i*nF+j+1] - a;
            fluxes_accumulated[i] += resm*(a+0.5*b);
        }
    }

    return;
}

/*
Old versions.  Delete when confident in new ones.

#ifdef __NVCC__
__global__ void calcFlux(unsigned int const nF,
                         unsigned int const numLayers,
                         REAL_t * const fluxDown,
                         REAL_t * const fluxUp,
                         REAL_t const * const T,
                         REAL_t const Tsurf,
                         REAL_t const * const tau,
                         REAL_t const w,
                         REAL_t const res)
{
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
        I_down = 0;
        I_up = planckFunc(Tsurf,wv);

        fluxDown[tid] = M_PI*I_down;
        fluxUp[numLayers*nF+tid] = M_PI*I_up;

#pragma unroll
        for (i=0;i<numLayers;++i)
        {
            tc_down = exp(coef*tau[i*nF+tid]);
            p_down = planckFunc(T[i],wv)*(1-tc_down);
            I_down = p_down + I_down*tc_down;
            fluxDown[(i+1)*nF+tid] = M_PI*I_down;

            up_index = numLayers - 1 - i;
            tc_up = exp(coef*tau[up_index*nF+tid]);
            p_up = planckFunc(T[up_index],wv)*(1-tc_up);
            I_up = p_up + I_up*tc_up;
            fluxUp[up_index*nF+tid] = M_PI*I_up;
        }
    }

    return;
}

#endif

void calcFlux_h(unsigned int const nF,
                unsigned int const numLayers,
                REAL_t * const fluxDown,
                REAL_t * const fluxUp,
                REAL_t const * const T,
                REAL_t const Tsurf,
                REAL_t const * const tau,
                REAL_t const w,
                REAL_t const res)
{
    unsigned int i;
    REAL_t I_down[nF];
    REAL_t I_up[nF];
    unsigned int j;
    REAL_t tc_down;
    REAL_t p_down;
    REAL_t tc_up;
    REAL_t p_up;
    int up_index;
    REAL_t const coef = -1.66;

    for (i=0;i<nF;++i)
    {
        I_down[i] = 0;
        I_up[i] = planckFunc(Tsurf,w+i*res);
        fluxDown[i] = M_PI*I_down[i];
        fluxUp[numLayers*nF+i] = M_PI*I_up[i];
    }

    for (i=0;i<numLayers;++i)
    {
        up_index = numLayers - 1 - i;

        for (j=0;j<nF;++j)
        {
            tc_down = exp(coef*tau[i*nF+j]);
            p_down = planckFunc(T[i],w+j*res)*(1-tc_down);
            I_down[j] = p_down + I_down[j]*tc_down;
            fluxDown[(i+1)*nF+j] = M_PI*I_down[j];

            tc_up = exp(coef*tau[up_index*nF+j]);
            p_up = planckFunc(T[up_index],w+j*res)*(1-tc_up);
            I_up[j] = p_up + I_up[j]*tc_up;
            fluxUp[up_index*nF+j] = M_PI*I_up[j];
        }
    }

    return;
}
*/
