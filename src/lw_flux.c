#include <math.h>
#include "floating_point_type.h"
#include "lw_flux.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#ifdef __NVCC__
__global__ void calc_lw_flux(unsigned int const nF,
                             int const numLayers,
                             fp_t * const fluxDown,
                             fp_t * const fluxUp,
                             fp_t const * const T,
                             fp_t const Tsurf,
                             fp_t const * const tau,
                             fp_t const w,
                             fp_t const res,
                             fp_t const emissivity,
                             fp_t const * const TLEV)
{
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    fp_t const wv = w + tid*res;
    fp_t c1[4];
    fp_t c2[4];
    int i;
    unsigned int j;
    fp_t I_down;
    fp_t tc;
    fp_t p;
    fp_t I_up;

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
                p = (1-tc)*effective_planck(T[i],
                                            TLEV[i+1],
                                            wv,
                                            tau[i*nF+tid]);
                I_down = p + I_down*tc;
                fluxDown[(i+1)*nF+tid] += c2[j]*I_down;
            }

            /*Upward pass.*/
            I_up = emissivity*planck_func(Tsurf,wv) + (1-emissivity)*I_down;
            fluxUp[numLayers*nF+tid] += c2[j]*I_up;

            for (i=numLayers-1;i>=0;--i)
            {
                tc = exp(c1[j]*tau[i*nF+tid]);
                p = (1-tc)*effective_planck(T[i],
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


void calc_lw_flux_h(unsigned int const nF,
                    int const numLayers,
                    fp_t * const fluxDown,
                    fp_t * const fluxUp,
                    fp_t const * const T,
                    fp_t const Tsurf,
                    fp_t const * const tau,
                    fp_t const w,
                    fp_t const res,
                    fp_t const emissivity,
                    fp_t const * const TLEV)
{
    /*Local variables*/
    fp_t c1[4];
    fp_t c2[4];
    int i;
    unsigned int j;
    unsigned int k;
    fp_t I_down[nF];
    fp_t tc;
    fp_t p;
    fp_t I_up[nF];

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
                p = (1-tc)*effective_planck(T[i],
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
            I_up[k] = emissivity*planck_func(Tsurf,w+res*k) +
                      (1-emissivity)*I_down[k];
            fluxUp[numLayers*nF+k] += c2[j]*I_up[k];
        }

        for (i=numLayers-1;i>=0;--i)
        {
            for (k=0;k<nF;k++)
            {
                tc = exp(c1[j]*tau[i*nF+k]);
                p = (1-tc)*effective_planck(T[i],
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
fp_t planck_func(fp_t const T,
                 fp_t const w)
{
    fp_t const h = 6.626070040E-34; /*Planck constant (J*s)*/
    fp_t const c = 299792458; /*Speed of light (m/s).*/
    fp_t const kB = 1.38064852E-23; /*Boltzmann constant (J/K).*/
    fp_t const MToCm = 100;
    fp_t const wm = w*MToCm; /*Wavenumber (1/m).*/

    return ((2*h*c*c*wm*wm*wm)/(exp(h*c*wm/(kB*T))-1));
}


#ifdef __NVCC__
__host__ __device__
#endif
fp_t effective_planck(fp_t const Tcenter,
                      fp_t const Tedge,
                      fp_t const w,
                      fp_t const tau)
{
    /*From equation 16 in Clough et al. (1992).*/
    fp_t const a = 0.193;
    fp_t const b = 0.013;

    return ((planck_func(Tcenter,w) + (a*tau + b*tau*tau)*
           planck_func(Tedge,w))/(1 + a*tau + b*tau*tau));
}
