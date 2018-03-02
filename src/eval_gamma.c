#include "eval_gamma.h"
#include "floating_point_type.h"
#include "LorentzFuncs.h"
#include "omp.h"


/*Compute the pressure broadened line halfwidth for each transition.

  Arguments:
      numLayers [in]      Size of the height dimension for the inputted
                              arrays.
      nL        [in]      Size of the line dimension for the inputted arrays.
      P         [in]      Array of pressures (atm).  This array is stored as
                              [height].
      T         [in]      Array of temperatures (K).  This array is stored as
                              [height].
      Ps        [in]      Array of partial pressures (atm).  This array is
                              stored as [height].
      Yself     [in]      Array of self-broadened halfwidth at half
                              maximum (cm^-1*atm^-1).  This array is stored as
                              [line].
      Yair      [in]      Array of air-broadened halfwidth at half
                              maximum (cm^-1*atm^-1).  This array is stored as
                              [line].
      n         [in]      Array of coefficients of temperature dependence of
                              the air-broadened halfwidth at half maximum.
                              This array is stored as [line].
      Gam       [in,out]  Array of pressure broadened line halfwidths (cm^-1).
                              This array is stored as [height][line].
*/
#ifdef __NVCC__
__global__
void eval_gamma(int const numLayers,
                unsigned int const nL,
                fp_t const * const P,
                fp_t const * const T,
                fp_t const * const Ps,
                float const * const Yself,
                float const * const Yair,
                float const * const n,
                fp_t * const Gam)
{
    int lyr;
    unsigned int ltid = blockIdx.x*blockDim.x + threadIdx.x;
    if (ltid < nL)
    {
#pragma unroll
        for (lyr=0;lyr<numLayers;++lyr)
        {
            Gam[lyr*nL+ltid] = lorentz_hwhm(P[lyr],
                                            T[lyr],
                                            Yself[ltid],
                                            Yair[ltid],
                                            n[ltid],
                                            Ps[lyr]);
        }
    }
    return;
}

#endif


void eval_gamma_h(int const numLayers,
                  unsigned int const nL,
                  fp_t const * const P,
                  fp_t const * const T,
                  fp_t const * const Ps,
                  float const * const Yself,
                  float const * const Yair,
                  float const * const n,
                  fp_t * const Gam)
{
    int lyr;
    unsigned int ltid;

#pragma omp parallel for schedule(static) \
                         collapse(2) \
                         default(none) \
                         private(ltid,lyr)
/*
                         shared(nL,numLayers,Gam,P,T,Yself, \
                                Yair,n,Ps) \
*/
    for (lyr=0;lyr<numLayers;++lyr)
    {
        for(ltid=0;ltid<nL;++ltid)
        {
            Gam[lyr*nL+ltid] = lorentz_hwhm(P[lyr],
                                            T[lyr],
                                            Yself[ltid],
                                            Yair[ltid],
                                            n[ltid],
                                            Ps[lyr]);
        }
    }
    return;
}
