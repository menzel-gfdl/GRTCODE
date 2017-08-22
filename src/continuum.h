#ifndef SET_CONTINUUM_H_
#define SET_CONTINUUM_H_

#include "myreal.h"

void parseCKD(const char fname[],
              REAL_t *AryPtr,
              const int maxwavenum,
              const int minw,
              REAL_t const res);

#ifdef __NVCC__
__global__
void calc_ctm_optdepth(unsigned int const nF,
                       unsigned int const numLayers,
                       REAL_t * const optdepth,
                       REAL_t const * const CS,
                       REAL_t const * const T,
                       REAL_t const * const PS_H2O,
                       REAL_t const * const Z,
                       REAL_t const * const T0,
                       REAL_t const * const CF,
                       REAL_t const * const P,
                       REAL_t const * const T0F);
#endif

void calc_ctm_optdepth_h(unsigned int const nF,
                         unsigned int const numLayers,
                         REAL_t * const optdepth,
                         REAL_t const * const CS,
                         REAL_t const * const T,
                         REAL_t const * const PS_H2O,
                         REAL_t const * const Z,
                         REAL_t const * const T0,
                         REAL_t const * const CF,
                         REAL_t const * const P,
                         REAL_t const * const T0F);

#endif
