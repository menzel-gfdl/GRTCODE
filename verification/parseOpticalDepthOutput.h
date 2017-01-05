#ifndef SET_PARSEOPTICALDEPTHOUTPUT_H_
#define SET_PARSEOPTICALDEPTHOUTPUT_H_

#include <stdlib.h>
#include "myreal.h"

typedef struct OpticalDepth_t
{
    size_t ntime;         /*Number of time grid points.*/
    size_t nlat;          /*Number of latitude grid points.*/
    size_t nlon;          /*Number of longitude grid points.*/
    size_t npfull;        /*Number of pressure layers.*/
    size_t nwavenumber;   /*Number of wavenumber grid points.*/
    REAL_t *OpticalDepth; /*Optical depth (time,lat,lon,pfull,wavenumber).*/
} OpticalDepth_t;

int readOpticalDepthFromGrtcodeFile(char fname[],
                                    OpticalDepth_t *GrtOutput);

int readOpticalDepthFromRfmFile(char fname[],
                                OpticalDepth_t *RfmOutput,
                                size_t num_layers,
                                size_t num_wavenumbers);

int freeOpticalDepth(OpticalDepth_t *OpticalDepth);

#endif
