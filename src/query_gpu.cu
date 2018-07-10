#include <stdio.h>
#include "debug.h"

#ifdef __NVCC__

#include "cuda_helpers.cuh"


__host__
int get_num_gpus(int * num_devices,
                 int const verbose)
{
    not_null(num_devices);
    HANDLE_ERROR(cudaGetDeviceCount(num_devices));
    if (verbose)
    {
        int i;
        for (i=0;i<(*num_devices);++i)
        {
            cudaDeviceProp prop;
            HANDLE_ERROR(cudaGetDeviceProperties(&prop,i));
            fprintf(stderr,
                    "Device number: %d\nDevice name: %s\n\n",
                    i,
                    prop.name);
        }
    }
    return SUCCESS;
}


#endif
