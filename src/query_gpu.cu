#include <stdio.h>
#include "debug.h"

#ifdef __NVCC__

#include "cuda_helpers.cuh"


__host__
int get_num_gpus(int * num_devices)
{
    not_null(num_devices);
    HANDLE_ERROR(cudaGetDeviceCount(num_devices));
    log_info("Found %d CUDA-enabled GPUS:",
             *num_devices);
    int i;
    for (i=0;i<(*num_devices);++i)
    {
        cudaDeviceProp prop;
        HANDLE_ERROR(cudaGetDeviceProperties(&prop,i));
        log_info("\tDevice number: %d (%s)",
                 i,
                 prop.name);
    }
    return SUCCESS;
}


#endif
