#include <stdio.h>
#include "cuda_helpers.cuh"

#ifndef __NVCC__
#error
#endif


cudaDeviceProp check_device_props()
{
    cudaDeviceProp p;
    int which_device;
    HANDLE_ERROR(cudaGetDevice(&which_device));
    HANDLE_ERROR(cudaGetDeviceProperties(&p,which_device));

    printf("\nFound Device %d : %s \n",
           which_device,
           p.name);
    printf("\tcompute capability : %d.%d \n",
           p.major,
           p.minor);
    printf("\tnumber of multiprocessors : %d \n",
           p.multiProcessorCount);

    int conc_kernel_support = p.concurrentKernels;
    if (conc_kernel_support == 0)
    {
        printf("\t%s does not support concurrent kernels\n",
               p.name);
    }
    else
    {
        printf("\t%s supports concurrent kernels\n",p.name);
    }

    return p;
}


#if defined(__CUDA_ARCH__)
#if __CUDA_ARCH__ < 600
__device__ double atomicAdd(double *address,
                            double val)
{
    unsigned long long int *address_as_ull = (unsigned long long int *)address;
    unsigned long long int old = *address_as_ull;
    unsigned long long int assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull,
                        assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif
#endif
