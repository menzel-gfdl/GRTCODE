#ifndef CUDA_HELPERS_H_
#define CUDA_HELPERS_H_

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef __NVCC__
#error
#endif

#define HANDLE_ERROR(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

cudaDeviceProp check_device_props();

#endif

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* address, double val);
#endif
