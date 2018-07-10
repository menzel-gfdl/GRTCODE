#ifndef QUERY_GPU_H_
#define QUERY_GPU_H_

#ifdef __NVCC__


__host__
int get_num_gpus(int * num_devices,
                 int const verbose);


#endif

#endif
