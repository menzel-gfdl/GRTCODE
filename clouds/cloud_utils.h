#ifndef CLOUD_UTILS_H
#define CLOUD_UTILS_H


#include "debug.h"


int grid_band_mapping(int const grid_size,
                      fp_t const w0,
                      fp_t const dw,
                      int const num_bands,
                      fp_t * const band_limits,
                      int * mapping);


#ifdef __NVCC__
__global__ void grid_band_mapping_d(int const grid_size,
                                    fp_t const w0,
                                    fp_t const dw,
                                    int const num_bands,
                                    fp_t * const band_limits,
                                    int * mapping);
#endif


int process_optics(int const grid_size,
                   int const num_layers,
                   int const num_bands,
                   int const * mapping,
                   fp_t const * tau,
                   fp_t const * omega,
                   fp_t const * g,
                   fp_t * optics_tau,
                   fp_t * optics_omega,
                   fp_t * optics_g);


#ifdef __NVCC__
__global__ void process_optics_d(int const grid_size,
                                 int const num_layers,
                                 int const num_bands,
                                 int const * mapping,
                                 fp_t const * tau,
                                 fp_t const * omega,
                                 fp_t const * g,
                                 fp_t * optics_tau,
                                 fp_t * optics_omega,
                                 fp_t * optics_g);
#endif


#endif
