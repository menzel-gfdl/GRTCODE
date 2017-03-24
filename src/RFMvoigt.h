#ifndef SET_RFMVOIGT_H
#define SET_RFMVOIGT_H

#ifdef __NVCC__
__host__ __device__
#endif
void voigt_shape_function(int const N,
                          double const * const DWNO,
                          double const WNOADJ,
                          float const DOPADJ,
                          float const WIDADJ,
                          float const STRADJ,
                          float *K);

#endif
