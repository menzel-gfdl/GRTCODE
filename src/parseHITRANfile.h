#ifndef PARSEHITRANFILE_H_
#define PARSEHITRANFILE_H_

#include "floating_point_type.h"

typedef struct line_params
{
    int *iso;
    fp_t *Vnn;
    fp_t *Snn_ref;
    float *Yair;
    float *Yself;
    float *En;
    float *n;
    float *d;
    int mol;
    unsigned int nLines;
} line_params_t;

typedef struct line_flags
{
    unsigned int cumemset_host_flags;
    unsigned int host : 1;
    unsigned int device : 1;
} line_flags_t;

int free_line_params_host(line_params_t **lineParams,
                          line_flags_t flags);

int parse_hitran_file(line_params_t **lineParams,
                      char *fileName,
                      line_flags_t flags,
                      double loWn,
                      double hiWn);

#ifdef __NVCC__
int alloc_line_params_device(line_params_t **lineParams,
                             unsigned int nLines);

int free_line_params_device(line_params_t **lineParams);
#endif

#endif
