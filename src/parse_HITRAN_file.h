#ifndef PARSE_HITRAN_FILE_H_
#define PARSE_HITRAN_FILE_H_

#include "floating_point_type.h"


typedef struct LineParams
{
    int *iso;
    fp_t *vnn;
    fp_t *snn_ref;
    float *yair;
    float *yself;
    float *en;
    float *n;
    float *d;
    int mol;
    unsigned int num_lines;
} LineParams_t;


typedef struct LineFlags
{
    unsigned int cumemset_host_flags;
    unsigned int host : 1;
    unsigned int device : 1;
} LineFlags_t;


int free_line_params_host(LineParams_t ** const line_params,
                          LineFlags_t const flags);


int parse_hitran_file(LineParams_t ** const line_params,
                      char const * const filename,
                      LineFlags_t const flags,
                      int const mol_id,
                      double const w0,
                      double const wn);


int alloc_line_params_device(LineParams_t ** const line_params,
                             unsigned int const num_lines);


int free_line_params_device(LineParams_t ** const line_params);


#endif
