#ifndef PARSE_HITRAN_FILE_H_
#define PARSE_HITRAN_FILE_H_

#include <stdint.h>
#include "floating_point_type.h"


typedef struct LineParams
{
    int *iso;
    fp_t *vnn;
    fp_t *snn;
    fp_t *yair;
    fp_t *yself;
    fp_t *en;
    fp_t *n;
    fp_t *d;
    uint64_t num_lines;
} LineParams_t;


int free_line_params(LineParams_t * const line_params);


int parse_hitran_file(LineParams_t * const line_params,
                      char const * const filename,
                      int const mol_id,
                      double const w0,
                      double const wn);


#endif
