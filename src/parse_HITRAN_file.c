#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "debug.h"
#include "floating_point_type.h"
#include "molecules.h"
#include "parse_HITRAN_file.h"
#include "utils.h"

#ifdef __NVCC__
#include "cuda_helpers.cuh"
#endif


typedef enum RefLinePtrIdx
{
    mol_pidx,
    iso_pidx,
    Vnn_pidx,
    Snn_ref_pidx,
    Yair_pidx,
    Yself_pidx,
    En_pidx,
    n_pidx,
    d_pidx
} RefLinePtrIdx_t;


typedef enum HITRAN2012_cols
{
    mol_c,
    iso_c,
    Vnn_c,
    Snn_c,
    A_c,
    Yair_c,
    Yself_c,
    Elo_c,
    n_c,
    del_c,
    Vu_c,
    Vl_c,
    Qu_c,
    Ql_c,
    Ierr_c,
    Iref_c,
    flag_c,
    gu_c,
    gl_c,
    NCOLS
} HITRAN2012_col_t;


typedef enum LookupCast
{
    NIL,
    I32,
    F32,
    F64
} LookupCast_t;


typedef union HITRAN2012_vals
{
    void *nil;
    int i;
    float f;
    double d;
} HITRAN2012_vals_t;


static unsigned int const HITRAN2012_recordLen = 160;
static unsigned int const HITRAN2012_pad = 2;
static unsigned int const HITRAN2012_fmt[NCOLS][2] =
{
    {2 , I32},
    {1 , I32},
    {12, F64},
    {10, F64},
    {10, NIL},
    {5 , F32},
    {5,  F32},
    {10, F32},
    {4,  F32},
    {8,  F32},
    {15, NIL},
    {15, NIL},
    {15, NIL},
    {15, NIL},
    {6,  NIL},
    {12, NIL},
    {1,  NIL},
    {7,  NIL},
    {7,  NIL}
};


static int alloc_line_params_host(LineParams_t ** const line_params,
                                  unsigned int const num_lines,
                                  LineFlags_t const flags)
{
    not_null(line_params);
    unsigned int const cuflags = flags.cumemset_host_flags;
    LineParams_t *self = NULL;
    check(malloc_ptr((void **)(&(self)),
                     sizeof(*self)));
    self->num_lines = num_lines;
    self->mol = -1;
    if (cuflags != ((unsigned int)-1))
    {
#ifdef __NVCC__
        HANDLE_ERROR(cudaHostAlloc(&(self->iso),
                                   num_lines*sizeof(*(self->iso)),
                                   cuflags));
        HANDLE_ERROR(cudaHostAlloc(&(self->vnn),
                                   num_lines*sizeof(*(self->vnn)),
                                   cuflags));
        HANDLE_ERROR(cudaHostAlloc(&(self->snn_ref),
                                   num_lines*sizeof(*(self->snn_ref)),
                                   cuflags));
        HANDLE_ERROR(cudaHostAlloc(&(self->yair),
                                   num_lines*sizeof(*(self->yair)),
                                   cuflags));
        HANDLE_ERROR(cudaHostAlloc(&(self->yself),
                                   num_lines*sizeof(*(self->yself)),
                                   cuflags));
        HANDLE_ERROR(cudaHostAlloc(&(self->en),
                                   num_lines*sizeof(*(self->en)),
                                   cuflags));
        HANDLE_ERROR(cudaHostAlloc(&(self->n),
                                   num_lines*sizeof(*(self->n)),
                                   cuflags));
        HANDLE_ERROR(cudaHostAlloc(&(self->d),
                                   num_lines*sizeof(*(self->d)),
                                   cuflags));
#endif
    }
    else
    {
        check(malloc_ptr((void **)(&(self->iso)),
                         sizeof(*(self->iso))*num_lines));
        check(malloc_ptr((void **)(&(self->vnn)),
                         sizeof(*(self->vnn))*num_lines));
        check(malloc_ptr((void **)(&(self->snn_ref)),
                         sizeof(*(self->snn_ref))*num_lines));
        check(malloc_ptr((void **)(&(self->yair)),
                         sizeof(*(self->yair))*num_lines));
        check(malloc_ptr((void **)(&(self->yself)),
                         sizeof(*(self->yself))*num_lines));
        check(malloc_ptr((void **)(&(self->en)),
                         sizeof(*(self->en))*num_lines));
        check(malloc_ptr((void **)(&(self->n)),
                         sizeof(*(self->n))*num_lines));
        check(malloc_ptr((void **)(&(self->d)),
                         sizeof(*(self->n))*num_lines));
    }
    *line_params = self;
    return SUCCESS;
}


int free_line_params_host(LineParams_t ** const line_params,
                          LineFlags_t const flags)
{
    not_null(line_params);
    not_null(*line_params);
    LineParams_t *self = *line_params;
    unsigned int const cuflags = flags.cumemset_host_flags;
    if (cuflags != ((unsigned int)-1))
    {
#ifdef __NVCC__
        HANDLE_ERROR(cudaFreeHost(self->iso));
        HANDLE_ERROR(cudaFreeHost(self->vnn));
        HANDLE_ERROR(cudaFreeHost(self->snn_ref));
        HANDLE_ERROR(cudaFreeHost(self->yair));
        HANDLE_ERROR(cudaFreeHost(self->yself));
        HANDLE_ERROR(cudaFreeHost(self->en));
        HANDLE_ERROR(cudaFreeHost(self->n));
        HANDLE_ERROR(cudaFreeHost(self->d));
#endif
    }
    else
    {
        free(self->iso);
        free(self->vnn);
        free(self->snn_ref);
        free(self->yair);
        free(self->yself);
        free(self->en);
        free(self->n);
        free(self->d);
    }
    free(self);
    *line_params = NULL;
    return SUCCESS;
}


static int realloc_line_params_host(LineParams_t ** const line_params,
                                    LineFlags_t const old_flags,
                                    LineFlags_t const new_flags)
{
    not_null(line_params);
    not_null(*line_params);
    LineParams_t *old_ptr = *line_params;
    LineParams_t *new_ptr = NULL;
    check(alloc_line_params_host(&new_ptr,
                                 old_ptr->num_lines,
                                 new_flags));
    new_ptr->mol = old_ptr->mol;
    new_ptr->num_lines = old_ptr->num_lines;
    unsigned int i;
    for (i=0;i<old_ptr->num_lines;++i)
    {
        new_ptr->iso[i] = old_ptr->iso[i];
        new_ptr->vnn[i] = old_ptr->vnn[i];
        new_ptr->snn_ref[i] = old_ptr->snn_ref[i];
        new_ptr->yair[i] = old_ptr->yair[i];
        new_ptr->yself[i] = old_ptr->yself[i];
        new_ptr->en[i] = old_ptr->en[i];
        new_ptr->n[i] = old_ptr->n[i];
        new_ptr->d[i] = old_ptr->d[i];
    }
    check(free_line_params_host(&old_ptr,
                                old_flags));
    *line_params = new_ptr;
    return SUCCESS;
}


static int HITRAN2012_cast(HITRAN2012_vals_t * const val,
                           int const col,
                           char const * const sval)
{
    not_null(val);
    not_null(sval);
    LookupCast_t const typ = (LookupCast_t)(HITRAN2012_fmt[col][1]);
    switch (typ)
    {
        case NIL:
            val->nil = NULL;
            break;
        case I32:
            check(to_int(sval,
                         &(val->i)));
            break;
        case F64:
        case F32:
            check(to_double(sval,
                            &(val->d)));
            if (typ == F32)
            {
                if (val->d >= -1.f*FLT_MAX && val->d <= FLT_MAX)
                {
                    val->f = val->d;
                }
                else
                {
                    fatal(VALUE_ERR,
                          "value %e from column %d cannot be safely"
                              " cast as a float.",
                          val->d,
                          col);
                }
            }
            break;
        default:
            fatal(VALUE_ERR,
                  "cast failed on col %d, LookupCast_t %d, sval: %s.",
                  HITRAN2012_fmt[col][0],
                  HITRAN2012_fmt[col][1],
                  sval);
    }
    return SUCCESS;
}


int parse_hitran_file(LineParams_t ** const line_params,
                      char const * const filename,
                      LineFlags_t const flags,
                      int const mol_id,
                      double const w0,
                      double const wn)
{
    not_null(line_params);
    not_null(filename);

    /*Open the file.*/
    log_info("Opening and reading HITRAN line parameters from file %s.",
             filename);
    FILE *fp = NULL;
    open_file(fp,
              filename,
              "r");

    /*Count the number of lines in the file.*/
    size_t const max_line = 163;
    char* buf = (char *)calloc(max_line,sizeof(*buf));
    not_null(buf);
    ssize_t ll = 0;
    size_t l = 0;
    unsigned int n = 0;
    while ((ll=getline(&buf,&l,fp)) != -1)
    {
        ++n;
    }
    rewind(fp);

    /*Malloc space.*/
    LineParams_t *lines = NULL;
    check(alloc_line_params_host(&lines,
                                 n,
                                 flags));
    lines->mol = mol_id;

    /*Parse out the line parameters.*/
    n = 0;
    while((ll=getline(&buf,&l,fp)) != -1)
    {
        if ((ll-HITRAN2012_recordLen) > HITRAN2012_pad)
        {
            fatal(VALUE_ERR,
                  "Found bad record at line %u (%zu exceeds max %zu chars)"
                      " in file %s.",
                  n,
                  (size_t)ll,
                  max_line,
                  filename);
        }
        else if (ll < HITRAN2012_recordLen)
        {
            fatal(VALUE_ERR,
                  "Found bad record at line %u (%zu less than %zu chars)"
                      " in file %s.",
                  n,
                  (size_t)ll,
                  max_line,
                  filename);
        }
        unsigned int val_idx = 0;
        size_t offset = 0;
        int col;
        for (col=0;col<NCOLS;++col)
        {
            size_t len = HITRAN2012_fmt[col][0];
            unsigned int t = HITRAN2012_fmt[col][1];
            char tmp[16];
            strncpy(tmp,&(buf[offset]),len);
            tmp[len]='\0';
            offset += len;
            if (t != NIL)
            {
                HITRAN2012_vals_t val;
                check(HITRAN2012_cast(&val,
                                      col,
                                      tmp));
                int m;
                switch (val_idx)
                {
                    case mol_pidx:
                        check(HITRAN_id_to_model_id(val.i,
                                                    &m));
                        if (m != lines->mol)
                        {
                            continue;
                        }
                        break;
                    case iso_pidx:
                        lines->iso[n] = val.i;
                        break;
                    case Vnn_pidx:
                        lines->vnn[n] = (fp_t)(val.d);
                        break;
                    case Snn_ref_pidx:
                        lines->snn_ref[n] = (fp_t)(val.d);
                        break;
                    case Yair_pidx:
                        lines->yair[n] = val.f;
                        break;
                    case Yself_pidx:
                        lines->yself[n] = val.f;
                        break;
                    case En_pidx:
                        lines->en[n] = val.f;
                        break;
                    case n_pidx:
                        lines->n[n] = val.f;
                        break;
                    case d_pidx:
                        lines->d[n] = val.f;
                        break;
                    default:
                        fatal(VALUE_ERR,
                              "Unknown column index (%d) on line %u in file"
                                  "%s.",
                              val_idx,
                              n,
                              filename);
                }
                ++val_idx;
            }
        }
        if ((w0 < 0) && (wn < 0))
        {
            /*Include the line for the calculation.*/
            ++n;
        }
        else if ((lines->vnn[n] >= w0) && (lines->vnn[n] <= wn))
        {
            /*Include the line for the calculation.*/
            ++n;
        }
    }
    free(buf);

    /*Set the number of lines that will be used in the calculation.*/
    lines->num_lines = n;

    /*Close the file.*/
    if (fclose(fp))
    {
        fatal(IO_ERR,
              "error closing file %s.",
              filename);
    }

    if ((w0 >= 0) || (wn >= 0))
    {
        check(realloc_line_params_host(&lines,
                                       flags,
                                       flags));
    }
    *line_params = lines;
    return SUCCESS;
}


int alloc_line_params_device(LineParams_t ** const line_params,
                             unsigned int const num_lines)
{
    not_null(line_params);
    LineParams_t *self = NULL;
    check(malloc_ptr((void **)(&self),
                     sizeof(*self)));
    self->num_lines = num_lines;
    self->mol = -1;
#ifdef __NVCC__
    HANDLE_ERROR(cudaMalloc(&(self->iso),
                            num_lines*sizeof(*(self->iso))));
    HANDLE_ERROR(cudaMalloc(&(self->vnn),
                            num_lines*sizeof(*(self->vnn))));
    HANDLE_ERROR(cudaMalloc(&(self->snn_ref),
                            num_lines*sizeof(*(self->snn_ref))));
    HANDLE_ERROR(cudaMalloc(&(self->yair),
                            num_lines*sizeof(*(self->yair))));
    HANDLE_ERROR(cudaMalloc(&(self->yself),
                            num_lines*sizeof(*(self->yself))));
    HANDLE_ERROR(cudaMalloc(&(self->en),
                            num_lines*sizeof(*(self->en))));
    HANDLE_ERROR(cudaMalloc(&(self->n),
                            num_lines*sizeof(*(self->n))));
    HANDLE_ERROR(cudaMalloc(&(self->d),
                            num_lines*sizeof(*(self->d))));
#endif
    *line_params = self;
    return SUCCESS;
}


int free_line_params_device(LineParams_t ** const line_params)
{
    not_null(line_params);
    not_null(*line_params);
    LineParams_t *self = *line_params;
#ifdef __NVCC__
    HANDLE_ERROR(cudaFree(self->iso));
    HANDLE_ERROR(cudaFree(self->vnn));
    HANDLE_ERROR(cudaFree(self->snn_ref));
    HANDLE_ERROR(cudaFree(self->yair));
    HANDLE_ERROR(cudaFree(self->yself));
    HANDLE_ERROR(cudaFree(self->en));
    HANDLE_ERROR(cudaFree(self->n));
    HANDLE_ERROR(cudaFree(self->d));
#endif
    free(self);
    *line_params = NULL;
    return SUCCESS;
}
