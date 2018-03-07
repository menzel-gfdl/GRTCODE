#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "debug.h"
#include "floating_point_type.h"
#include "molecules.h"
#include "parseHITRANfile.h"
#include "utils.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
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

static int alloc_line_params_host(line_params_t **lineParams,
                                  unsigned int nLines,
                                  line_flags_t flags)
{
    not_null(lineParams);
    is_null(*lineParams);
    unsigned int const cuFlags = flags.cumemset_host_flags;
    line_params_t *self = NULL;
    check(malloc_ptr((void **)(&(self)),
                     sizeof(*self)));
    self->nLines = nLines;
    self->mol = -1;
    if (cuFlags != ((unsigned int)-1))
    {
        using_gpu();
#ifdef __NVCC__
        HANDLE_ERROR(cudaHostAlloc(&(self->iso),
                                   nLines*sizeof(*(self->iso)),
                                   cuFlags));
        HANDLE_ERROR(cudaHostAlloc(&(self->Vnn),
                                   nLines*sizeof(*(self->Vnn)),
                                   cuFlags));
        HANDLE_ERROR(cudaHostAlloc(&(self->Snn_ref),
                                   nLines*sizeof(*(self->Snn_ref)),
                                   cuFlags));
        HANDLE_ERROR(cudaHostAlloc(&(self->Yair),
                                   nLines*sizeof(*(self->Yair)),
                                   cuFlags));
        HANDLE_ERROR(cudaHostAlloc(&(self->Yself),
                                   nLines*sizeof(*(self->Yself)),
                                   cuFlags));
        HANDLE_ERROR(cudaHostAlloc(&(self->En),
                                   nLines*sizeof(*(self->En)),
                                   cuFlags));
        HANDLE_ERROR(cudaHostAlloc(&(self->n),
                                   nLines*sizeof(*(self->n)),
                                   cuFlags));
        HANDLE_ERROR(cudaHostAlloc(&(self->d),
                                   nLines*sizeof(*(self->d)),
                                   cuFlags));
#endif
    }
    else
    {
        check(malloc_ptr((void **)(&(self->iso)),
                         sizeof(*(self->iso))*nLines));
        check(malloc_ptr((void **)(&(self->Vnn)),
                         sizeof(*(self->Vnn))*nLines));
        check(malloc_ptr((void **)(&(self->Snn_ref)),
                         sizeof(*(self->Snn_ref))*nLines));
        check(malloc_ptr((void **)(&(self->Yair)),
                         sizeof(*(self->Yair))*nLines));
        check(malloc_ptr((void **)(&(self->Yself)),
                         sizeof(*(self->Yself))*nLines));
        check(malloc_ptr((void **)(&(self->En)),
                         sizeof(*(self->En))*nLines));
        check(malloc_ptr((void **)(&(self->n)),
                         sizeof(*(self->n))*nLines));
        check(malloc_ptr((void **)(&(self->d)),
                         sizeof(*(self->n))*nLines));
    }
    *lineParams = self;
    return SUCCESS;
}

int free_line_params_host(line_params_t **lineParams,
                          line_flags_t flags)
{
    not_null(lineParams);
    not_null(*lineParams);
    line_params_t *self = *lineParams;
    unsigned int const cuFlags = flags.cumemset_host_flags;
    if (cuFlags != ((unsigned int)-1))
    {
        using_gpu();
#ifdef __NVCC__
        HANDLE_ERROR(cudaFreeHost(self->iso));
        HANDLE_ERROR(cudaFreeHost(self->Vnn));
        HANDLE_ERROR(cudaFreeHost(self->Snn_ref));
        HANDLE_ERROR(cudaFreeHost(self->Yair));
        HANDLE_ERROR(cudaFreeHost(self->Yself));
        HANDLE_ERROR(cudaFreeHost(self->En));
        HANDLE_ERROR(cudaFreeHost(self->n));
        HANDLE_ERROR(cudaFreeHost(self->d));
#endif
    }
    else
    {
        free(self->iso);
        free(self->Vnn);
        free(self->Snn_ref);
        free(self->Yair);
        free(self->Yself);
        free(self->En);
        free(self->n);
        free(self->d);
    }
    free(self);
    *lineParams = NULL;
    return SUCCESS;
}

static int realloc_line_params_host(line_params_t **lineParams,
                                    line_flags_t old_flags,
                                    line_flags_t new_flags)
{
    not_null(lineParams);
    not_null(*lineParams);
    line_params_t *old_ptr = *lineParams;
    line_params_t *new_ptr = NULL;
    check(alloc_line_params_host(&new_ptr,
                                 old_ptr->nLines,
                                 new_flags));
    new_ptr->mol = old_ptr->mol;
    new_ptr->nLines = old_ptr->nLines;
    unsigned int i;
    for (i=0;i<old_ptr->nLines;++i)
    {
        new_ptr->iso[i] = old_ptr->iso[i];
        new_ptr->Vnn[i] = old_ptr->Vnn[i];
        new_ptr->Snn_ref[i] = old_ptr->Snn_ref[i];
        new_ptr->Yair[i] = old_ptr->Yair[i];
        new_ptr->Yself[i] = old_ptr->Yself[i];
        new_ptr->En[i] = old_ptr->En[i];
        new_ptr->n[i] = old_ptr->n[i];
        new_ptr->d[i] = old_ptr->d[i];
    }
    check(free_line_params_host(&old_ptr,
                                old_flags));
    *lineParams = new_ptr;
    return SUCCESS;
}

static int HITRAN2012_cast(HITRAN2012_vals_t *val,
                           int const col,
                           char *sval)
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
                    fatal("value %e from column %d cannot be safely"
                              " cast as a float.",
                          val->d,
                          col);
                }
            }
            break;
        default:
            fatal("cast failed on col %d, LookupCast_t %d, sval: %s.",
                  HITRAN2012_fmt[col][0],
                  HITRAN2012_fmt[col][1],
                  sval);
    }
    return SUCCESS;
}

int parse_hitran_file(line_params_t **lineParams,
                      char *fileName,
                      line_flags_t flags,
                      double loWn,
                      double hiWn)
{
    not_null(lineParams);
    is_null(*lineParams);
    not_null(fileName);

    /*Open the file.*/
    log_mesg("Opening and reading HITRAN line parameters from file %s.",
             fileName);
    FILE *fp = NULL;
    open_file(fp,
              fileName,
              "r");

    /*Count the number of lines in the file.*/
    size_t const maxLine = 163;
    char* buf = (char *)calloc(maxLine,sizeof(*buf));
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
    line_params_t *lines = NULL;
    check(alloc_line_params_host(&lines,
                                 n,
                                 flags));

    /*Parse out the line parameters.*/
    n = 0;
    while((ll=getline(&buf,&l,fp)) != -1)
    {
        if ((ll-HITRAN2012_recordLen) > HITRAN2012_pad)
        {
            fatal("Found bad record at line %u (%zu exceeds max %zu chars)"
                      " in file %s.",
                  n,
                  (size_t)ll,
                  maxLine,
                  fileName);
        }
        else if (ll < HITRAN2012_recordLen)
        {
            fatal("Found bad record at line %u (%zu less than %zu chars)"
                      " in file %s.",
                  n,
                  (size_t)ll,
                  maxLine,
                  fileName);
        }
        unsigned int valIdx = 0;
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
                int mol_id;
                switch (valIdx)
                {
                    case mol_pidx:
                        check(HITRAN_id_to_model_id(val.i,
                                                    &mol_id));
                        if (n == 0)
                        {
                            lines->mol = mol_id;
                        }
                        else if (lines->mol != mol_id)
                        {
                            fatal("molecule changed from %d to %d at line %u"
                                      " in file %s.",
                                  lines->mol,
                                  mol_id,
                                  n,
                                  fileName);
                        }
                        break;
                    case iso_pidx:
                        lines->iso[n] = val.i;
                        break;
                    case Vnn_pidx:
                        lines->Vnn[n] = (fp_t)(val.d);
                        break;
                    case Snn_ref_pidx:
                        lines->Snn_ref[n] = (fp_t)(val.d);
                        break;
                    case Yair_pidx:
                        lines->Yair[n] = val.f;
                        break;
                    case Yself_pidx:
                        lines->Yself[n] = val.f;
                        break;
                    case En_pidx:
                        lines->En[n] = val.f;
                        break;
                    case n_pidx:
                        lines->n[n] = val.f;
                        break;
                    case d_pidx:
                        lines->d[n] = val.f;
                        break;
                    default:
                        fatal("Unknown column index (%d) on line %u in file"
                                  "%s.",
                              valIdx,
                              n,
                              fileName);
                }
                ++valIdx;
            }
        }
        if ((loWn < 0) && (hiWn < 0))
        {
            /*Include the line for the calculation.*/
            ++n;
        }
        else if ((lines->Vnn[n] >= loWn) && (lines->Vnn[n] <= hiWn))
        {
            /*Include the line for the calculation.*/
            ++n;
        }
    }
    free(buf);

    /*Set the number of lines that will be used in the calculation.*/
    lines->nLines = n;

    /*Close the file.*/
    if (fclose(fp))
    {
        fatal("error closing file %s.",
              fileName);
    }

    if ((loWn >=0) || (hiWn>=0))
    {
        check(realloc_line_params_host(&lines,
                                       flags,
                                       flags));
    }
    *lineParams = lines;
    return SUCCESS;
}

#ifdef __NVCC__
int alloc_line_params_device(line_params_t **lineParams,
                             unsigned int nLines)
{
    not_null(lineParams);
    is_null(*lineParams);
    line_params_t *self = NULL;
    check(malloc_ptr((void **)(&self),
                     sizeof(*self)));
    self->nLines = nLines;
    self->mol = -1;
    HANDLE_ERROR(cudaMalloc(&(self->iso),
                            nLines*sizeof(*(self->iso))));
    HANDLE_ERROR(cudaMalloc(&(self->Vnn),
                            nLines*sizeof(*(self->Vnn))));
    HANDLE_ERROR(cudaMalloc(&(self->Snn_ref),
                            nLines*sizeof(*(self->Snn_ref))));
    HANDLE_ERROR(cudaMalloc(&(self->Yair),
                            nLines*sizeof(*(self->Yair))));
    HANDLE_ERROR(cudaMalloc(&(self->Yself),
                            nLines*sizeof(*(self->Yself))));
    HANDLE_ERROR(cudaMalloc(&(self->En),
                            nLines*sizeof(*(self->En))));
    HANDLE_ERROR(cudaMalloc(&(self->n),
                            nLines*sizeof(*(self->n))));
    HANDLE_ERROR(cudaMalloc(&(self->d),
                            nLines*sizeof(*(self->d))));
    *lineParams = self;
    return SUCCESS;
}

int free_line_params_device(line_params_t **lineParams)
{
    not_null(lineParams);
    not_null(*lineParams);
    line_params_t *self = *lineParams;
    HANDLE_ERROR(cudaFree(self->iso));
    HANDLE_ERROR(cudaFree(self->Vnn));
    HANDLE_ERROR(cudaFree(self->Snn_ref));
    HANDLE_ERROR(cudaFree(self->Yair));
    HANDLE_ERROR(cudaFree(self->Yself));
    HANDLE_ERROR(cudaFree(self->En));
    HANDLE_ERROR(cudaFree(self->n));
    HANDLE_ERROR(cudaFree(self->d));
    free(self);
    *lineParams = NULL;
    return SUCCESS;
}
#endif
