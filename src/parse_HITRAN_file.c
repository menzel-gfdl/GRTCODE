#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "debug.h"
#include "floating_point_type.h"
#include "parse_HITRAN_file.h"
#include "tips2017.h"
#include "utils.h"


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


static int alloc_line_params(LineParams_t * const line_params,
                             uint64_t const num_lines)
{
    not_null(line_params);
    line_params->num_lines = num_lines;
    check(malloc_ptr((void **)(&(line_params->iso)),
                     sizeof(*(line_params->iso))*num_lines));
    check(malloc_ptr((void **)(&(line_params->vnn)),
                     sizeof(*(line_params->vnn))*num_lines));
    check(malloc_ptr((void **)(&(line_params->snn)),
                     sizeof(*(line_params->snn))*num_lines));
    check(malloc_ptr((void **)(&(line_params->yair)),
                     sizeof(*(line_params->yair))*num_lines));
    check(malloc_ptr((void **)(&(line_params->yself)),
                     sizeof(*(line_params->yself))*num_lines));
    check(malloc_ptr((void **)(&(line_params->en)),
                     sizeof(*(line_params->en))*num_lines));
    check(malloc_ptr((void **)(&(line_params->n)),
                     sizeof(*(line_params->n))*num_lines));
    check(malloc_ptr((void **)(&(line_params->d)),
                     sizeof(*(line_params->n))*num_lines));
    return SUCCESS;
}


int free_line_params(LineParams_t * const line_params)
{
    not_null(line_params);
    check(free_ptr((void **)&(line_params->iso)));
    check(free_ptr((void **)&(line_params->vnn)));
    check(free_ptr((void **)&(line_params->snn)));
    check(free_ptr((void **)&(line_params->yair)));
    check(free_ptr((void **)&(line_params->yself)));
    check(free_ptr((void **)&(line_params->en)));
    check(free_ptr((void **)&(line_params->n)));
    check(free_ptr((void **)&(line_params->d)));
    return SUCCESS;
}


static int realloc_line_params(LineParams_t * const line_params)
{
    not_null(line_params);
    LineParams_t t;
    t.iso = line_params->iso;
    t.vnn = line_params->vnn;
    t.snn = line_params->snn;
    t.yair = line_params->yair;
    t.yself = line_params->yself;
    t.en = line_params->en;
    t.n = line_params->n;
    t.d = line_params->d;
    check(alloc_line_params(line_params,
                            line_params->num_lines));
    memcpy(line_params->iso,
           t.iso,
           sizeof(*(t.iso))*line_params->num_lines);
    memcpy(line_params->vnn,
           t.vnn,
           sizeof(*(t.vnn))*line_params->num_lines);
    memcpy(line_params->snn,
           t.snn,
           sizeof(*(t.snn))*line_params->num_lines);
    memcpy(line_params->yair,
           t.yair,
           sizeof(*(t.yair))*line_params->num_lines);
    memcpy(line_params->yself,
           t.yself,
           sizeof(*(t.yself))*line_params->num_lines);
    memcpy(line_params->en,
           t.en,
           sizeof(*(t.en))*line_params->num_lines);
    memcpy(line_params->n,
           t.n,
           sizeof(*(t.n))*line_params->num_lines);
    memcpy(line_params->d,
           t.d,
           sizeof(*(t.d))*line_params->num_lines);
    check(free_line_params(&t));
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


int parse_hitran_file(LineParams_t * const line_params,
                      char const * const filename,
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
    char* buf;
    check(malloc_ptr((void **) &buf,
                     sizeof(*buf)*max_line));
    memset(buf,
           0,
           sizeof(*buf)*max_line);
    ssize_t ll = 0;
    size_t l = 0;
    uint64_t n = 0;
    while ((ll=getline(&buf,&l,fp)) != -1)
    {
        ++n;
    }
    rewind(fp);

    /*Malloc space.*/
    check(alloc_line_params(line_params,
                            n));

    /*Parse out the line parameters.*/
    n = 0;
    size_t line_count = 0;
    while((ll=getline(&buf,&l,fp)) != -1)
    {
        line_count++;
        if ((ll-HITRAN2012_recordLen) > HITRAN2012_pad)
        {
            fatal(VALUE_ERR,
                  "Found bad record at line %zu (%zu exceeds max %zu chars)"
                      " in file %s.",
                  line_count,
                  (size_t)ll,
                  max_line,
                  filename);
        }
        else if (ll < HITRAN2012_recordLen)
        {
            fatal(VALUE_ERR,
                  "Found bad record at line %zu (%zu less than %zu chars)"
                      " in file %s.",
                  line_count,
                  (size_t)ll,
                  max_line,
                  filename);
        }
        unsigned int val_idx = 0;
        size_t offset = 0;
        int col;
        int go_to_next_line = 0;
        for (col=0;col<NCOLS;++col)
        {
            if (go_to_next_line)
            {
                break;
            }
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

                /*Column indices for parsing.*/
                enum RefLinePtrIdx
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
                };

                switch (val_idx)
                {
                    case mol_pidx:
                        if (val.i != mol_id)
                        {
                            go_to_next_line = 1;
                            continue;
                        }
                        break;
                    case iso_pidx:
                        line_params->iso[n] = val.i;
                        break;
                    case Vnn_pidx:
                        line_params->vnn[n] = (fp_t)(val.d);
                        break;
                    case Snn_ref_pidx:
                        line_params->snn[n] = (fp_t)(val.d);
                        break;
                    case Yair_pidx:
                        line_params->yair[n] = val.f;
                        break;
                    case Yself_pidx:
                        line_params->yself[n] = val.f;
                        break;
                    case En_pidx:
                        line_params->en[n] = val.f;
                        break;
                    case n_pidx:
                        line_params->n[n] = val.f;
                        break;
                    case d_pidx:
                        line_params->d[n] = val.f;
                        break;
                    default:
                        fatal(VALUE_ERR,
                              "Unknown column index (%d) on line %zu in file"
                                  "%s.",
                              val_idx,
                              n,
                              filename);
                }
                ++val_idx;
            }
        }
        if (!go_to_next_line)
        {
            if ((w0 < 0 && wn < 0) || (line_params->vnn[n] >= w0 &&
                line_params->vnn[n] <= wn))
            {
                /*Include the line for the calculation.*/
                ++n;
            }
        }
    }
    free(buf);

    /*Close the file.*/
    if (fclose(fp))
    {
        fatal(IO_ERR,
              "error closing file %s.",
              filename);
    }

    /*Reallocate if necessary.*/
    if (line_params->num_lines != n)
    {
        line_params->num_lines = n;
        check(realloc_line_params(line_params));
    }

    /*Adjust the raw read-in line strengths.*/
    fp_t const tref = 296.f;
    fp_t const c2 = -1.4387686f;
    fp_t *snn = line_params->snn;
    int *iso = line_params->iso;
    fp_t *en = line_params->en;
    fp_t *vnn = line_params->vnn;
    uint64_t i;
#pragma omp parallel for default(none) private(i) shared(snn,iso,en,vnn,n)
    for (i=0;i<n;++i)
    {
        snn[i] *= Q(mol_id,tref,iso[i])/(EXP(c2*en[i]/tref)*
                  (1.f - EXP(c2*vnn[i]/tref)));
    }
    return SUCCESS;
}
