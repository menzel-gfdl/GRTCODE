#include <stdlib.h>
#include "debug.h"
#include "model_fields.h"
#include "utils.h"


int init_req_model_fields(req_model_fields_t * const fields,
                          int const nt,
                          int const nlon,
                          int const nlat,
                          int const nlev,
                          int const nmol)
{
    not_null(fields);
    fields->ntime = nt;
    fields->nlon = nlon;
    fields->nlat = nlat;
    fields->nlevel = nlev;
    fields->nmol = nmol;
    int n = nt*nlon*nlat;
    check(malloc_ptr((void **)(&(fields->TSURF)),
                     sizeof(*(fields->TSURF))*n));
    check(malloc_ptr((void **)(&(fields->EMIS)),
                     sizeof(*(fields->EMIS))*n));
    check(malloc_ptr((void **)(&(fields->SFC_DIR_ALB)),
                     sizeof(*(fields->SFC_DIR_ALB))*n));
    check(malloc_ptr((void **)(&(fields->SFC_DIF_ALB)),
                     sizeof(*(fields->SFC_DIF_ALB))*n));
    check(malloc_ptr((void **)(&(fields->COS_SOL_ZEN_ANG)),
                     sizeof(*(fields->COS_SOL_ZEN_ANG))*n));
    check(malloc_ptr((void **)(&(fields->TOTAL_SOL_FLUX)),
                     sizeof(*(fields->TOTAL_SOL_FLUX))*n));
    check(malloc_ptr((void **)(&(fields->x)),
                     sizeof(*(fields->x))*nmol));
    n *= nlev;
    check(malloc_ptr((void **)(&(fields->P)),
                     sizeof(*(fields->P))*n));
    check(malloc_ptr((void **)(&(fields->T)),
                     sizeof(*(fields->T))*n));
    int i;
    for (i=0;i<nmol;++i)
    {
        check(malloc_ptr((void **)(&(fields->x[i])),
                         sizeof(*(fields->x[i]))*n));
    }
    return SUCCESS;
}

int free_req_model_fields(req_model_fields_t * const fields)
{
    not_null(fields);
    free(fields->P);
    free(fields->T);
    free(fields->TSURF);
    free(fields->EMIS);
    free(fields->SFC_DIR_ALB);
    free(fields->SFC_DIF_ALB);
    free(fields->COS_SOL_ZEN_ANG);
    int i;
    for (i=0;i<fields->nmol;++i)
    {
        free(fields->x[i]);
    }
    free(fields->x);
    return SUCCESS;
}
