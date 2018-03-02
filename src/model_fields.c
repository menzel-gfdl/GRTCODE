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
    malloc_ptr(fields->TSURF,n);
    malloc_ptr(fields->EMIS,n);
    malloc_ptr(fields->x,nmol);
    n *= nlev;
    malloc_ptr(fields->P,n);
    malloc_ptr(fields->T,n);
    int i;
    for (i=0;i<nmol;++i)
    {
        malloc_ptr(fields->x[i],n);
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
    int i;
    for (i=0;i<fields->nmol;++i)
    {
        free(fields->x[i]);
    }
    free(fields->x);
    return SUCCESS;
}
