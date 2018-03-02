#include "debug.h"
#include "output_fields.h"
#include "utils.h"


int alloc_output_fields(OutputFields_t * const var,
                        int const nws,
                        int const nlevels)
{
    not_null(var);
    malloc_ptr(var->tau,nlevels*nws);
    malloc_ptr(var->lw_flux_down,nlevels);
    malloc_ptr(var->lw_flux_up,nlevels);
    return SUCCESS;
}


int free_output_fields(OutputFields_t * const var)
{
    not_null(var);
    free(var->tau);
    free(var->lw_flux_down);
    free(var->lw_flux_up);
    return SUCCESS;
}
