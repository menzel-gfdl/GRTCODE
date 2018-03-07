#include "debug.h"
#include "output_fields.h"
#include "utils.h"


int alloc_output_fields(OutputFields_t * const var,
                        int const nws,
                        int const nlevels,
                        int const device_launch)
{
    not_null(var);
    check(malloc_ptr((void **)(&(var->tau)),
                     sizeof(*(var->tau))*nlevels*nws));
    check(malloc_ptr((void **)(&(var->lw_flux_down)),
                     sizeof(*(var->lw_flux_down))*nlevels));
    check(malloc_ptr((void **)(&(var->lw_flux_up)),
                     sizeof(*(var->lw_flux_up))*nlevels));
    if (device_launch)
    {
        check(malloc_ptr((void **)(&(var->lw_flux_down_per_w)),
                         sizeof(*(var->lw_flux_down_per_w))*nws*nlevels));
        check(malloc_ptr((void **)(&(var->lw_flux_up_per_w)),
                         sizeof(*(var->lw_flux_up_per_w))*nws*nlevels));
    }
    return SUCCESS;
}


int free_output_fields(OutputFields_t * const var,
                       int const device_launch)
{
    not_null(var);
    free(var->tau);
    free(var->lw_flux_down);
    free(var->lw_flux_up);
    if (device_launch)
    {
        free(var->lw_flux_down_per_w);
        free(var->lw_flux_up_per_w);
    }
    return SUCCESS;
}
