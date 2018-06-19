#include "debug.h"
#include "output_fields.h"
#include "utils.h"


int alloc_output_fields(OutputFields_t * const var,
                        int const nws,
                        int const nlevels,
                        int const device_launch)
{
    not_null(var);
    int nlayers = nlevels - 1;
    check(malloc_ptr((void **)(&(var->tau_gas)),
                     sizeof(*(var->tau_gas))*nlayers*nws));
    check(malloc_ptr((void **)(&(var->tau_scatter)),
                     sizeof(*(var->tau_scatter))*nlayers*nws));
    check(malloc_ptr((void **)(&(var->lw_flux_down)),
                     sizeof(*(var->lw_flux_down))*nlevels));
    check(malloc_ptr((void **)(&(var->lw_flux_up)),
                     sizeof(*(var->lw_flux_up))*nlevels));
    check(malloc_ptr((void **)(&(var->sw_flux_down)),
                     sizeof(*(var->sw_flux_down))*nlevels));
    check(malloc_ptr((void **)(&(var->sw_flux_up)),
                     sizeof(*(var->sw_flux_up))*nlevels));
    if (device_launch)
    {
        check(malloc_ptr((void **)(&(var->lw_flux_down_per_w)),
                         sizeof(*(var->lw_flux_down_per_w))*nws*nlevels));
        check(malloc_ptr((void **)(&(var->lw_flux_up_per_w)),
                         sizeof(*(var->lw_flux_up_per_w))*nws*nlevels));
        check(malloc_ptr((void **)(&(var->sw_flux_down_per_w)),
                         sizeof(*(var->sw_flux_down_per_w))*nws*nlevels));
        check(malloc_ptr((void **)(&(var->sw_flux_up_per_w)),
                         sizeof(*(var->sw_flux_up_per_w))*nws*nlevels));
    }
    return SUCCESS;
}


int free_output_fields(OutputFields_t * const var,
                       int const device_launch)
{
    not_null(var);
    free(var->tau_gas);
    free(var->tau_scatter);
    free(var->lw_flux_down);
    free(var->lw_flux_up);
    free(var->sw_flux_down);
    free(var->sw_flux_up);
    if (device_launch)
    {
        free(var->lw_flux_down_per_w);
        free(var->lw_flux_up_per_w);
        free(var->sw_flux_down_per_w);
        free(var->sw_flux_up_per_w);
    }
    return SUCCESS;
}
