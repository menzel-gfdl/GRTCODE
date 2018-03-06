#include "debug.h"
#include "output_fields.h"
#include "utils.h"


int alloc_output_fields(OutputFields_t * const var,
                        int const nws,
                        int const nlevels,
                        int const device_launch)
{
    not_null(var);
    malloc_fp_ptr(var->tau,nlevels*nws);
    malloc_fp_ptr(var->lw_flux_down,nlevels);
    malloc_fp_ptr(var->lw_flux_up,nlevels);
    if (device_launch)
    {
        malloc_fp_ptr(var->lw_flux_down_per_w,nws*nlevels);
        malloc_fp_ptr(var->lw_flux_up_per_w,nws*nws);
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
