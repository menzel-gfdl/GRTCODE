#include "debug.h"
#include "device.h"
#include "extern.h"
#include "rs_config.h"


/** @brief Determine the number of CUDA-enabled GPUs on the system.
    @return RS_SUCCESS or an error code.*/
static int get_num_gpus(int * num_devices, /**< Number of CUDA-enabled devices found.*/
                        int const verbose /**< Verbosity flag.*/
                       )
{
    not_null(num_devices);
#ifdef __NVCC__
    gpu_catch(cudaGetDeviceCount(num_devices));
    if (verbose)
    {
        char const *mesg = "Found %d GPU devices:";
        log_mesg(mesg, *num_devices);
        int i;
        for (i=0; i<(*num_devices); ++i)
        {
            cudaDeviceProp prop;
            gpu_catch(cudaGetDeviceProperties(&prop, i));
            mesg = "\tDevice #%d: %s";
            log_mesg(mesg, i, prop.name);
        }
    }
#else
    *num_devices = 0;
#endif
    return RS_SUCCESS;
}


/*Set the device identifier.*/
EXTERN int create_device(Device_t * const device, int const * const id)
{
    not_null(device);
    int num_devices;
    catch(get_num_gpus(&num_devices, 1));
    if (id != NULL)
    {
        if (*id != HOST_ONLY)
        {
            in_range(*id, 0, num_devices);
        }
        *device = *id;
    }
    else if (num_devices > 0)
    {
        *device = DEFAULT_GPU;
    }
    else
    {
        *device = HOST_ONLY;
    }
    return RS_SUCCESS;
}
