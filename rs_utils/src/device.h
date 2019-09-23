#ifndef DEVICE_H_
#define DEVICE_H_

#include "extern.h"


/** @brief Device object.*/
typedef int Device_t;


/** @brief Determine the number of CUDA-enabled GPUs on the system.
    @return RS_SUCCESS or an error code.*/
EXTERN int get_num_gpus(int * num_devices, /**< Number of CUDA-enabled devices found.*/
                        int const verbose /**< Verbosity flag.*/
                       );


/** @brief Set the device identifier.
    @return RS_SUCCESS or an error code.*/
EXTERN int create_device(Device_t * const device, /**< Device object.*/
                         int const * const id /**< Device identifier.*/
                        );


#endif
