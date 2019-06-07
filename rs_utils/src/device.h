#ifndef DEVICE_H_
#define DEVICE_H_


/** @brief Device object.*/
typedef int Device_t;


/** @brief Set the device identifier.
    @return RS_SUCCESS or an error code.*/
int create_device(Device_t * const device, /**< Device object.*/
                  int const * const id /**< Device identifier.*/
                 );

#endif
