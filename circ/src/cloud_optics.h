#ifndef CLOUD_OPTICS_H_
#define CLOUD_OPTICS_H_

#include "extern.h"
#include "floating_point_type.h"
#include "optics.h"


EXTERN int cloud_optics(Optics_t * const optics,
                        fp_t * const liquid_water_path,
                        fp_t * const droplet_equivalent_radius
                       );


#endif
