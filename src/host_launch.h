#ifndef HOST_LAUNCH_H_
#define HOST_LAUNCH_H_

#include "floating_point_type.h"
#include "molecular_lines.h"


int launch_h(GrtContext_t * const context,
             fp_t const * const p,
             fp_t const * const t,
             fp_t * const tau);


#endif
