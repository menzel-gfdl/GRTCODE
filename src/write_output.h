#ifndef WRITE_OUTPUT_H_
#define WRITE_OUTPUT_H_

#include "floating_point_type.h"


int init_output_file(char const * const filename,
                     int * const ncid,
                     int const nlons,
                     int const nlats,
                     int const nlevels,
                     int const nws,
                     int const output_spectra);


int close_output_file(int const ncid);


int write_data_column(int const ncid,
                      fp_t *lw_flux_down,
                      fp_t *lw_flux_up,
                      fp_t *sw_flux_down,
                      fp_t *sw_flux_up,
                      fp_t *tau_gas,
                      fp_t *tau_scatter,
                      int const time,
                      int const lon,
                      int const lat,
                      int const nlevels,
                      int const nws,
                      int const output_spectra);


#endif
