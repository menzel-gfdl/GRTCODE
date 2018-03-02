#include "constants.h"
#include "floating_point_type.h"

fp_t const PA_TO_ATM = 9.86923e-6;
fp_t const FROM_PPMV = 1.e-6;
int const MISSING_CONC = -1;
int const CONC_FROM_FILE = -2;
int const MIN_TIME = 0;
int const MIN_LON = 0;
int const MIN_LAT = 0;
int const MIN_WVN = 1;
int const MAX_WVN = 50000;
int const DEFAULT_WVN = 3000;
int const DEFAULT_DEVICE = 0;
double const RES_MIN = 1.e-4;
double const RES_MAX = 100.;
int const HOST_LAUNCH = 0;
int const DEVICE_LAUNCH = 1;
int const MAX_NUM_LINES = 524288; /*2^19.*/
