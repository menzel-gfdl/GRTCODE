#include "extern.h"
#include "verbosity.h"


static int verbosity = RS_NONE;


EXTERN void rs_set_verbosity(int const level)
{
    verbosity = level;
}


EXTERN int rs_get_verbosity()
{
    return verbosity;
}
