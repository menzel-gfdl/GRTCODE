#include "verbosity.h"


static int verbosity = RS_NONE;


void rs_set_verbosity(int const level)
{
    verbosity = level;
}


int rs_get_verbosity()
{
    return verbosity;
}
