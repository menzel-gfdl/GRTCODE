#include "verbosity.h"


static int verbosity = NONE;


void set_verbosity(int const level)
{
    verbosity = level;
}


int get_verbosity()
{
    return verbosity;
}
