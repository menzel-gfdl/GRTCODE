#ifndef VERBOSITY_H_
#define VERBOSITY_H_

#include "extern.h"


enum verbosity
{
    RS_NONE,
    RS_ERROR,
    RS_WARN,
    RS_INFO
};


EXTERN void rs_set_verbosity(int const level);


EXTERN int rs_get_verbosity();


#endif
