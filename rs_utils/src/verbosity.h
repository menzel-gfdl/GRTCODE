#ifndef VERBOSITY_H_
#define VERBOSITY_H_


enum verbosity
{
    RS_NONE,
    RS_ERROR,
    RS_WARN,
    RS_INFO
};


void rs_set_verbosity(int const level);


int rs_get_verbosity();


#endif
