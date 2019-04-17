#ifndef VERBOSITY_H_
#define VERBOSITY_H_


enum verbosity
{
    NONE,
    INFO,
    WARN,
    ERROR
};


void set_verbosity(int const level);


int get_verbosity();


#endif
