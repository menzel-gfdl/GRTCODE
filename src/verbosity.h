#ifndef VERBOSITY_H_
#define VERBOSITY_H_


enum verbosity
{
    ERROR,
    WARN,
    INFO,
    NONE,
};


void set_verbosity(int const level);


int get_verbosity();


#endif
