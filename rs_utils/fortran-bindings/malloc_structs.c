#include <stdlib.h>
#include "debug.h"
#include "molecular_lines.h"
#include "optics.h"
#include "spectral_grid.h"


enum StructTypes_t {
    GRID = 0,
    OPTICS,
    MOLECULAR_LINES
};


/** @brief Malloc struct and associate input pointer.*/
/** @return RS_SUCCESS or an error code.*/
int malloc_struct(void **p, int type)
{
    size_t s;
    switch (type)
    {
        case GRID:
            s = sizeof(SpectralGrid_t);
            break;
        case OPTICS:
            s = sizeof(Optics_t);
            break;
        case MOLECULAR_LINES:
            s = sizeof(MolecularLines_t);
            break;
        default:
            {char *mesg = "unrecognized structure type %d.";
            raise(RS_VALUE_ERR, mesg, type);}
    }
    not_null(p);
    is_null(*p);
    *p = malloc(s);
    not_null(*p);
    return RS_SUCCESS;
}


/** @brief Free struct associated with input pointer.*/
/** @return RS_SUCCESS or an error code.*/
int free_struct(void **p)
{
    gfree(*p, HOST_ONLY);
    *p = NULL;
    return RS_SUCCESS;
}
