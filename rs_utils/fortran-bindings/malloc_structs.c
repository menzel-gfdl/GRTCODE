#include <stdlib.h>
#include <string.h>
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


/** @brief Retrieve arrays from optics structure.*/
/** @return RS_SUCCESS or an error code.*/
int optical_properties(Optics_t const * const optics, fp_t * const tau, fp_t * const omega,
                       fp_t * const g)
{
    not_null(optics);
    size_t n = (optics->num_layers)*(optics->grid.n);
    if (tau != NULL)
    {
        memcpy(tau, optics->tau, sizeof(*(optics->tau))*n);
    }
    if (omega != NULL)
    {
        memcpy(omega, optics->omega, sizeof(*(optics->omega))*n);
    }
    if (g != NULL)
    {
        memcpy(g, optics->g, sizeof(*(optics->g))*n);
    }
    return RS_SUCCESS;
}


/** @brief Get the spectral grid properties.
    @return RS_SUCCESS or an error code.*/
int spectral_grid_properties(SpectralGrid_t const * const grid, /**< Spectral grid.*/
                             double * const w0, /**< Grid lower bound.*/
                             uint64_t * const n, /**< Spectral grid size.*/
                             double * const dw /**< Grid spacing.*/
                            )
{
    not_null(grid);
    if (w0 != NULL)
    {
        *w0 = grid->w0;
    }
    if (n != NULL);
    {
        *n = grid->n;
    }
    if (dw != NULL)
    {
        *dw = grid->dw;
    }
    return RS_SUCCESS;
}
