#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "molecules.h"

static int const MIN_NAME_LEN = 4;

int get_mol_name(int const id,
                 char *name,
                 int const name_len)
{
    if (name_len < MIN_NAME_LEN)
    {
        fatal("input name buffer must be at least %d characters long.",
              MIN_NAME_LEN);
    }
    switch (id)
    {
        case H2O:
            snprintf(name,
                     name_len,
                     "h2o");
            break;
        case CO2:
            snprintf(name,
                     name_len,
                     "co2");
            break;
        case O3:
            snprintf(name,
                     name_len,
                     "o3");
            break;
        case N2O:
            snprintf(name,
                     name_len,
                     "n2o");
            break;
        case CO:
            snprintf(name,
                     name_len,
                     "co");
            break;
        case CH4:
            snprintf(name,
                     name_len,
                     "ch4");
            break;
        case O2:
            snprintf(name,
                     name_len,
                     "o2");
            break;
        default:
            fatal("unrecognized molecule id %d.",
                  id);
    }
    return SUCCESS;
}

int HITRAN_id_to_model_id(int hitran_id,
                          int *model_id)
{
    not_null(model_id);
    *model_id = hitran_id - 1;
    if (*model_id < 0 || *model_id >= NUM_MOL)
    {
        fatal("the model id (%d) calculated from the HITRAN id (%d) must be"
                  " >= 0 and <= %d.",
              *model_id,
              hitran_id,
              NUM_MOL);
    }
    return SUCCESS;
}
