#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "molecules.h"


static int const MIN_NAME_LEN = 4;


int get_mol_name(int const id,
                 char * const name,
                 int const name_len)
{
    if (name_len < MIN_NAME_LEN)
    {
        fatal(VALUE_ERR,
              "input name buffer must be at least %d characters long.",
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
            fatal(VALUE_ERR,
                  "unrecognized molecule id %d.",
                  id);
    }
    return SUCCESS;
}


int HITRAN_id_to_model_id(int const hitran_id,
                          int * const model_id)
{
    not_null(model_id);
    switch (hitran_id)
    {
        case Hitran_H2O:
            *model_id = H2O;
            break;
        case Hitran_CO2:
            *model_id = CO2;
            break;
        case Hitran_O3:
            *model_id = O3;
            break;
        case Hitran_N2O:
            *model_id = N2O;
            break;
        case Hitran_CO:
            *model_id = CO;
            break;
        case Hitran_CH4:
            *model_id = CH4;
            break;
        case Hitran_O2:
            *model_id = O2;
            break;
        default:
            fatal(VALUE_ERR,
                  "unrecognized HITRAN molecule id %d.",
                  hitran_id);
    }
    return SUCCESS;
}


int molecule_hash(int const mol_id,
                  int * const hash)
{
    not_null(hash);
    *hash = 0;
    unsigned int a = mol_id;
    if (mol_id < 1)
    {
        fatal(VALUE_ERR,
              "input molecule id (%d) must be at least 1.",
              mol_id);
    }
    while (a != 1)
    {
        a = a >> 1;
        (*hash)++;
    }
    if (*hash >= NUM_MOLS)
    {
        fatal(VALUE_ERR,
              "hash (%d) of molecule id (%d) cannot be >= %d.",
              *hash,
              mol_id,
              NUM_MOLS);
    }
    return SUCCESS;
}


int is_molecule_active(int const molecule_bit_field,
                       int const mol_id)
{
    return molecule_bit_field & mol_id;
}
