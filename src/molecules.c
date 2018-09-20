#include <stdint.h>
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


int molecule_hash(int const mol_id,
                  int * const hash)
{
    not_null(hash);
    if (mol_id < H2O || mol_id > COCl2)
    {
        fatal(VALUE_ERR,
              "unrecognized molecule id %d.",
              mol_id);
    }
    *hash = mol_id - 1;
    return SUCCESS;
}


int activate_molecule(uint64_t * const molecule_bit_field,
                      int const mol_id)
{
    not_null(molecule_bit_field);
    uint64_t const one = 1;
    *molecule_bit_field = (*molecule_bit_field) | (one << (mol_id-1));
    return SUCCESS;
}


int is_molecule_active(uint64_t const molecule_bit_field,
                       int const mol_id)
{
    uint64_t const one = 1;
    return molecule_bit_field & (one << (mol_id-1));
}
