#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "molecules.h"


static int const MIN_NAME_LEN = 8;


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
    switch(id)
    {
        case H2O:
            snprintf(name,name_len,"H2O");
            break;
        case CO2:
            snprintf(name,name_len,"CO2");
            break;
        case O3:
            snprintf(name,name_len,"O3");
            break;
        case N2O:
            snprintf(name,name_len,"N2O");
            break;
        case CO:
            snprintf(name,name_len,"CO");
            break;
        case CH4:
            snprintf(name,name_len,"CH4");
            break;
        case O2:
            snprintf(name,name_len,"O2");
            break;
        case NO:
            snprintf(name,name_len,"NO");
            break;
        case SO2:
            snprintf(name,name_len,"SO2");
            break;
        case NO2:
            snprintf(name,name_len,"NO2");
            break;
        case NH3:
            snprintf(name,name_len,"NH3");
            break;
        case HNO3:
            snprintf(name,name_len,"HNO3");
            break;
        case OH:
            snprintf(name,name_len,"OH");
            break;
        case HF:
            snprintf(name,name_len,"HF");
            break;
        case HCl:
            snprintf(name,name_len,"HCl");
            break;
        case HBr:
            snprintf(name,name_len,"HBr");
            break;
        case HI:
            snprintf(name,name_len,"HI");
            break;
        case ClO:
            snprintf(name,name_len,"ClO");
            break;
        case OCS:
            snprintf(name,name_len,"OCS");
            break;
        case H2CO:
            snprintf(name,name_len,"H2CO");
            break;
        case HOCl:
            snprintf(name,name_len,"HOCl");
            break;
        case N2:
            snprintf(name,name_len,"N2");
            break;
        case HCN:
            snprintf(name,name_len,"HCN");
            break;
        case CH3Cl:
            snprintf(name,name_len,"CH3Cl");
            break;
        case H2O2:
            snprintf(name,name_len,"H2O2");
            break;
        case C2H2:
            snprintf(name,name_len,"C2H2");
            break;
        case C2H6:
            snprintf(name,name_len,"C2H6");
            break;
        case PH3:
            snprintf(name,name_len,"PH3");
            break;
        case COF2:
            snprintf(name,name_len,"COF2");
            break;
        case SF6:
            snprintf(name,name_len,"SF6");
            break;
        case H2S:
            snprintf(name,name_len,"H2S");
            break;
        case HCOOH:
            snprintf(name,name_len,"HCOOH");
            break;
        case HO2:
            snprintf(name,name_len,"HO2");
            break;
        case O:
            snprintf(name,name_len,"O");
            break;
        case ClONO2:
            snprintf(name,name_len,"ClONO2");
            break;
        case NOp:
            snprintf(name,name_len,"NO+");
            break;
        case HOBr:
            snprintf(name,name_len,"HOBr");
            break;
        case C2H4:
            snprintf(name,name_len,"C2H4");
            break;
        case CH3OH:
            snprintf(name,name_len,"CH3OH");
            break;
        case CH3Br:
            snprintf(name,name_len,"CH3Br");
            break;
        case CH3CN:
            snprintf(name,name_len,"CH3CN");
            break;
        case CF4:
            snprintf(name,name_len,"CF4");
            break;
        case C4H2:
            snprintf(name,name_len,"C4H2");
            break;
        case HC3N:
            snprintf(name,name_len,"HC3N");
            break;
        case H2:
            snprintf(name,name_len,"H2");
            break;
        case CS:
            snprintf(name,name_len,"CS");
            break;
        case SO3:
            snprintf(name,name_len,"SO3");
            break;
        case C2N2:
            snprintf(name,name_len,"C2N2");
            break;
        case COCl2:
            snprintf(name,name_len,"COCl2");
            break;
        case SO:
            snprintf(name,name_len,"SO");
            break;
        case C3H4:
            snprintf(name,name_len,"C3H4");
            break;
        case CH3:
            snprintf(name,name_len,"CH3");
            break;
        case CS2:
            snprintf(name,name_len,"CS2");
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
    if (mol_id < H2O || mol_id > NUM_MOLS)
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
