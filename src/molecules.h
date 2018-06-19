#ifndef MOLECULES_H_
#define MOLECULES_H_

/*Note: the MoleculeNumber enum is structured so that the molecule ids used
  in this model are equal to the hitran ids - 1.  If this changes, then the
  body of the HITRAN_id_to_model_id function must also change.*/

/*The constants defined in this enum are used as indices into an array, so
  make sure they start at zero and only increment by one.*/
typedef enum MoleculeNumber
{
    H2O = 0,
    CO2,
    O3,
    N2O,
    CO,
    CH4,
    O2,
    NUM_MOL
} MoleculeNumber_t;

int get_mol_name(int const id,
                 char *name,
                 int const name_len);

int HITRAN_id_to_model_id(int hitran_id,
                          int *model_id);

#endif
