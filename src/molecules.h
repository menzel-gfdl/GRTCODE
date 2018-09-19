#ifndef MOLECULES_H_
#define MOLECULES_H_


/*Note: these numbers are assigned by the HITRAN database.*/
typedef enum HitranMoleculeId
{
    Hitran_H2O = 1,
    Hitran_CO2 = 2,
    Hitran_O3 = 3,
    Hitran_N2O = 4,
    Hitran_CO = 5,
    Hitran_CH4 = 6,
    Hitran_O2 = 7
} HitranMoleculeId_t;


/*Note: The must be powers of 2 starting at 1, since they are used to turn
  on unique bits in a bit field.*/
typedef enum MoleculeNumber
{
    H2O = 1,
    CO2 = 2,
    O3 = 4,
    N2O = 8,
    CO = 16,
    CH4 = 32,
    O2 = 64,
    NUM_MOLS = 7
} MoleculeNumber_t;


int get_mol_name(int const id,
                 char * const name,
                 int const name_len);


int HITRAN_id_to_model_id(int const hitran_id,
                          int * const model_id);


int molecule_hash(int const mol_id,
                  int * const hash);


int is_molecule_active(int const molecule_bit_field,
                       int const mol_id);


#endif
