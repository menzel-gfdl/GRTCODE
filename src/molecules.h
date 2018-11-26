#ifndef MOLECULES_H_
#define MOLECULES_H_

#include <stdint.h>
#include "parse_HITRAN_file.h"


#define MOL_NAME_LEN 8


/*Note: these numbers are assigned by the HITRAN database.*/
typedef enum HitranMoleculeId
{
    H2O = 1,
    CO2 = 2,
    O3 = 3,
    N2O = 4,
    CO = 5,
    CH4 = 6,
    O2 = 7,
    NO = 8,
    SO2 = 9,
    NO2 = 10,
    NH3 = 11,
    HNO3 = 12,
    OH = 13,
    HF = 14,
    HCl = 15,
    HBr = 16,
    HI = 17,
    ClO = 18,
    OCS = 19,
    H2CO = 20,
    HOCl = 21,
    N2 = 22,
    HCN = 23,
    CH3Cl = 24,
    H2O2 = 25,
    C2H2 = 26,
    C2H6 = 27,
    PH3 = 28,
    COF2 = 29,
    SF6 = 30,
    H2S = 31,
    HCOOH = 32,
    HO2 = 33,
    O = 34,
    ClONO2 = 35,
    NOp = 36,
    HOBr = 37,
    C2H4 = 38,
    CH3OH = 39,
    CH3Br = 40,
    CH3CN = 41,
    CF4 = 42,
    C4H2 = 43,
    HC3N = 44,
    H2 = 45,
    CS = 46,
    SO3 = 47,
    C2N2 = 48,
    COCl2 = 49,
    SO = 50,
    C3H4 = 51,
    CH3 = 52,
    CS2 = 53,
    NUM_MOLS = 53
} HitranMoleculeId_t;


/** @brief Molecule type.*/
typedef struct Molecule
{
    char name[MOL_NAME_LEN]; /**< Name.*/
    int id; /**< HITRAN id.*/
    fp_t mass; /**< Mass [g].*/
    int num_isotopologues; /**< Number of isotopologues.*/
    LineParams_t line_params; /**< Line parameters.*/
} Molecule_t;


/** @brief Initialize a molecule.  Read in its line parameters from the
           input HITRAN database file.
    @return SUCCESS or an error code.*/
int molecule(Molecule_t * const mol, /**< Molecule object.*/
             int const id, /**< HITRAN molecule id.*/
             char const * const hitran_path, /**< Path to HITRAN database
                                                  file.*/
             double const min_line_center, /**< Lower bound [1/cm] for
                                                spectral line centers.*/
             double const max_line_center /**< Upper bound [1/cm] for
                                               spectral line centers.*/
            );


/** @brief Free memory stored by the molecule object.
    @return SUCCESS or an error code.*/
int free_molecule(Molecule_t * const mol /**< Molecule.*/
                 );


int molecule_hash(int const mol_id,
                  int * const hash);


int activate_molecule(uint64_t * const molecule_bit_field,
                      int const mol_id);


int is_molecule_active(uint64_t const molecule_bit_field,
                       int const mol_id);


#endif
