#ifndef COLLISION_INDUCED_ABSORPTION_H_
#define COLLISION_INDUCED_ABSORPTION_H_

#include <stdint.h>
#include "floating_point_type.h"


#define CIA_NAME_LEN 8

/*Number of possible CIA combinations (i.e. N2-N2, N2-O2, O2-O2).*/
#define MAX_NUM_CIAS 3


/*Collision-induced absorption identifiers.*/
typedef enum CiaId
{
    CIA_N2 = 0,
    CIA_O2,
    NUM_CIAS,
} CiaId_t;


/** @brief Collision-induced absorption object.*/
typedef struct CollisionInducedAbsorption
{
    int id[2]; /**< Id of molecules.*/
    char *name[2]; /**< Molecule names.*/
    char name_buf[2*CIA_NAME_LEN]; /**< Buffer for molecule names.*/
    fp_t *cross_section; /**< Cross-section [cm^4].*/
    uint64_t num_wpoints; /**< Spectral grid size.*/
    int gpu_id; /**< Device id.*/
} CollisionInducedAbsorption_t;


/** @brief Read in the collision-induced absorption cross sections.
    @return RS_SUCCESS or an error code.*/
int get_collision_induced_cross_sections(CollisionInducedAbsorption_t * const cia, /**< Collision-induced absorption object.*/
                                         int const id[2], /**< Ids of species.*/
                                         char const * const filepath, /**< Cross section csv file.*/
                                         uint64_t const num_wpoints, /**< Spectral grid size.*/
                                         double const w0, /**< Spectral grid lower bound [1/cm].*/
                                         double const res, /**< Spectral resolution [1/cm].*/
                                         int const gpu_id /**< Device id.*/
                                        );


/** @brief Free memory reserved for collision-induced absorption cross sections.
    @return RS_SUCCESS or an error code.*/
int free_collision_induced_cross_sections(CollisionInducedAbsorption_t * const cia /**< Collision-induced absorption object.*/
                                         );


#endif
