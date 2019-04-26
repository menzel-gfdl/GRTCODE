#ifndef CFCS_H_
#define CFCS_H_

#include <stdint.h>
#include "floating_point_type.h"


#define CFC_NAME_LEN 8


/*CFC identifiers.*/
typedef enum CfcId
{
    F11 = 0,
    F12,
    NUM_CFCS
} CfcId_t;


/*Container for CFC cross section values.*/
typedef struct CfcCrossSection
{
    fp_t *cross_section; /*CFC cross-section [cm^2].*/
    int id; /*CFC identifier.*/
    char name[CFC_NAME_LEN]; /*CFC name.*/
    uint64_t num_wpoints; /*Spectral grid size.*/
    int gpu_id; /*Device id.*/
} CfcCrossSection_t;


/*Read in the cfc cross sections.*/
int get_cfc_cross_sections(CfcCrossSection_t *xsc, /*CFC cross section object.*/
                           int const id, /*CFC identifier.*/
                           char const * const filepath, /*CSV file with cross section values.*/
                           uint64_t const num_wpoints, /*Spectral grid size.*/
                           double const w0, /*Spectral grid lower bound [1/cm].*/
                           double const res, /*Spectral grid resolution [1/cm].*/
                           int const gpu_id /*Device id.*/
                          );


/*Free memory.*/
int free_cfc_cross_sections(CfcCrossSection_t *xsc /*CFC cross section object.*/
                           );


#endif
