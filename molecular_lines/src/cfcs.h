/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef CFCS_H_
#define CFCS_H_

#include <stdint.h>
#include "floating_point_type.h"


#define CFC_NAME_LEN 16


/*CFC identifiers.*/
typedef enum CfcId
{
    CFC11 = 0,
    CFC12,
    CFC113,
    CFC114,
    CFC115,
    HCFC22,
    HCFC141b,
    HCFC142b,
    HFC23,
    HFC125,
    HFC134a,
    HFC143a,
    HFC152a,
    HFC227ea,
    HFC245fa,
    CCl4,
    C2F6,
    CF4,
    CH2Cl2,
    NF3,
    SF6,
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
