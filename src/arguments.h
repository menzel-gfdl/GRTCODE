#ifndef ARGUMENTS_H_
#define ARGUMENTS_H_

#include "molecules.h"

#define minNhitfiles 1
#define maxNhitfiles NUM_MOL


struct arguments
{
    char *atmosInputFile; /*Input atmosphere netCDF file.*/
    char *hitFiles[maxNhitfiles]; /*Input HITRAN .par files.*/
    int nHitFiles; /*Total number of inputted HITRAN files.*/
    int host; /*Flag for host-only execution.*/
    int device; /*Specific device id to run on.*/
    int t; /*Starting time dimension index, inclusive.*/
    int T; /*Ending time dimension index, inclusive.*/
    int x; /*Starting longitude index, inclusive.*/
    int X; /*Ending longitude index, inclusive.*/
    int y; /*Starting latitude index, inclusive.*/
    int Y; /*Ending longitude index, inclusive.*/
    int w; /*Wavenumber lower bound (1/cm), inclusive.*/
    int W; /*Wavenumber upper bound (1/cm), inclusive.*/
    int h2o_ctm; /*Flag for including the water vapor continuum.*/
    int o3_ctm; /*Flag for including the ozone continuum.*/
    double res; /*Wavenumber resolution (1/cm).*/
    int wingBreadth; /*Wings cutoff (1/cm).*/
    double molConc[NUM_MOL]; /*Molecular concentrations (ppmv).  Each spot
                               in the array corresponds to a molecule listed
                               in the MoleculeNumber_t enum defined in
                               molecules.h.*/
    char *outputFile; /*Output netCDF file.*/
    int workers; /*Number of workers (threads if host only, or devices if
                   run on GPUs.*/
    int write_spectra; /*Flag for writing out optical depth values.*/
};


void parse_options(int argc,
                   char **argv,
                   struct arguments *arguments);


#endif
