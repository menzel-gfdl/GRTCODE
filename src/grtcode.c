/*
    GRTCODE is a GPU-able Radiative Transfer Code
    Copyright (C) 2016  Garrett Wright

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; version 2.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
    USA.
*/

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include <argp.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "continuum.h"
#include "eval_gamma.h"
#include "eval_profile.h"
#include "eval_pShift.h"
#include "eval_Snn_correction.h"
#include "GasProps.h"
#include "GaussianFuncs.h"
#include "grtcode.h"
#include "IdaVoigtFuncs.h"
#include "LineShapeUtils.h"
#include "LorentzFuncs.h"
#include "outputNetcdfSpec.h"
#include "parseHITRANfile.h"
#include "parseNetcdfRadiation.h"
#include "pre_eval_Snn.h"
#include "RfmVoigtFuncs.h"
#include "TIPS_2011.h"

/*---------------------------------------------------------------------------*/
/*Helper data structures.*/

/*These are equal to HITRAN_MOLID-1.*/
typedef enum MoleculeNumber_t
{
    H2O = 0,
    CO2 = 1,
    O3  = 2,
    N2O = 3,
    CO  = 4,
    CH4 = 5,
    O2  = 6,
    NUM_MOL = 7
} MoleculeNumber_t;

/*---------------------------------------------------------------------------*/
/*Set argp variables.*/

const char *argp_program_version = "lbl-dev 0.1";
const char *argp_program_bug_address = "<raymond.menzel@noaa.gov>";
static char doc[] = "GFDL style documentation goes >/\n\n\\"
                        "<^here.\n\v"
                        "Other Documentation goes here.";
#define minNhitfiles 1
#define maxNhitfiles NUM_MOL
static const unsigned int minNargs=minNhitfiles;
static const unsigned int maxNargs=maxNhitfiles;
static char args_doc[] = "-aINPUT.nc -fFORMAT -oOUT.nc"
                             " [molecule concentration specifications]"
                             " HITFILES";

/*Command line options.*/
static struct argp_option options[] =
{
    {"verbose",
     'v',
     0,
     0,
     "Opens the elevator door"},

    {"quiet",
     'q',
     0,
     0,
     "Closes the elevator door"},

    {"silent",
     's',
     0,
     OPTION_ALIAS},

    {"device",
     'd',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Use gpu implementation on specifed DEVICE."
         "\n\tDefault DEVICE is simply GPU0"
         "\n\tIncompatible with --host."
         "\n\t --mpi modifies this flag to prescribe numDevices per node."},

    {"host",
     'h',
     0,
     0,
     "Use HOST cpu implementation. \n\t(incompatible with --device)"},

    {"mpi",
     'M',
     0,
     0,
     "Use MPI: Ranks taken from MPI_Comm_World. (Must be compiled for MPI!)"},

    {"output",
     'o',
     "FILE",
     0,
     "Output to FILE"},

    {"atmos",
     'a',
     "INPUT.NC",
     0,
     "NC file containing model atmosphere."},

    {"atmos_format",
     'f',
     "FORMAT",
     0,
     "Format of the input atmosphere file.  Allowed values are 'rfmip' and"
         " 'gfdl'."},

    {"minw",
     'w',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "minimum Wavenumber (lower bound, inclusive), defaults 1",
     -3},

    {"maxw",
     'W',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "maximum Wavenumber (upper bound, inclusive), defaults 50000",
     -3},

    {"mint",
     't',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "minimum time (lower bound, inclusive), defaults 0",
     -3},

    {"maxt",
     'T',
     "VAL",
     OPTION_ARG_OPTIONAL, 
     "maximum time (upper bound, inclusive), defaults 0",
     -3},

    {"res",
     'r',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Resolution (wavenumber), defaults 1.0",
     -3},

    {"wings",
     'c',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Wings cutoff (+/- integer wavenumber), defaults 25",
     -3},

    {"h2o",
     '1',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Water Concentration.  Default reads from INPUT.NC, else supply global"
         " value (ppmv).  Layer partial pressure = (layer pressure)*"
         "(water concentration/10^6).",
     -2},

    {"co2",
     '2',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Carbon Dioxide Concentration.  Default reads from INPUT.NC, else"
         " supply global value (ppmv).  Layer partial pressure ="
         " (layer pressure)*(carbon dioxide concentration/10^6).",
     -2},

    {"o3",
     '3',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Ozone Concentration.  Default reads from INPUT.NC, else supply global"
         " value (ppmv).  Layer partial pressure = (layer pressure)*"
         "(ozone concentration/10^6).",
     -2},

    {"n2o",
     '4',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Nitrous Oxide.  Default reads from INPUT.NC, else supply global value"
         " (ppmv).  Layer partial pressure = (layer pressure)*"
         "(nitrous oxide concentration/10^6).",
     -2},

    {"co",
     '5',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Carbon Monoxide Concentration.  Default read from INPUT.NC, else"
         " supply global value (ppmv).  Layer partial pressure ="
         " (layer pressure)*(carbon monoxide concentration/10^6).",
     -2},

    {"ch4",
     '6',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Methane Conentration.  Default read from INPUT.NC, else supply global"
         " value (ppmv).  Layer partial pressure = (layer pressure)*"
         "(methane concentration/10^6).",
     -2},

    {"o2",
     '7',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Oxygen Concentration.  Default read from INPUT.nc, else supply global"
         " value (ppmv).  Layer partial pressure = (layer_pressure)*"
         "(oxygen concentration/10^6).",
     -2},

    {"ctm",
     'C',
     0,
     0,
     "Enables the continuum codes for testing",
     -1},

    {0}
};

/*---------------------------------------------------------------------------*/
/*Command line arguments structure.*/
struct arguments
{
    char *atmos;                  /*Input atmosphere netCDF file.*/
    char *atmos_format;           /*Format of the inputted netCDF file.*/
    char *hitfiles[maxNhitfiles]; /*Input HITRAN .par files.*/
    int nhitfiles;                /*Total number of inputted HITRAN files.*/
    int nmolConc;                 /*Number of inputted molecular concentrations.*/
    int nmolConcOver;             /*Number of molecular concentrations that will be taken from the netCDF file.*/
    int silent;                   /*Use silent mode.*/
    int verbose;                  /*Use verbose mode.*/
    int host;                     /*Flag for host-only execution.*/
    int device;                   /*Specific device id to run on.*/
    int mpi;                      /*Flag for using mpi.*/
    int t;                        /*Starting time (seconds?), inclusive.*/
    int T;                        /*Ending time (seconds?), inclusive.*/
    int w;                        /*Wavenumber lower bound (cm?), inclusive.*/
    int W;                        /*Wavenumber upper bound (cm?), inclusive.*/
    int ctm;                      /*Flag for including continuum.*/
    double res;                   /*Wavenumber resolution (cm?).*/
    int wingBreadth;              /*Wings cutoff (wavenumber).  Must be an integer.*/
    double h2o;                   /*Water concentration (atm).*/
    double co2;                   /*Carbon dioxide concentration (atm).*/
    double o3;                    /*Ozone concentration (atm).*/
    double n2o;                   /*Nitrous oxide concentration (atm).*/
    double co;                    /*Carbon monoxide concentration (atm).*/
    double ch4;                   /*Methane concentration (atm).*/
    double o2;                    /*Oxygen concentration (atm).*/
    char *output_file;            /*Output netCDF file.*/
};

/*---------------------------------------------------------------------------*/
/*Helper parsing function for molecular concentrations.*/
static double parse_MolecConc(char *arg)
{
    /*Local variables*/
    double res = 0; /*molecular concentration.*/

    /*Make sure that the inputted pointer is not null.*/
    if (arg == NULL)
    {
        fprintf(stderr,
                "Error(parse_MolecConc): the inputted pointer is null.\n");
        exit(EXIT_FAILURE);
    }

    /*Get the inputted molecular concentration.*/
    if (arg[0] == 'a')
    {
        /*A leading 'a' character specifies that the concentration should be
         taken from the "nc" file.*/
        res = -1;
    }
    else if (isalpha(arg[0]))
    {
        fprintf(stderr,
                "Error(parse_MolecConc): the supplied character (%c) for"
                    " overriding a molecule concentration is not understood."
                    "  Review args.\n",
                arg[0]);
        exit(EXIT_FAILURE);
    }
    else
    {
        /*This should probably be changed to a strtod call with err checks
          later.*/
        res = atof(arg);
    }

    return res;
}

/*---------------------------------------------------------------------------*/
/*Argp options parser function.*/
static error_t parse_opt(int key,
                         char *arg,
                         struct argp_state *state)
{
    /*Local variables.*/
    struct arguments *arguments = NULL; /*Arguments pointer.*/

    /*Point to the inputted argument from argp_parse.*/
    arguments = (struct arguments*)(state->input);

    /*Store the inputted arguments.*/
    switch(key)
    {
        case 'q':case 's':
            arguments->silent = 1;
            break;
        case 'v':
            arguments->verbose = 1;
            break;
        case 'o':
            arguments->output_file = arg;
            break;
        case 'h':
            arguments->host = 1;
            break;
        case 'd':
            arguments->device = atoi(arg);
            break;
        case 'M':
            arguments->mpi = 1;
            break;
        case 'a':
            arguments->atmos = arg;
            break;
        case 'f':
            arguments->atmos_format = arg;
            break;
        case 'w':
            arguments->w = atoi(arg);
            break;
        case 'W':
            arguments->W = atoi(arg);
            break;
        case 't':
            arguments->t = atoi(arg);
            break;
        case 'T':
            arguments->T = atoi(arg);
            break;
        case 'r':
            arguments->res = atof(arg);
            break;      
        case 'c':
            arguments->wingBreadth = atoi(arg);
            break;
        case '1':      
            arguments->h2o = parse_MolecConc(arg);
            arguments->nmolConc++;
            break;
        case '2':
            arguments->co2 = parse_MolecConc(arg);
            arguments->nmolConc++;
            break;
        case '3':
            arguments->o3 = parse_MolecConc(arg);
            arguments->nmolConc++;
            break;
        case '4':
            arguments->n2o = parse_MolecConc(arg);
            arguments->nmolConc++;
            break;
        case '5':
            arguments->co = parse_MolecConc(arg);
            arguments->nmolConc++;
            break;
        case '6':
            arguments->ch4 = parse_MolecConc(arg);
            arguments->nmolConc++;
            break;
        case '7':
            arguments->o2 = parse_MolecConc(arg);
            arguments->nmolConc++;
            break;
        case 'C':
            arguments->ctm = 1;
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num >= maxNargs)
            {
                fprintf(stderr,
                        "Error(parse_opt): there are too many command line"
                            " arguments.\n");
                argp_usage(state);
            }
            arguments->hitfiles[state->arg_num] = arg;
            arguments->nhitfiles++;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < minNargs )
            {
                fprintf(stderr,
                        "Error(parse_opt): there are too few command line"
                            " arguments.\n");
                argp_usage( state );
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

/*---------------------------------------------------------------------------*/
/*Necessary argp struct.*/
static struct argp argp = {options,parse_opt,args_doc,doc};

/*---------------------------------------------------------------------------*/
/*Calculate parital pressures and number densities for the molecule designated
  by the inputted molid.

  Arguments:
      args      [in]      Pointer to the command line argument structure.
      molid     [in]      Id of the molecule whose partial pressures and number
                              densities will be calculated.
      atmosData [in,out]  Pointer to a structured containing molecule
                              specific atmospheric data.
      time      [in]      Size of the time dimension for the atmospheric
                              data arrays.
*/
static void checkMolConfig(struct arguments *args,
                           unsigned int const molid,
                           radiationOutputFields_t *atmosData,
                           int const time)
{
    /*Local variables*/
    REAL_t *PS = atmosData->PS;
    const size_t nlat = atmosData->nlat;
    const size_t nlon = atmosData->nlon;
    const size_t nlvl = atmosData->npfull;
    int in_input_file = 0;
    int is_rfmip = 0;
    double concentration;

    /*Determine if the input file was a rfmip formatted netcdf file.*/
    if (strcmp(args->atmos_format,"rfmip") == 0)
    {
        is_rfmip = 1;
    }

    if (molid == 1 || molid == 3)
    {
        /*Currently water vapor and ozone concentrations are contained
          in both kinds (rfmip and gfdl) input netcdf files.*/
        in_input_file = 1;
    }
    else
    {
        if (is_rfmip)
        {
            /*Currently co2, n2o, co, ch4, and o2 are only contained in
              gfdl formatted input netcdf files.*/
            in_input_file = 1;
        }
    }

    /*Store the inputted molecular concentration of the inputted molecule.*/
    switch(molid)
    {
        case 1:
            concentration = args->h2o;
            break;
        case 2:
            concentration = args->co2;
            break;
        case 3:
            concentration = args->o3;
            break;
        case 4:
            concentration = args->n2o;
            break;
        case 5:
            concentration = args->co;
            break;
        case 6:
            concentration = args->ch4;
            break;
        case 7:
            concentration = args->o2;
            break;
        default:
            fprintf(stderr,
                    "Error(checkMolConfig): this Hitfiles MolId (%d) does not"
                        " appear to be supported yet.\n",
                    molid);
            exit(EXIT_FAILURE);
    }

    if (concentration == 0)
    {
        fprintf(stderr,
                "Error(checkMolConfig): an inputted molecular concentration"
                    " of zero is not currently supported for this"
                    " Hitfiles MolId (%d).\n",
                molid);
        exit(EXIT_FAILURE);
    }
    else if (concentration < 0)
    {
        if (!in_input_file)
        {
            fprintf(stderr,
                    "Error(checkMolConfig): Hitfiles MolId's (%d)"
                        " concentration is not contained in gfdl formatted"
                        " input atmosphere files.  Please provide a value"
                        " on the command line.\n.",
                    molid);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        /*Calculate the molecule's partial pressure from the inputted
          concentration (ppmv) value.*/
        setGlobalPartialPres(concentration,
                             PS,
                             atmosData->P,
                             molid,
                             time,
                             nlat,
                             nlon,
                             (size_t)NUM_MOL,
                             nlvl);
    }

    /*Calculate the number densities for the molecule.*/
    setGlobalNumberDensity(atmosData->N,
                           PS,
                           atmosData->T,
                           molid,
                           time,
                           nlat,
                           nlon,
                           (size_t)NUM_MOL,
                           nlvl);

    return;
}

/*---------------------------------------------------------------------------*/
/*Set some constants.*/

/*Set the maximum number of spectra lines to a number divisible by 32 to
  keep global device memory accesses aligned.  Here we choose 2^19.*/
#define MAX_NUM_SPECTRAL_LINES 524288

/*Set the maximum number of CUDA streams.*/
#ifndef MAXNSTREAMS
#define MAXNSTREAMS 2
#endif

/*---------------------------------------------------------------------------*/
/*Include some GPU helper functions.*/

#ifdef __NVCC__

#undef FORCE_KERNEL_CHECK
/* #define FORCE_KERNEL_CHECK */

#undef EVENTS
/* #define EVENTS */

#include "cudaHelpers.cuh"

#endif

/*---------------------------------------------------------------------------*/
/*Calculate the dimensionless optical depth values at all appropriate heights
  and frequencies for the molecule associated with the inputted molId.

  Arguments:
      molId        [in]
      nL           [in]
      loWn         [in]
      nF           [in]
      resolution   [in]
      breadth      [in]
      numLayers    [in]
      out_h        [in,out]
      T_h          [in]
      P_h          [in]
      iso          [in]
      Vnn          [in]
      Snn_ref      [in]
      Yair         [in]
      Yself        [in]
      En           [in]
      n            [in]
      d            [in]
      TauU_h       [in]
      pathLength_h [in]
      PS_h         [in]
      Lines_h      [in,out]
      OptBuf_h     [in,out]

  Return:
      EXIT_SUCCESS if the function completes normally.
*/
int host_optics_perMol(uint8_t const molId,
                       unsigned int const nL,
                       REAL_t const loWn,
                       unsigned int const nF,
                       REAL_t const resolution,
                       int const breadth,
                       unsigned int const numLayers,
                       REAL_t *out_h,
                       REAL_t const * const T_h,
                       REAL_t const * const P_h,
                       uint8_t const * const iso,
                       REAL_t const * const Vnn,
                       REAL_t const * const Snn_ref,
                       float const * const Yair,
                       float const * const Yself,
                       float const * const En,
                       float const * const n,
                       float const * const d,
                       REAL_t const * const TauU_h,
                       REAL_t const * const pathLength_h,
                       REAL_t const * const PS_h,
                       RefLinePtrs_t* const Lines_h,
                       OpticsBufPtrs_t* const OptBuf_h)
{

    /*Initilize TIPS.*/
    initTIPS();

    /*Print out the current molecules HITRAN molecule id.*/
    printf("HITRAN molId:=%d\n",
           molId);

    /*Point at the arrays located in the inputted Lines_h structure.*/
    uint8_t *iso_h = Lines_h->iso;
    REAL_t *Vnn_h = Lines_h->Vnn;
    REAL_t *Snn_ref_h = Lines_h->Snn_ref;
    float *Yair_h = Lines_h->Yair;
    float *Yself_h = Lines_h->Yself;
    float *En_h = Lines_h->En;
    float *n_h = Lines_h->n;
    float *d_h = Lines_h->d;

    /*Copy the inputted arrays into the arrays located in the inputted
      Lines_h structure.*/
    memcpy(iso_h,
           iso,
           nL*sizeof(uint8_t));
    memcpy(Vnn_h,
           Vnn,
           nL*sizeof(REAL_t));
    memcpy(Snn_ref_h,
           Snn_ref,
           nL*sizeof(REAL_t));
    memcpy(Yair_h,
           Yair,
           nL*sizeof(float));
    memcpy(Yself_h,
           Yself,
           nL*sizeof(float));
    memcpy(En_h,
           En,
           nL*sizeof(float));
    memcpy(n_h,
           n,
           nL*sizeof(float));
    memcpy(d_h,
           d,
           nL*sizeof(float));

    /*Point at arrays located in the inputted OptBuf_h structure.*/
    REAL_t *Gam_h = OptBuf_h->Gam;
    REAL_t *PShift_h = OptBuf_h->PShift;
    REAL_t *S_h = OptBuf_h->S;

    /*Execute the pre_eval_Snn kernel.*/
    printf("%s Kernel Launch.. %d lines ... ",
           "pre_eval_Snn_h",
           nL);
    pre_eval_Snn_h(nL,
                   molId,
                   iso_h,
                   Vnn_h,
                   En_h,
                   Snn_ref_h);
    printf("..done!\n");

    /*Execute the gamma kernel.*/
    printf("%s Kernel Launch.. %d lines, %d spatial points ... ",
           "eval_gamma_h",
           nL,
           numLayers);
    eval_gamma_h(numLayers,
                 nL,
                 P_h,
                 T_h,
                 PS_h,
                 Yself_h,
                 Yair_h,
                 n_h,
                 Gam_h);
    printf("..done!\n");

    /*Execute the pshift kernel.*/
    printf("%s Kernel Launch.. %d lines, %d spatial points ... ",
           "eval_pShift_h",
           nL,
           numLayers);
    eval_pShift_h(numLayers,
                  nL,
                  P_h,
                  Vnn_h,
                  d_h,
                  PShift_h);
    printf("..done!\n");

    /*Execute the Snn kernel.*/
    printf("%s Kernel Launch.. %d lines, %d spatial points ... ",
           "eval_Snn_correction_h",
           nL,
           numLayers);
    eval_Snn_correction_h(numLayers,
                          nL,
                          molId,
                          T_h,
                          iso_h,
                          Vnn_h,
                          En_h,
                          Snn_ref_h,
                          S_h);
    printf("..done!\n");

    /*Execute the profile kernel.*/
    printf("%s Kernel Launch.. %d lines, %d fsamples ... ",
           "eval_profile_h",
           nL,
           nF);
    eval_profile_h(molId,
                   nL,
                   nF,
                   loWn,
                   resolution,
                   numLayers,
                   breadth,
                   T_h,
                   Gam_h,
                   PShift_h,
                   S_h,
                   pathLength_h,
                   TauU_h,
                   out_h);
    printf("..done!\n");

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Allocate space for the arrays contained in the RefLinePtrs_t and
  OpticsBufPtrs_t structures.

  Arguments:
      numLayers   [in]
      numBufLines [in]
      L_h         [in,out]
      OpticBuf_h  [in,out]
      flags       [in]

  Return:
      EXIT_SUCCESS if the function completes normally.
*/
int host_optics_init(unsigned int const numLayers,
                     unsigned int const numBufLines,
                     RefLinePtrs_t * L_h,
                     OpticsBufPtrs_t * OpticBuf_h,
                     RefLine_flags_t const flags)
{
    /*Allocate arrays in the RefLinePtrs_t structure.*/
    *L_h = allocHost(numBufLines,
                     flags);

    /*Allocate arrays in the OpticsBufPtrs_t structure.*/
    OpticBuf_h->Gam = (REAL_t *)malloc(numLayers*numBufLines*sizeof(REAL_t));
    OpticBuf_h->PShift = (REAL_t *)malloc(numLayers*numBufLines*sizeof(REAL_t));
    OpticBuf_h->S = (REAL_t *)malloc(numLayers*numBufLines*sizeof(REAL_t));

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Free arrays contained in the RefLinePtrs_t and OpticsBufPtrs_t structures.

  Arguments:
      L_h        [in,out]
      OpticBuf_h [in,out]
      flags      [in]

  Return:
      EXIT_SUCCESS if the function completes normally.
*/
int host_optics_free(RefLinePtrs_t *L_h,
                     OpticsBufPtrs_t *OpticBuf_h,
                     RefLine_flags_t const flags)
{
    /*Free arrays in the RefLinePtrs_t structure.*/
    freeHost(*L_h,
             flags);

    /*Free arrays in the OpticsBufPtrs_t structure*/
    free(OpticBuf_h->Gam);
    free(OpticBuf_h->PShift);
    free(OpticBuf_h->S);

    return EXIT_SUCCESS;
}
/*---------------------------------------------------------------------------*/
/*Loop through each molecule, and calculate the optical depth values at
  each height and frequency.

  Arguments:
      numMols     [in]
      L           [in]
      loWn        [in]
      nF          [in]
      resolution  [in]
      wingBreadth [in]
      atmosData   [in]
      out         [in,out]
      time        [in]
      lat         [in]
      lon         [in]

  Return:
      EXIT_SUCCESS if the function completes normally.
*/
int host_launch(unsigned int const numMols,
                RefLinePtrs_t const * const L,
                REAL_t const loWn,
                unsigned int const nF,
                REAL_t const resolution,
                unsigned int const wingBreadth,
                radiationOutputFields_t *atmosData,
                REAL_t * const out,
                int const time,
                int const lat,
                int const lon)
{
    /*Local variables*/
    unsigned int mol;
    unsigned int numLayers=atmosData->npfull;
    OpticsBufPtrs_t OpticsBuf_h;
    RefLinePtrs_t LinesBuf_h;

    /*Sanity check.*/
    assert(numMols<=NUM_MOL);

    /*Point the beginning of the appropriate "column" of data.*/
    const size_t idx = (time*atmosData->nlon*atmosData->nlat +
                        lat*atmosData->nlon + lon)*atmosData->npfull;
    REAL_t *T = &(atmosData->T[idx]);
    REAL_t *P = &(atmosData->P[idx]);
    REAL_t *DELTAZ = &(atmosData->DELTAZ[idx]);
    REAL_t *N = &(atmosData->N[idx*NUM_MOL]);
    REAL_t *PS = &(atmosData->PS[idx*NUM_MOL]);

    /*Set the flags used to malloc/free arrays in a RefLinePtr_t structure.
      {(unsigned int)-1,1,0} = host cuda malloc default, host=True, device=false.*/
    RefLine_flags_t flags= {((unsigned int) -1),1,0};

    /*Allocate arrays.*/
    host_optics_init(numLayers,
                     MAX_NUM_SPECTRAL_LINES,
                     &LinesBuf_h,
                     &OpticsBuf_h,
                     flags);

    /*Loop over the molecules.*/
    for (mol=0;mol<numMols;++mol)
    {
        /*Calcuate the optical depth values.*/
        host_optics_perMol(L[mol].mol,
                           L[mol].nLines,
                           loWn,
                           nF,
                           resolution,
                           wingBreadth,
                           numLayers,
                           out,
                           T,
                           P,
                           L[mol].iso,
                           L[mol].Vnn,
                           L[mol].Snn_ref,
                           L[mol].Yair,
                           L[mol].Yself,
                           L[mol].En,
                           L[mol].n,
                           L[mol].d,
                           &(N[(L[mol].mol-1)*atmosData->npfull]),
                           DELTAZ,
                           &(PS[(L[mol].mol-1)*atmosData->npfull]),
                           &(LinesBuf_h),
                           &(OpticsBuf_h)); 
    }

    /*Free arrays.*/
    host_optics_free(&LinesBuf_h,
                     &OpticsBuf_h,
                     flags);

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/

#ifdef __NVCC__
/*---------------------------------------------------------------------------*/
/*Malloc device arrays and copy in data from the host.

  Arguments:
      numLayers [in]      Size of the height dimension of the arrays.
      numMols   [in]      Size of the molecule dimension of the arrays.
      T_h       [in]      Array of temperatures (K).  This array is stored as
                              [height].
      P_h       [in]      Array of pressures (atm).  This array is stored as
                              [height].
      N_h       [in]      Array of number densities (cm^-3).  This array is
                              stored as [molecule][height].
      Z_h       [in]      Array of layer thicknesses (cm).  This array is
                              stored as [height].
      Ps_h      [in]      Array of partial pressures (atm).  This array is
                              stored as [molecule][height].
      T_d       [in,out]  Array of temperatures (K).  This array is stored as
                              [height].
      P_d       [in,out]  Array of pressures (atm).  This array is stored as
                              [height].
      N_d       [in,out]  Array of number densities (cm^-3).  This array is
                              stored as [molecule][height].
      Z_d       [in,out]  Array of layer thicknesses (cm).  This array is
                              stored as [height].
      Ps_d      [in,out]  Array of partial pressures (atm).  This array is
                              stored as [molecule][height].

  Return:
      EXIT_SUCCESS if the function completes normally.
*/
__host__
int device_atmos_init(unsigned int const numLayers,
                      unsigned int const numMols,
                      REAL_t const * const T_h,
                      REAL_t const * const P_h,
                      REAL_t const * const N_h,
                      REAL_t const * const Z_h,
                      REAL_t const * const Ps_h,
                      REAL_t **T_d,
                      REAL_t **P_d,
                      REAL_t **N_d,
                      REAL_t **Z_d,
                      REAL_t **Ps_d)
{
    printf("\nInitializing atmosphere on device:\n");

    /*Malloc device arrays.*/
    printf("\tMallocing atmospheric data on device..");
    HANDLE_ERROR(cudaMalloc(T_d,
                            (numLayers)*sizeof(REAL_t)));
    HANDLE_ERROR(cudaMalloc(P_d,
                            numLayers*sizeof(REAL_t)));
    HANDLE_ERROR(cudaMalloc(N_d,
                            numMols*numLayers*sizeof(REAL_t)));
    HANDLE_ERROR(cudaMalloc(Z_d,
                            numLayers*sizeof(REAL_t)));
/*
    HANDLE_ERROR(cudaMalloc(Ps_d,
                            (numLayers+1)*numMols*sizeof(REAL_t)));
*/
    HANDLE_ERROR(cudaMalloc(Ps_d,
                            (numLayers)*numMols*sizeof(REAL_t)));
    printf(".done!\n");

    /*Copy data from the host to the device.*/
    printf("\tCopying atmospheric data to device..\n");
    printf("\tMemcpy..T_h -> T_d");
    HANDLE_ERROR(cudaMemcpy(*T_d,
                            T_h,
                            numLayers*sizeof(REAL_t),
                            cudaMemcpyHostToDevice));
    printf(".done!\n");
    printf("\tMemcpy..P_h -> P_d");
    HANDLE_ERROR(cudaMemcpy(*P_d,
                            P_h,
                            numLayers*sizeof(REAL_t),
                            cudaMemcpyHostToDevice));
    printf(".done!\n");
    printf("\tMemcpy.. N_h -> N_d");
    HANDLE_ERROR(cudaMemcpy(*N_d,
                            N_h,
                            numMols*numLayers*sizeof(REAL_t),
                            cudaMemcpyHostToDevice));
    printf(".done!\n");
    printf("\tMemcpy.. Z_h -> Z_d");
    HANDLE_ERROR(cudaMemcpy(*Z_d,
                            Z_h,
                            numLayers*sizeof(REAL_t),
                            cudaMemcpyHostToDevice));
    printf(".done!\n");
    printf("\tMemcpy.. Ps_h -> Ps_d");
    HANDLE_ERROR(cudaMemcpy(*Ps_d,
                            Ps_h,
                            numLayers*NUM_MOL*sizeof(REAL_t),
                            cudaMemcpyHostToDevice));
    printf(".done!\n");

    printf("Initializing atmosphere on device successful.\n\n");

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Free device arrays.

  Arguments:
      T_d  [in,out]  Array of temperatures (K).  This array is stored as
                         [height].
      P_d  [in,out]  Array of pressures (atm).  This array is stored as
                         [height].
      N_d  [in,out]  Array of number densities (cm^-3).  This array is
                         stored as [molecule][height].
      Z_d  [in,out]  Array of layer thicknesses (cm).  This array is stored
                         as [height].
      Ps_d [in,out]  Array of partial pressures (atm).  This array is stored
                         as [molecule][height].

  Return:
      EXIT_SUCCESS if the function completes normally.
*/
__host__
int device_atmos_free(REAL_t *T_d,
                      REAL_t *P_d,
                      REAL_t *N_d,
                      REAL_t *Z_d,
                      REAL_t *Ps_d)
{
    /*Free arrays.*/
    HANDLE_ERROR(cudaFree(T_d));
    HANDLE_ERROR(cudaFree(P_d));
    HANDLE_ERROR(cudaFree(N_d));
    HANDLE_ERROR(cudaFree(Z_d));
    HANDLE_ERROR(cudaFree(Ps_d));

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Malloc space for the device arrays contained in the RefLinePtrs_t and
  OpticsBufPtrs_t structures.

  Arguments:
      numLayers   [in]      Size of the height dimension of the arrays.
      numBufLines [in]      Size of the line dimension of the arrays.
      L_d         [in,out]  Pointer to a structure containing device
                                arrays stored as [line].
      OpticBuf_d  [in,out]  Pointer to a structure containing device
                                arrays stored as [height][line].

  Return:
      EXIT_SUCCESS if the function completes normally.
*/
__host__
int device_optics_init(unsigned int const numLayers,
                       unsigned int const numBufLines,
                       RefLinePtrs_t *L_d,
                       OpticsBufPtrs_t *OpticBuf_d)
{
#ifdef EVENTS
    /*Record the start time for a cuda event.*/
    cudaEvent_t start;
    cudaEvent_t stop;
    float elapsed;
    HANDLE_ERROR(cudaEventCreate(&start));
    HANDLE_ERROR(cudaEventCreate(&stop));
    HANDLE_ERROR(cudaEventRecord(start)); 
#endif

    /*Malloc device arrays contained in the RefLinesPtrs_t structure.*/
    *L_d = allocDevice(numBufLines);

    /*Malloc device arrays contained in the OpticsBufPtrs_t structure.*/
    HANDLE_ERROR(cudaMalloc(&(OpticBuf_d->Gam),
                            numLayers*numBufLines*sizeof(REAL_t)));
    HANDLE_ERROR(cudaMalloc(&(OpticBuf_d->PShift),
                            numLayers*numBufLines*sizeof(REAL_t)));
    HANDLE_ERROR(cudaMalloc(&(OpticBuf_d->S),
                            numLayers*numBufLines*sizeof(REAL_t)));

#ifdef EVENTS
    /*Record the stop time for the cuda event and calculate and print the
      elapsed time.*/
    HANDLE_ERROR(cudaEventRecord(stop));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    HANDLE_ERROR(cudaEventElapsedTime(&elapsed,
                                      start,
                                      stop));
    printf("Time to device_optics_init: %3.1f ms\n",
           elapsed);
#endif

#ifdef EVENTS
    /*Clean up cuda event objects.*/
    HANDLE_ERROR(cudaEventDestroy(start));
    HANDLE_ERROR(cudaEventDestroy(stop));
#endif

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Free device arrays contained in the RefLinePtrs_t and
  OpticsBufPtrs_t structures.

  Arguments:
      L_d         [in,out]  Pointer to a structure containing device
                                arrays stored as [line].
      OpticBuf_d  [in,out]  Pointer to a structure containing device
                                arrays stored as [height][line].

  Return:
      EXIT_SUCCESS if the function completes normally.
*/
__host__
int device_optics_free(RefLinePtrs_t* L_d,
                       OpticsBufPtrs_t* OpticBuf_d)
{
#ifdef EVENTS
    /*Record the start time for a cuda event.*/
    cudaEvent_t start;
    cudaEvent_t stop;
    float elapsed;
    HANDLE_ERROR(cudaEventCreate(&start));
    HANDLE_ERROR(cudaEventCreate(&stop));
    HANDLE_ERROR(cudaEventRecord(start)); 
#endif

    /*Free device arrays contained in the RefLinesPtrs_t structure.*/
    HANDLE_ERROR(cudaFree(L_d->iso));
    HANDLE_ERROR(cudaFree(L_d->Vnn));
    HANDLE_ERROR(cudaFree(L_d->Snn_ref));
    HANDLE_ERROR(cudaFree(L_d->Yair));
    HANDLE_ERROR(cudaFree(L_d->Yself));
    HANDLE_ERROR(cudaFree(L_d->En));
    HANDLE_ERROR(cudaFree(L_d->n));
    HANDLE_ERROR(cudaFree(L_d->d));

    /*Free device arrays contained in the OpticsBufPtrs_t structure.*/
    HANDLE_ERROR(cudaFree(OpticBuf_d->Gam));
    HANDLE_ERROR(cudaFree(OpticBuf_d->PShift));
    HANDLE_ERROR(cudaFree(OpticBuf_d->S));

#ifdef EVENTS
    /*Record the stop time for the cuda event and calculate and print the
      elapsed time.*/
    HANDLE_ERROR(cudaEventRecord(stop));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    HANDLE_ERROR(cudaEventElapsedTime(&elapsed,
                                      start,
                                      stop));
    printf("Time to device_optics_free: %3.1f ms\n",
           elapsed);
#endif

#ifdef EVENTS
    /*Clean up cuda event objects.*/
    HANDLE_ERROR(cudaEventDestroy(start));
    HANDLE_ERROR(cudaEventDestroy(stop));
#endif

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
__host__
int device_optics_perMol(cudaStream_t stream,
                         uint8_t const molId,
                         unsigned int const nL,
                         REAL_t const loWn,
                         unsigned int const nF,
                         REAL_t const resolution,
                         int const breadth,
                         unsigned int const numLayers,
                         REAL_t* out_d,
                         const REAL_t* const T_d,
                         const REAL_t* const P_d,
                         const uint8_t* const iso,
                         const REAL_t* const Vnn,
                         const REAL_t* const Snn_ref,
                         const float* const Yair,
                         const float* const Yself,
                         const float* const En,
                         const float* const n,
                         const float* const d,
                         REAL_t* const TauU_d,
                         REAL_t* const pathLength_d,
                         const REAL_t* const PS_d,
                         RefLinePtrs_t* const Lines_d,
                         OpticsBufPtrs_t* const OptBuf_d)
{
    /*Point at the device arrays located in the inputted Lines_d structure.*/
    uint8_t *iso_d = Lines_d->iso;
    REAL_t *Vnn_d = Lines_d->Vnn;
    REAL_t *Snn_ref_d = Lines_d->Snn_ref;
    float *Yair_d = Lines_d->Yair;
    float *Yself_d = Lines_d->Yself;
    float *En_d = Lines_d->En;
    float *n_d = Lines_d->n;
    float *d_d = Lines_d->d;

    /*Point at the device arrays located in the inputted OptBuf_d structure.*/
    REAL_t *Gam_d = OptBuf_d->Gam;
    REAL_t *PShift_d = OptBuf_d->PShift;
    REAL_t *S_d = OptBuf_d->S;

    /*Print out the inputted molecule id.*/
    printf("Hitran molId:=%d\n",
           molId);

    /*Initialize the launch configurator returned block size (dimBlock),
      the minimum grid size needed to achieve the maximum occupancy for a full
      device launch (minGridSize), and the actual grid size needed, based on
      the input data size (dimGrid).*/
    int dimBlock = 0;
    int minGridSize = 0;
    int dimGrid = 0;

#ifdef EVENTS
    /*Record the start time for a cuda event.*/
    cudaEvent_t start;
    cudaEvent_t stop;
    float elapsed;
    HANDLE_ERROR(cudaEventCreate(&start));
    HANDLE_ERROR(cudaEventCreate(&stop));
    HANDLE_ERROR(cudaEventRecord(start,
                                 stream));
#endif

    /*Copy the inputted arrays associated with the Lines_d structure to the
      device.*/
    printf("Memcpy..");
    HANDLE_ERROR(cudaMemcpyAsync(iso_d,
                                 iso,
                                 nL*sizeof(uint8_t),
                                 cudaMemcpyHostToDevice,
                                 stream));
    printf(".done!\n");
    printf("Memcpy..");
    HANDLE_ERROR(cudaMemcpyAsync(Vnn_d,
                                 Vnn,
                                 nL*sizeof(REAL_t),
                                 cudaMemcpyHostToDevice,
                                 stream));
    printf(".done!\n");
    printf("Memcpy..");
    HANDLE_ERROR(cudaMemcpyAsync(Snn_ref_d,
                                 Snn_ref,
                                 nL*sizeof(REAL_t),
                                 cudaMemcpyHostToDevice,
                                 stream));
    printf(".done!\n");
    printf("Memcpy..");
    HANDLE_ERROR(cudaMemcpyAsync(Yair_d,
                                 Yair,
                                 nL*sizeof(float),
                                 cudaMemcpyHostToDevice,
                                 stream));
    printf(".done!\n");
    printf("Memcpy..");
    HANDLE_ERROR(cudaMemcpyAsync(Yself_d,
                                 Yself,
                                 nL*sizeof(float),
                                 cudaMemcpyHostToDevice,
                                 stream));
    printf(".done!\n");
    printf("Memcpy..");
    HANDLE_ERROR(cudaMemcpyAsync(En_d,
                                 En,
                                 nL*sizeof(float),
                                 cudaMemcpyHostToDevice,
                                 stream));
    printf(".done!\n");
    printf("Memcpy..");
    HANDLE_ERROR(cudaMemcpyAsync(n_d,
                                 n,
                                 nL*sizeof(float),
                                 cudaMemcpyHostToDevice,
                                 stream));
    printf(".done!\n");
    printf("Memcpy..");
    HANDLE_ERROR(cudaMemcpyAsync(d_d,
                                 d,
                                 nL*sizeof(float),
                                 cudaMemcpyHostToDevice,
                                 stream));
    printf(".done!\n");

#ifdef EVENTS
    /*Record the stop time for the cuda event and calculate and print the
      elapsed time.*/
    HANDLE_ERROR(cudaEventRecord(stop,
                                 stream));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    HANDLE_ERROR(cudaEventElapsedTime(&elapsed,
                                      start,
                                      stop));
    printf("Time to cudaMemcpy: %3.1f ms\n",
           elapsed);

    /*Record the start time for a cuda event.*/
/*
    HANDLE_ERROR(cudaEventRecord(start,
                                 stream));
*/
#endif

    /*Calculate the thread blocksize and number of thread blocks that maximize
      the occupancy on the device.  Round up to make sure that all input data
      is used.  The CUDA API may produce awarning that can be safely ignored
      depending on the sdk version and -W flags.*/
    HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&minGridSize,
                                                    &dimBlock,
                                                    pre_eval_Snn,
                                                    0,
                                                    ((int)nL)));
    dimGrid = (((int)nL) + dimBlock - 1)/dimBlock;

#ifdef EVENTS
    /*Record the start time for a cuda event.*/
    HANDLE_ERROR(cudaEventRecord(start,
                                 stream));
#endif

    /*Execute the pre_Snn kernel.*/
    printf("%s Kernel Launch.. %d lines ... ",
           "pre_eval_Snn",
           nL);
    printf("using dimBlock = %d and dimGrid = %d  ... ",
           dimBlock,
           dimGrid);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    pre_eval_Snn<<<((unsigned int)dimGrid),((unsigned int)dimBlock),0,stream>>>(nL,
                                                                                molId,
                                                                                iso_d,
                                                                                Vnn_d,
                                                                                En_d,
                                                                                Snn_ref_d);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaPeekAtLastError());
    printf("..peek-OK..");
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    printf("..done!\n");

#ifdef EVENTS
    /*Record the stop time for the cuda event and calculate and print the
      elapsed time.*/
    HANDLE_ERROR(cudaEventRecord(stop,
                                 stream));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    HANDLE_ERROR(cudaEventElapsedTime(&elapsed,
                                      start,
                                      stop));
    printf("Time in Kernel: %3.1f ms\n",
           elapsed);
#endif

    /*Calculate the thread blocksize and number of thread blocks that maximize
      the occupancy on the device.  Round up to make sure that all input data
      is used.*/
    HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&minGridSize,
                                                    &dimBlock,
                                                    eval_gamma,
                                                    0,
                                                    ((int)nL)));
    dimGrid = (((int)nL) + dimBlock - 1)/dimBlock;

#ifdef EVENTS
    /*Record the start time for a cuda event.*/
    HANDLE_ERROR(cudaEventRecord(start,
                                 stream));
#endif

    /*Execute the gamma kernel */
    printf("%s Kernel Launch.. %d lines, %d spatial points ... ",
           "eval_gamma",
           nL,
           numLayers);
    printf("using dimBlock = %d and dimGrid = %d  ... ",
           dimBlock,
           dimGrid);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    eval_gamma<<<((unsigned int)dimGrid),((unsigned int)dimBlock),0,stream>>>(numLayers,
                                                                              nL,
                                                                              P_d,
                                                                              T_d,
                                                                              PS_d,
                                                                              Yself_d,
                                                                              Yair_d,
                                                                              n_d,
                                                                              Gam_d);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaPeekAtLastError());
    printf("..peek-OK..");
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    printf("..done!\n");

#ifdef EVENTS
    /*Record the stop time for the cuda event and calculate and print the
      elapsed time.*/
    HANDLE_ERROR(cudaEventRecord(stop,stream));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    HANDLE_ERROR(cudaEventElapsedTime(&elapsed,
                                      start,
                                      stop));
    printf("Time in Kernel: %3.1f ms\n",
           elapsed);
#endif

    /*Calculate the thread blocksize and number of thread blocks that maximize
      the occupancy on the device.  Round up to make sure that all input data
      is used.*/
    HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&minGridSize,
                                                    &dimBlock,
                                                    eval_pShift,
                                                    0,
                                                    0));
    dimGrid = (((int)nL) + dimBlock - 1)/dimBlock;

#ifdef EVENTS
    /*Record the start time for a cuda event.*/
    HANDLE_ERROR(cudaEventRecord(start,
                                 stream));
#endif

    /*Execute the pshift kernel.*/
    printf("%s Kernel Launch.. %d lines, %d spatial points ... ",
           "eval_pShift",
           nL,
           numLayers);
    printf("using dimBlock = %d and dimGrid = %d  ... ",
           dimBlock,
           dimGrid);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    eval_pShift<<<((unsigned int)dimGrid),((unsigned int)dimBlock),0,stream>>>(numLayers,
                                                                               nL,
                                                                               P_d,
                                                                               Vnn_d,
                                                                               d_d,
                                                                               PShift_d);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaPeekAtLastError());
    printf("..peek-OK..");
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    printf("..done!\n");

#ifdef EVENTS
    /*Record the stop time for the cuda event and calculate and print the
      elapsed time.*/
    HANDLE_ERROR(cudaEventRecord(stop,
                                 stream));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    HANDLE_ERROR(cudaEventElapsedTime(&elapsed,
                                      start,
                                      stop));
    printf("Time in Kernel: %3.1f ms\n",
           elapsed);
#endif

    /*Calculate the thread blocksize and number of thread blocks that maximize
      the occupancy on the device.  Round up to make sure that all input data
      is used.*/
    HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&minGridSize,
                                                    &dimBlock,
                                                    eval_Snn_correction,
                                                    0,
                                                    0));
    dimGrid = (((int)nL) + dimBlock - 1)/dimBlock;

#ifdef EVENTS
    /*Record the start time for a cuda event.*/
    HANDLE_ERROR(cudaEventRecord(start,
                                 stream));
#endif
    /*Execute the Snn kernel.*/
    printf("%s Kernel Launch.. %d lines, %d spatial points ... ",
           "eval_Snn_correction",
           nL,
           numLayers);
    printf("using dimBlock = %d and dimGrid = %d  ... ",
           dimBlock,
           dimGrid);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    eval_Snn_correction<<<((unsigned int)dimGrid),((unsigned int)dimBlock),0,stream>>>(numLayers,
                                                                                       nL,
                                                                                       molId,
                                                                                       T_d,
                                                                                       iso_d,
                                                                                       Vnn_d,
                                                                                       En_d,
                                                                                       Snn_ref_d,
                                                                                       S_d);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaPeekAtLastError());
    printf("..peek-OK..");
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    printf("..done!\n");

#ifdef EVENTS
    /*Record the stop time for the cuda event and calculate and print the
      elapsed time.*/
    HANDLE_ERROR(cudaEventRecord(stop,
                                 stream));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    HANDLE_ERROR(cudaEventElapsedTime(&elapsed,
                                      start,
                                      stop));
    printf("Time in Kernel: %3.1f ms\n",
           elapsed);
#endif

    /*Calculate the thread blocksize and number of thread blocks that maximize
      the occupancy on the device.  Round up to make sure that all input data
      is used.*/
    HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(&minGridSize,
                                                    &dimBlock,
                                                    eval_profile,
                                                    0,
                                                    0));
    dimGrid = (((int)nL) + dimBlock - 1)/dimBlock;

#ifdef EVENTS
    /*Record the start time for a cuda event.*/
    HANDLE_ERROR(cudaEventRecord(start,
                                 stream));
#endif
    /*Execute the profile kernel.*/
    printf("%s Kernel Launch.. %d lines, %d fsamples ... ",
           "eval_profile",
           nL,
           nF);
    printf("using dimBlock = %d and dimGrid = %d  ... ",
           dimBlock,
           dimGrid);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    eval_profile<<<((unsigned int)dimGrid),((unsigned int)dimBlock),0,stream>>>(molId,
                                                                                nL,
                                                                                nF,
                                                                                loWn,
                                                                                resolution,
                                                                                numLayers,
                                                                                breadth,
                                                                                T_d,
                                                                                Gam_d,
                                                                                PShift_d,
                                                                                S_d,
                                                                                pathLength_d,
                                                                                TauU_d,
                                                                                out_d);
#ifdef FORCE_KERNEL_CHECK
    HANDLE_ERROR(cudaPeekAtLastError());
    printf("..peek-OK..");
    HANDLE_ERROR(cudaDeviceSynchronize());
#endif
    printf("..done!\n");

#ifdef EVENTS
    /*Record the stop time for the cuda event and calculate and print the
      elapsed time.*/
    HANDLE_ERROR(cudaEventRecord(stop,
                                 stream));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    HANDLE_ERROR(cudaEventElapsedTime(&elapsed,
                                      start,
                                      stop));
    printf("Time in Kernel: %3.1f ms\n",
           elapsed);

/*
  HANDLE_ERROR( cudaEventRecord(start,stream) );
*/
#endif

#ifdef EVENTS
    /*Clean up cuda event objects.*/
    HANDLE_ERROR(cudaEventDestroy(start));
    HANDLE_ERROR(cudaEventDestroy(stop));
#endif

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
__host__
static void initStreams(unsigned int const numMols,
                        cudaStream_t **streams,
                        int *nStreams)
{
    /*Local variables*/
    int s;

    if (*nStreams == -1)
    {
        if (*streams == NULL)
        {
            cudaDeviceProp p = checkDeviceProps();
            if (p.concurrentKernels)
            {
                *nStreams = (numMols <= MAXNSTREAMS ? numMols : MAXNSTREAMS);
            }
            else
            {
                *nStreams=1;
            }
            fprintf(stderr,
                    "\nUsing %d streams with %d molecules.\n",
                    *nStreams,
                    numMols);
            (*streams) = (cudaStream_t*)malloc(sizeof(cudaStream_t)*(*nStreams));
            if (*nStreams > 1)
            {
                for (s=0;s<(*nStreams);++s)
                {
                    HANDLE_ERROR(cudaStreamCreate(&((*streams)[s])));
                }
            }
            else
            {
                (*streams)[0] = NULL;
            }
        }
    }

    return;
}

/*---------------------------------------------------------------------------*/
__host__
int device_launch(int *nStreams,
                  cudaStream_t **streams,
                  unsigned int const numMols,
                  RefLinePtrs_t L[],
                  REAL_t const loWn,
                  unsigned int const nF,
                  REAL_t const resolution,
                  unsigned int const wingBreadth,
                  radiationOutputFields_t *atmosData,
                  REAL_t * const out,
                  int const time,
                  int const lat,
                  int const lon)
{
    /*Local variables*/
    unsigned int m;
    unsigned int mol;
    unsigned int numLayers=atmosData->npfull;
    int s;
    REAL_t *T_d;
    REAL_t *P_d;
    REAL_t *N_d;
    REAL_t *Z_d;
    REAL_t *PS_d;
    REAL_t *out_d;

    /*Initialize TIPS.*/
    initTIPS_d();
    initStreams(numMols,
                streams,
                nStreams);

    /*Sanity check.*/
    assert(numMols<=NUM_MOL);

    /*Point the beginning of the appropriate "column" of data.*/
    const size_t idx = (time*atmosData->nlon*atmosData->nlat +
                        lat*atmosData->nlon + lon)*atmosData->npfull;
    REAL_t *T = &(atmosData->T[idx]);
    REAL_t *P = &(atmosData->P[idx]);
    REAL_t *DELTAZ = &(atmosData->DELTAZ[idx]);
    REAL_t *N = &(atmosData->N[idx*NUM_MOL]);
    REAL_t *PS = &(atmosData->PS[idx*NUM_MOL]);

    /*Print out values for debugging.  Delete this later.*/
/*
    if (time == 0 && lat == 0 && lon == 0)
    {
        printf("time = %d, lat = %d, lon = %d, mol = %d\n",
               time,
               lat,
               lon,
               0);
        printf("layer n(cm^-3)   delz(cm)   delz*n(cm^-2)   P(atm)"
                   "   Ps(atm)\n");
        for (m=0;m<numLayers;m++)
        {
            printf("%u %e %e %e %e %e\n",
                   m,
                   N[m+atmosData->npfull],
                   DELTAZ[m],
                   N[m+atmosData->npfull]*DELTAZ[m],
                   P[m],
                   PS[m+atmosData->npfull]);
        }
    }
*/

    /*Malloc and copy data to the device arrays.*/
    device_atmos_init(numLayers,
                      NUM_MOL,
                      T,
                      P,
                      N,
                      DELTAZ,
                      PS,
                      &T_d,
                      &P_d,
                      &N_d,
                      &Z_d,
                      &PS_d);

    /*Malloc the out_d array and set its values to all zeros by memcpying
      out_h to out_d.  For this to work correctly, out_h should be all zeros.*/
/*
    HANDLE_ERROR(cudaMalloc(&out_d,
                            (numLayers+1)*nF*sizeof(REAL_t)));
*/
    HANDLE_ERROR(cudaMalloc(&out_d,
                            (numLayers)*nF*sizeof(REAL_t)));
    printf("Memcpy.. out_h -> out_d");
    printf("Memcpy.. out_h -> out_d");
    HANDLE_ERROR(cudaMemcpy(out_d,
                            out,
                            numLayers*nF*sizeof(REAL_t),
                            cudaMemcpyHostToDevice));
    printf(".done!\n");

    /*Setup buffers for re-use under streaming.*/
    OpticsBufPtrs_t OpticsBuf_d[*nStreams];
    RefLinePtrs_t LinesBuf_d[*nStreams];
    for (s=0;s<*nStreams;++s)
    {
        device_optics_init(numLayers,
                           MAX_NUM_SPECTRAL_LINES,
                           &(LinesBuf_d[s]),
                           &(OpticsBuf_d[s]));
    }

    /*Loop throught the molecules and calculate the optical depth values.*/
    for (m=0;m<(numMols+(*nStreams-1));m+=*nStreams)
    {
        for (s=0;s<*nStreams;++s)
        {
            mol = m+s;
            if (mol<numMols)
            {
                device_optics_perMol((*streams)[s],
                                     L[mol].mol,
                                     L[mol].nLines,
                                     loWn,
                                     nF,
                                     resolution,
                                     wingBreadth,
                                     numLayers,
                                     out_d,
                                     T_d,
                                     P_d,
                                     L[mol].iso,
                                     L[mol].Vnn,
                                     L[mol].Snn_ref,
                                     L[mol].Yair,
                                     L[mol].Yself,
                                     L[mol].En,
                                     L[mol].n,
                                     L[mol].d,
                                     &(N_d[(L[mol].mol-1)*atmosData->npfull]),
                                     Z_d,
                                     &(PS_d[(L[mol].mol-1)*atmosData->npfull]),
                                     &(LinesBuf_d[s]),
                                     &(OpticsBuf_d[s]));
            }
        }
    }

    /*Synchronize all streams.*/
    HANDLE_ERROR(cudaDeviceSynchronize());

    /*Free device buffers.*/
    for (s=0;s<*nStreams;++s)
    {
        device_optics_free(&(LinesBuf_d[s]),
                           &(OpticsBuf_d[s]));
    }

    /*Memcpy the optical depths from the device back to the host.*/
    printf("Memcpy opticalDepth (out) to host..");
    HANDLE_ERROR(cudaMemcpy(out,
                            out_d,
                            numLayers*nF*sizeof(REAL_t),
                            cudaMemcpyDeviceToHost));
    printf("..done\n");

    /*Free device arrays.*/
    device_atmos_free(T_d,
                      P_d,
                      N_d,
                      Z_d,
                      PS_d);
    HANDLE_ERROR(cudaFree(out_d));

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/

#else
/*---------------------------------------------------------------------------*/
int device_launch(int* nStreams,
                  void** streams,
                  const unsigned int numMols,
                  const char* const molFnames[],
                  const REAL_t loWn,
                  const unsigned int nF,
                  const REAL_t resolution,
                  const unsigned int wingBreadth,
                  radiationOutputFields_t* atmosData,
                  REAL_t* const out,
                  int time,
                  int lat,
                  int lon)
{
    /*Prevent compiler warnings.*/
    (void) nStreams;
    (void) streams;
    (void) numMols;
    (void) molFnames;
    (void) loWn;
    (void) nF;
    (void) resolution;
    (void) wingBreadth;
    (void) atmosData;
    (void) out;
    (void) time;
    (void) lat;
    (void) lon;

    printf("\n\nYou've not compiled with NVCC, device_launch does"
               " nothing...\n\n");

    return -1;
}

/*---------------------------------------------------------------------------*/
#endif

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Main part of the program.*/
int main(int argc,
         char* argv[])
{
    /*Local variables*/
    int world_size = -1; /*Number of ranks in MPI_COMM_WORLD.*/
    int world_rank = -1; /*Process rank id in MPI_COMM_WORLD.*/
#ifdef MPI_ENABLED
    int ierr = 0;        /*MPI error code.*/
#endif

    int time;
    unsigned int lat;
    unsigned int compute_lat_beg;
    unsigned int compute_lat_end;
    unsigned int lon;
    unsigned int compute_lon_beg;
    unsigned int compute_lon_end;
    unsigned int mol;
    size_t idx;
    REAL_t *out = NULL;

    struct arguments arguments;
    static char default_output_fname[] ="didyouforgettospecifyoutfile.nc";

    /*If necessary, initialize MPI.*/
#ifdef MPI_ENABLED
    ierr = MPI_Init(&argc,
                    &argv);
    if (ierr != MPI_SUCCESS)
    {
        fprintf(stderr,
                "Error(main): MPI_Init returned error code: %d.\n",
                ierr);
        exit(EXIT_FAILURE);
    }
    ierr = MPI_Comm_size(MPI_COMM_WORLD,
                         &world_size);
    if (ierr != MPI_SUCCESS)
    {
        fprintf(stderr,
                "Error(main): MPI_Comm_size returned error code: %d.\n",
                ierr);
        exit(EXIT_FAILURE);
    }
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,
                         &world_rank);
    if (ierr != MPI_SUCCESS)
    {
        fprintf(stderr,
                "Error(main): MPI_Comm_rank returned error code: %d.\n",
                ierr);
        exit(EXIT_FAILURE);
    }
#endif

    /*Set default argument values.*/
    arguments.silent = 0;
    arguments.verbose = 0;
    arguments.device = 0;
    arguments.mpi = 0;
    arguments.host = 0;
    arguments.nhitfiles = 0;
    arguments.nmolConc = 0;
    arguments.nmolConcOver = 0;
    arguments.atmos = NULL;
    arguments.atmos_format = NULL;
    arguments.output_file = default_output_fname;
    arguments.wingBreadth = 25;
    arguments.ctm = 0;
    arguments.w = 1;
    arguments.W = 3000;
    arguments.t = 0;
    arguments.T = 0;
    arguments.res = 1.0;
    arguments.h2o = 0;
    arguments.co2 = 0;
    arguments.o3 = 0;
    arguments.n2o = 0;
    arguments.co = 0;
    arguments.ch4 = 0;
    arguments.o2 = 0;

    /*Parse the program's arguments using arg_parse.*/
    argp_parse(&argp,
               argc,
               argv,
               0,
               0,
               &arguments);

    /*If a water concentration was not specified in the program's arguments,
      then use the value from the inputted netCDF file if it exists.*/
    if (arguments.h2o == 0)
    {
        arguments.h2o = -1;
        arguments.nmolConcOver++;
    }

    /*If a carbon dioxide concentration was not specified in the program's
      arguments, then use the value from the inputted netCDF file if it
      exists.*/
    if (arguments.co2 == 0)
    {
        arguments.co2 = -1;
        arguments.nmolConcOver++;
    }

    /*If a ozone concentration was not specified in the program's arguments,
      then use the value from the inputted netCDF file if it exists.*/
    if (arguments.o3 == 0)
    {
        arguments.o3 = -1;
        arguments.nmolConcOver++;
    }

    /*If a nitrous oxide concentration was not specified in the program's
      arguments, then use the value from the inputted netCDF file if it
      exists.*/
    if (arguments.n2o == 0)
    {
        arguments.n2o = -1;
        arguments.nmolConcOver++;
    }

    /*If a carbon monoxide concentration was not specified in the program's
      arguments, then use the value from the inputted netCDF file if it
      exists.*/
    if (arguments.co == 0)
    {
        arguments.co = -1;
        arguments.nmolConcOver++;
    }

    /*If a methane concentration was not specified in the program's
      arguments, then use the value from the inputted netCDF file if it
      exists.*/
    if (arguments.ch4 == 0)
    {
        arguments.ch4 = -1;
        arguments.nmolConcOver++;
    }

    /*If an oxygen concentration was not specified in the program's
      arguments, then use the value from the inputted netCDF file if it
      exists.*/
    if (arguments.o2 == 0)
    {
        arguments.o2 = -1;
        arguments.nmolConcOver++;
    }

    /*Set the wavenumber "grid" size (i.e., the number of different wavenumber
      points at which the spectra will be calculated).*/
    const unsigned int nF = (arguments.W-arguments.w)/arguments.res + 1;

    /*Make sure either only host or device has been targeted.*/
    if (arguments.device != 0 && arguments.host != 0)
    {
        fprintf(stderr,
                "Error(main): more than one target specified.  Please use"
                    " either --host or --device or neither flag to just"
                    " default to device 0.\n");
        exit(EXIT_FAILURE);
    }
    else if (arguments.mpi != 0 && (arguments.device == 0 &&
             arguments.host == 0))
    {
        fprintf(stderr,
                "Error(main): when specifying mpi you must specify the number"
                    " of devices per node or --host.\n");
        exit(EXIT_FAILURE);
    }

    /*Set launchType = 0 for host, = 1 for device.*/
    const int launchType = arguments.host == 1 ? 0 : 1;

#ifdef MPI_ENABLED
    /*Determine the number of devices.*/
    const int device_number = world_size % arguments.device;

    assert(device_number >= 0); /* if this ever trips try (a+n) % n */

    /*Make sure that mpirun is used to execute the program if mpi is turned
      on.*/
    if (arguments.mpi != 0 && world_size < 1)
    {
        fprintf(stderr,
                "Error(main): you have invoked with the program with --mpi"
                    " but did not used an mpirun style executer.\n");
        exit(EXIT_FAILURE);
    }
#else
    /*Make sure that MPI_ENABLED was included if MPI is turned on.*/
    if (arguments.mpi != 0)
    {
        fprintf(stderr,
                "Error(main): you must build with -DMPI_ENABLED in order to"
                    " use MPI.\n");
        exit(EXIT_SUCCESS);
    }
#endif

    /*Set the GPU device number.*/
    if (launchType == 1)
    {
#ifdef __NVCC__
        const int device_number = arguments.device;
        HANDLE_ERROR(cudaSetDevice(device_number));
#endif
    }

    /*Read in atmospheric data from the inputted netCDF file.*/
    char *atmosFile = arguments.atmos;
    char *atmosFormat = arguments.atmos_format;
    radiationOutputFields_t atmosData;
    getAndSetAtmosFieldsFromFile(atmosFile,
                                 atmosFormat,
                                 &atmosData);
    const unsigned int numLayers = atmosData.npfull;

    /*Check to make sure that the number of inputted molecular
      concentrations matches the number of inputted HITRAN files.*/
    if (arguments.nmolConc != arguments.nhitfiles)
    {
        fprintf(stderr,
                "Warning(main): the number of hitfiles (%d) does not match"
                    " the number of prescribed concentrations (%d). Checking"
                    " for overrides...\n",
                arguments.nhitfiles,
                arguments.nmolConc);
        if (arguments.nmolConc+arguments.nmolConcOver == arguments.nhitfiles)
        {
            fprintf(stderr,
                    "\t...found %d overrides, okay.\n",
                    arguments.nmolConcOver);
        }
        else
        {
            fprintf(stderr,
                    "\t...found %d overrides.\nError(main): the number of"
                        " inputted hitfiles does not match the number of"
                        " inputted + overridden molecular concentrations.\n",
                    arguments.nmolConcOver);
            exit(EXIT_FAILURE);
        }
    }
    const unsigned int nMols = arguments.nhitfiles;
    char** hitFnameList = arguments.hitfiles;
    printf("\nSubmitted %d molecules.\n",nMols);

    /*Initialize the output file.*/
    int ncid;
    int varid;
    char *OUTPUT_FNAME=NULL;
    compute_lat_beg = 0;
    compute_lat_end = atmosData.nlat;
    compute_lon_beg = 0;
    compute_lon_end = atmosData.nlon;

    /* output file and compute setup is very different for mpi */
    if (arguments.mpi != 0)
    {
        OUTPUT_FNAME = (char *)malloc(strlen(arguments.output_file) + 9);
        if (OUTPUT_FNAME == NULL)
        {
            fprintf(stderr,
                    "Error(main): malloc failed for %zu bytes of"
                        " OUTPUTFNAME.\n",
                    strlen(arguments.output_file) + 9);
            exit(EXIT_FAILURE);
        }
        lat = atmosData.nlat / world_size;
        if(lat*world_size != atmosData.nlat)
        {
            fprintf(stderr, 
                    "Warning(main): specified %d global lats across ranks=%zu"
                        " yields between %d and %d lats per rank.  This will"
                        " result in idle hardware, suggest a different work"
                        " share.\n",
                    world_size,
                    atmosData.nlat,
                    lat,
                    lat+1);
        }
        compute_lat_beg = world_rank*lat;
        compute_lat_end = compute_lat_beg+lat;
        if (compute_lat_end>atmosData.nlat)
        {
            compute_lat_end = atmosData.nlat;
        }
        compute_lon_beg = 0;
        compute_lon_end = atmosData.nlon;
    /* copy existing name
     * cat .rankN */
        sprintf(OUTPUT_FNAME,
                "%s.rank%d",
                arguments.output_file,
                world_rank);
    }
    else
    {
        OUTPUT_FNAME = arguments.output_file;
    }
    fprintf(stderr,
            "Opening output file %s.\n",
            OUTPUT_FNAME);
    openOpticalDepthOutput(&ncid,
                           &varid,
                           OUTPUT_FNAME,
                           compute_lat_end - compute_lat_beg,
                           compute_lon_end - compute_lon_beg,
                           numLayers,
                           nF);

    /*Declare stream parameters.*/
#ifdef __NVCC__
    int nstreams = -1;
    cudaStream_t* streams = NULL;
#endif

    /*Setup HITRAN lines.*/
    RefLinePtrs_t HitLines[nMols];

    /*Get filename and parse in lines*/
    RefLine_flags_t flags= {((unsigned int) -1),1,0}; /* host cuda malloc default, host=True, device=false */
    arguments.T = atmosData.ntime;
    time = arguments.T;
/*
    time = arguments.T - arguments.t + 1;
*/
    for(mol=0;mol<nMols;++mol)
    {
        HitLines[mol] = parseHITRANfile(hitFnameList[mol],
                                        flags,
                                        arguments.w,
                                        arguments.W);

        /*Check the molecular configurations.  For all molecules whose
          partial pressure is not taken from the input NetCDF file, calculate
          the partial pressure from the concentrations inputted on the
          command line.*/
        checkMolConfig(&arguments,
                       HitLines[mol].mol,
                       &atmosData,
                       time);
    }

    /*Compute the spectra.*/
    for (time=arguments.t;time<arguments.T;++time)
    {
        for (lat=compute_lat_beg;lat<compute_lat_end;++lat)
        {
            for (lon=compute_lon_beg;lon<compute_lon_end;++lon)
            {
                if (launchType == 0)
                {
                    if (out == NULL)  /* malloc if needed, otherwise pass */
                    {
                        out = (REAL_t*)calloc(nF*numLayers,
                                              sizeof(REAL_t));
                    }
                    host_launch(nMols,
                                HitLines,
                                ((REAL_t)arguments.w),
                                nF,
                                arguments.res,
                                arguments.wingBreadth,
                                &atmosData,
                                out,
                                time,
                                lat,
                                lon);
                }
                else if (launchType == 1)
                {
#ifdef __NVCC__
                    if(out == NULL)  /* malloc if needed, otherwise pass */
                    {
                        HANDLE_ERROR(cudaHostAlloc(&out,
                                                   nF*numLayers*sizeof(REAL_t),
                                                   cudaHostAllocDefault));
                    }
                    device_launch(&nstreams,
                                  &streams,
                                  nMols,
                                  HitLines,
                                  ((REAL_t)arguments.w),
                                  nF,
                                  arguments.res,
                                  arguments.wingBreadth,
                                  &atmosData,
                                  out,
                                  time,
                                  lat,
                                  lon);
#else
                    fprintf(stderr,
                            "Error(main): requested cuda launch type (%d),"
                                " but compiled host only.\n",
                            launchType);
                    exit(EXIT_FAILURE);
#endif
                }
                else
                {
                    fprintf(stderr,
                            "Error(main): unkown launch type (%d)"
                                " requested.\n",
                            launchType);
                    exit(EXIT_FAILURE);
                }

                /*Calculate the continuum spectra.*/
                if (arguments.ctm == 1)
                {
                    printf("Computing Continuum\n");
                    assert(arguments.res==1.);  /* presently the continuum code is only safe for widths of one wavenumber */
                    idx = (time*atmosData.nlon*atmosData.nlat + 
                           lat*atmosData.nlon + lon)*atmosData.npfull;
                    /* dbg print */
                    printf("%zu: T=%g P=%g DELTAZ=%g PS[%zu]=%g \n",
                           idx,
                           atmosData.T[idx],
                           atmosData.P[idx],
                           atmosData.DELTAZ[idx],
                           NUM_MOL*idx + H2O*atmosData.npfull,
                           atmosData.PS[NUM_MOL*idx + H2O*atmosData.npfull]);
                    /* getchar(); */
                    get_CTM(out,
                            &(atmosData.T[idx]),
                            &(atmosData.P[idx]),
                            &(atmosData.DELTAZ[idx]),
                            &(atmosData.PS[NUM_MOL*idx + H2O*atmosData.npfull ]),
                            nF,
                            atmosData.npfull);
                }

                fprintf(stderr,
                        "Writing hyperslab of %d samples "
                        "@{t=%d, lat=%d, lon=%d, layers=0:%d} to output"
                        " file %s \n",
                        nF,
                        time,
                        lat,
                        lon,
                        numLayers,
                        OUTPUT_FNAME);
                writeOpticalDepthOutputByColumn(ncid,
                                                varid,
                                                time,
                                                lat-compute_lat_beg,
                                                lon-compute_lon_beg,
                                                numLayers,
                                                nF,
                                                out);

                memset(out, 0, nF*numLayers*sizeof(REAL_t));
            }  /* nlat */
        }  /* nlon */
    }  /* time */

    /*Close the output file.*/
    fprintf(stderr,
            "Closing output file %s\n",
            OUTPUT_FNAME);
    closeOpticalDepthOutput(ncid);

    /*Cleanup all, this is dirty,  into earlier stage later */
    for (mol=0;mol<nMols;++mol)
    {
        freeHost(HitLines[mol],flags);
    }

///wrap up something like this in a function, then key off of outputfile extension for csv output
/* #undef WRITEOUT */
/*   ///#define WRITEOUT */
/* #ifdef WRITEOUT */
/*   printf("\n\n\t Attempting Result Write.\n\n"); */
/*   /\* output *\/ */
/*   REAL_t wv; */
/*   FILE* ofp; */
/*   ofp = fopen("opticaldepth.testout.csv","w"); */
/*   if(ofp==NULL){ */
/*     fprintf(stderr,"\nopening output file for writing failed, aborting.\n"); */
/*     exit(1); */
/*   } */
/*   /\* header *\/ */
/*   fprintf(ofp, "z,wavenumber,val\n"); */
/*   /\* data *\/ */
/*   for(unsigned int l=0; l<numLayers; ++l){ */
/*     /\* for(unsigned int iter=0; iter<nF; iter+=(((double)1)/arguments.res) ){ /\\* output every one wavenumber *\\/ *\/ */
/*     for(unsigned int iter=0; iter<nF; ++iter ){ /\* output every sample *\/ */
/*       wv = iter*arguments.res + arguments.w; */
/*       assert(wv <= arguments.W); */
/*       fprintf(ofp,"%u %.17f %.17f\n", l, wv , out[ l*nF + iter ] ); */
/*     } */
/*   } */
  
/*   if( fclose(ofp)!=0 ){ */
/*     fprintf(stderr,"\nclosing output file for writing failed, aborting.\n"); */
/*     exit(1); */
/*   } */
/* #endif */

    /* cleanup */
    if (launchType==1)
    {
#ifdef __NVCC__
        cudaFreeHost(out);
#endif
    }
    else
    {
        free(out);
    }

#ifdef MPI_ENABLED
    MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}

