#include <argp.h>
#include <ctype.h>
#include "arguments.h"
#include "constants.h"
#include "debug.h"
#include "molecules.h"
#include "utils.h"

/*Variables required by argparse.*/
char const *argp_program_version = "lbl-dev 0.1";
char const *argp_program_bug_address = "<raymond.menzel@noaa.gov>";
static char doc[] = "GFDL style documentation goes >/\n\n\\"
                    "<^here.\n\vOther Documentation goes here.";

static unsigned int const minNargs = minNhitfiles;
static unsigned int const maxNargs = maxNhitfiles;
static char args_doc[] = "-aINPUT.nc -oOUT.nc"
                         " [molecule concentration specifications]"
                         " HITFILES";

enum arg_group_types
{
    LAUNCH_GROUP,
    IO_GROUP,
    BOUNDS_GROUP,
    PPMV_GROUP
};

static struct argp_option options[] =
{
    {"device",
     'd',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Use gpu implementation on specifed DEVICE."
         "\n\tDefault DEVICE is simply GPU0"
         "\n\tIncompatible with --host."
         "\n\t --mpi modifies this flag to prescribe numDevices per node.",
     LAUNCH_GROUP},

    {"host",
     'h',
     0,
     0,
     "Use HOST cpu implementation. \n\t(incompatible with --device)",
     LAUNCH_GROUP},

    {"output",
     'o',
     "FILE",
     0,
     "Output to FILE",
     IO_GROUP},

    {"atmos",
     'a',
     "INPUT.NC",
     0,
     "NC file containing model atmosphere.",
     IO_GROUP},

    {"minw",
     'w',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "minimum Wavenumber (lower bound, inclusive), defaults 1",
     BOUNDS_GROUP},

    {"maxw",
     'W',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "maximum Wavenumber (upper bound, inclusive), defaults 50000",
     BOUNDS_GROUP},

    {"mint",
     't',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "minimum time (lower bound, inclusive), defaults 0",
     BOUNDS_GROUP},

    {"maxt",
     'T',
     "VAL",
     OPTION_ARG_OPTIONAL, 
     "maximum time (upper bound, inclusive), defaults to maximum time level"
         " in the input file.",
     BOUNDS_GROUP},

    {"minlon",
     'x',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "minimum longitude index (lower bound, inclusive), defaults 0",
     BOUNDS_GROUP},

    {"maxlon",
     'X',
     "VAL",
     OPTION_ARG_OPTIONAL, 
     "maximum longitude index (upper bound, inclusive), defaults to"
         " maximum longitude index in the input file.",
     BOUNDS_GROUP},

    {"minlat",
     'y',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "minimum latitude index (lower bound, inclusive), defaults 0",
     BOUNDS_GROUP},

    {"maxlat",
     'Y',
     "VAL",
     OPTION_ARG_OPTIONAL, 
     "maximum latitude index (upper bound, inclusive), defaults to"
         " maximum latitude index in the input file.",
     BOUNDS_GROUP},

    {"res",
     'r',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Resolution (wavenumber), defaults 1.0",
     BOUNDS_GROUP},

    {"wings",
     'c',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Wings cutoff (+/- integer wavenumber), defaults 25",
     BOUNDS_GROUP},

    {"h2o",
     '1',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Water Concentration.  Default reads from INPUT.NC, else supply global"
         " value (ppmv).  Layer partial pressure = (layer pressure)*"
         "(water concentration/10^6).",
     PPMV_GROUP},

    {"co2",
     '2',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Carbon Dioxide Concentration.  Default reads from INPUT.NC, else"
         " supply global value (ppmv).  Layer partial pressure ="
         " (layer pressure)*(carbon dioxide concentration/10^6).",
     PPMV_GROUP},

    {"o3",
     '3',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Ozone Concentration.  Default reads from INPUT.NC, else supply global"
         " value (ppmv).  Layer partial pressure = (layer pressure)*"
         "(ozone concentration/10^6).",
     PPMV_GROUP},

    {"n2o",
     '4',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Nitrous Oxide.  Default reads from INPUT.NC, else supply global value"
         " (ppmv).  Layer partial pressure = (layer pressure)*"
         "(nitrous oxide concentration/10^6).",
     PPMV_GROUP},

    {"co",
     '5',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Carbon Monoxide Concentration.  Default read from INPUT.NC, else"
         " supply global value (ppmv).  Layer partial pressure ="
         " (layer pressure)*(carbon monoxide concentration/10^6).",
     PPMV_GROUP},

    {"ch4",
     '6',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Methane Conentration.  Default read from INPUT.NC, else supply global"
         " value (ppmv).  Layer partial pressure = (layer pressure)*"
         "(methane concentration/10^6).",
     PPMV_GROUP},

    {"o2",
     '7',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "Oxygen Concentration.  Default read from INPUT.nc, else supply global"
         " value (ppmv).  Layer partial pressure = (layer_pressure)*"
         "(oxygen concentration/10^6).",
     PPMV_GROUP},

    {"ctm",
     'C',
     0,
     0,
     "Enables the water vapor continuum.",
     PPMV_GROUP},

    {0}
};

/*Helper parsing function for molecular concentrations.*/
static int parse_mol_conc(char *arg,
                          double *res)
{
    not_null(arg);
    if (arg[0] == 'a')
    {
        /*A leading 'a' character specifies that the concentration should be
         taken from the input netcdf file.*/
        *res = CONC_FROM_FILE;
    }
    else if (isalpha(arg[0]))
    {
        fatal("the supplied character (%c) for overriding a molecule"
                  " concentration is not understood.  Review args.",
              arg[0]);
    }
    else
    {
        check(to_double(arg,
                        res));
        if (*res <= 0.0)
        {
            fatal("the supplied molecular concentration (%e) must be > 0.",
                  *res);
        }
    }
    return SUCCESS;
}

#define check_usage(e) {if (e != SUCCESS) {argp_usage(state);}}

static error_t parse_opt(int key,
                         char *arg,
                         struct argp_state *state)
{
    not_null(state);
    struct arguments *arguments = (struct arguments*)(state->input);

    /*Store the inputted arguments.*/
    switch(key)
    {
        case 'a':
            arguments->atmosInputFile = arg;
            break;
        case 'c':
            check_usage(to_int(arg,
                               &(arguments->wingBreadth)));
            break;
        case 'C':
            arguments->ctm = 1;
            break;
        case 'd':
            check_usage(to_int(arg,
                               &(arguments->device)));
            break;
        case 'h':
            arguments->host = 1;
            break;
        case 'o':
            arguments->outputFile = arg;
            break;
        case 'r':
            check_usage(to_double(arg,
                                  &(arguments->res)));
            break;
        case 't':
            check_usage(to_int(arg,
                               &(arguments->t)));
            break;
        case 'T':
            check_usage(to_int(arg,
                               &(arguments->T)));
            break;
        case 'w':
            check_usage(to_int(arg,
                               &(arguments->w)));
            break;
        case 'W':
            check_usage(to_int(arg,
                               &(arguments->W)));
            break;
        case 'x':
            check_usage(to_int(arg,
                               &(arguments->x)));
            break;
        case 'X':
            check_usage(to_int(arg,
                               &(arguments->X)));
            break;
        case 'y':
            check_usage(to_int(arg,
                               &(arguments->y)));
            break;
        case 'Y':
            check_usage(to_int(arg,
                               &(arguments->Y)));
            break;
        case '1':
            check_usage(parse_mol_conc(arg,
                                       &(arguments->molConc[H2O])));
            break;
        case '2':
            check_usage(parse_mol_conc(arg,
                                       &(arguments->molConc[CO2])));
            break;
        case '3':
            check_usage(parse_mol_conc(arg,
                                       &(arguments->molConc[O3])));
            break;
        case '4':
            check_usage(parse_mol_conc(arg,
                                       &(arguments->molConc[N2O])));
            break;
        case '5':
            check_usage(parse_mol_conc(arg,
                                       &(arguments->molConc[CO])));
            break;
        case '6':
            check_usage(parse_mol_conc(arg,
                                       &(arguments->molConc[CH4])));
            break;
        case '7':
            check_usage(parse_mol_conc(arg,
                                       &(arguments->molConc[O2])));
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num >= maxNargs)
            {
                log_err("there are too many (%d) command line arguments"
                            " (only %d allowed).",
                        state->arg_num,
                        maxNargs);
                argp_usage(state);
            }
            arguments->hitFiles[state->arg_num] = arg;
            arguments->nHitFiles++;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < minNargs)
            {
                log_err("there are too few (%d) command line arguments"
                            " (%d are required.).",
                        state->arg_num,
                        minNargs);
                argp_usage(state);
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options,parse_opt,args_doc,doc,NULL,NULL,NULL};

void parse_options(int argc,
                   char **argv,
                   struct arguments *arguments)
{
    argp_parse(&argp,
               argc,
               argv,
               0,
               0,
               arguments);
}
