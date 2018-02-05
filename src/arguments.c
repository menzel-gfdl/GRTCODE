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
     "maximum time (upper bound, inclusive), defaults to maximum time level"
         " in the input file.",
     -3},

    {"minlon",
     'x',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "minimum longitude index (lower bound, inclusive), defaults 0",
     -3},

    {"maxlon",
     'X',
     "VAL",
     OPTION_ARG_OPTIONAL, 
     "maximum longitude index (upper bound, inclusive), defaults to"
         " maximum longitude index in the input file.",
     -3},

    {"minlat",
     'y',
     "VAL",
     OPTION_ARG_OPTIONAL,
     "minimum latitude index (lower bound, inclusive), defaults 0",
     -3},

    {"maxlat",
     'Y',
     "VAL",
     OPTION_ARG_OPTIONAL, 
     "maximum latitude index (upper bound, inclusive), defaults to"
         " maximum latitude index in the input file.",
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
    int x; /*Starting longitude index, inclusive.*/
    int X; /*Ending longitude index, inclusive.*/
    int y; /*Starting latitude index, inclusive.*/
    int Y; /*Ending longitude index, inclusive.*/
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
