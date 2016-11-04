#include <argp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "parseOpticalDepthOutput.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Argp variables and functions.*/

/*---------------------------------------------------------------------------*/
/*Set argp variables.*/

const char *argp_program_version = "verification-dev 0.0";
const char *argp_program_bug_address = "<raymond.menzel@noaa.gov>";
static char doc[] = "Comments go here";
static char args_doc[] = "-rREFERENCE.spc -oOUTPUT GRTCODE_OUTPUT";

/*Command line options.*/
static struct argp_option options[] =
{
    {"reference",
     'r',
     "REFERENCE.spc",
     0,
     "Text file containing RFM reference data."},

    {"output",
     'o',
     "OUTPUT",
     0,
     "File where output will be written."},

    {0}
};

struct arguments
{
    char *reference_file; /*Input RFM reference file.*/
    char *output_file;    /*Output text file.*/
    char *grtcode_file;   /*NetCDF file produced by the GRTcode.*/
};

/*---------------------------------------------------------------------------*/
/*Argp options parser function.*/
static error_t parse_opt(int key,
                         char *arg,
                         struct argp_state *state)
{
    /*Local variables*/
    struct arguments *arguments = NULL; /*Arguments pointer.*/

    /*Point to the inputted argument from argp_parse.*/
    arguments = (struct arguments*)(state->input);

    /*Store the inputted arguments.*/
    switch(key)
    {
        case 'r':
            arguments->reference_file = arg;
            break;
        case 'o':
            arguments->output_file = arg;
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num > 1)
            {
                fprintf(stderr,
                        "Error(parse_opt): there are too many command line"
                        " arguments.\n");
                argp_usage(state);
            }
            arguments->grtcode_file = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 1)
            {
                fprintf(stderr,
                        "Error(parse_opt): there are too few command line"
                        " arguments.\n");
                argp_usage(state);
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

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*Main part of the program.*/
int main(int argc,
         char* argv[])
{
    /*Local variables*/
    struct arguments arguments;                 /*Arg arguments structure.*/
    char *grt_filename;                         /*File containing GRTcode output.*/
    OpticalDepth_t grt_optical_depth;           /*Optical depth values from GRTcode.*/
    char *rfm_filename;                         /*File containing RFM output.*/
    OpticalDepth_t rfm_optical_depth;           /*Optical depth values from RFM.*/
    size_t time_end;                            /*Number of times at which the data will be compared.*/
    size_t lat_end;                             /*Number of latitudes at which the data will be compared.*/
    size_t lon_end;                             /*Number of longitudes at which the data will be compared.*/
    size_t pfull_end;                           /*Number of pressure layers at which the data will be compared.*/
    size_t wavenumber_end;                      /*Number of wavenumbers at which the data will be compared.*/
    size_t wavenumber_start;                    /*Number of wavenumbers at which the data will be compared.*/
    int grt_shift;                              /*Shift in the GRT index needed because wavenumbers start at 1.*/
    int rfm_shift;                              /*Shift in the RFM index needed because wavenumbers start at 0.*/
    size_t grt_index = 0;                       /*GRT optical depth index variable for comparison.*/
    size_t rfm_index = 0;                       /*RFM optical depth index variable for comparison.*/
    float diff;                                 /*Difference between the optical depth values.*/
    float percent_error;                        /*Precent error for the GRT optical depth compared to RFM.*/
    float max_error = -1.0;                     /*Maximum percent error.*/
    size_t time;                                /*Loop variable.*/
    size_t lat;                                 /*Loop variable.*/
    size_t lon;                                 /*Loop variable.*/
    size_t pfull;                               /*Loop variable.*/
    size_t wavenumber;                          /*Loop variable.*/
    char *output_file_name;                     /*Output file name.*/
    FILE *output_file = NULL;                   /*Output file pointer.*/
    char *output_file_name_gnu_plot;            /*Output file name for gnu_plot heat map.*/
    FILE *output_file_gnu_plot = NULL;          /*Output file pointer for gnu_plot heat map.*/
    int ioerr;                                  /*I/O error code.*/
    size_t num_rfm_opt_depth_non_zero;          /*Number of RFM optical depth values that are not zero.*/
    size_t num_errors_over_one_percent;         /*Number of optical depth differences that are greater than or eqaul to one percent.*/
    size_t num_errors_over_ten_percent;         /*Number of optical depth differences that are greater than or eqaul to ten percent.*/
    size_t num_errors_over_fifty_percent;       /*Number of optical depth differences that are greater than or eqaul to fifty percent.*/
    size_t num_errors_over_one_hundred_percent; /*Number of optical depth differences that are greater than or eqaul to one hundred percent.*/

    /*Set default argument values.*/
    arguments.reference_file = NULL;
    arguments.output_file = NULL;
    arguments.grtcode_file = NULL;

    /*Parse the program's arguments using arg_parse.*/
    argp_parse(&argp,
               argc,
               argv,
               0,
               0,
               &arguments);

    /*Make sure that the arguments are not null.*/
    if (arguments.reference_file == NULL)
    {
        fprintf(stderr,
                "Error(main): the RFM reference file was not passed in"
                "correctly.\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        rfm_filename = arguments.reference_file;
    }
    if (arguments.output_file == NULL)
    {
        fprintf(stderr,
                "Error(main): the output file was not passed in"
                "correctly.\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        output_file_name = arguments.output_file;
        output_file_name_gnu_plot = (char *)malloc(sizeof(char)*(strlen(output_file_name)+10));
        if (output_file_name_gnu_plot == NULL)
        {
            fprintf(stderr,
                    "Error(main): the malloc for output_file_name_gnu_plot"
                    " failed.\n");
            exit(EXIT_FAILURE);
        }
        snprintf(output_file_name_gnu_plot,
                 strlen(output_file_name)+10,
                 "%s.gnuplot",
                 output_file_name);
    }
    if (arguments.grtcode_file == NULL)
    {
        fprintf(stderr,
                "Error(main): the GRTcode file was not passed in"
                "correctly.\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        grt_filename = arguments.grtcode_file;
    }

    /*Read in the optical depth from the grtcode output file.*/
    printf("\nReading in the GRT output file %s.\n",
           grt_filename);
    readOpticalDepthFromGrtcodeFile(grt_filename,
                                    &grt_optical_depth);
    printf("Read complete.\n");

    /*Read in the RFM reference data.*/
    printf("\nReading in the RFM reference data from file %s.\n",
           rfm_filename);
    readOpticalDepthFromRfmFile(rfm_filename,
                                &rfm_optical_depth,
                                48,
                                3001);
    printf("Read complete.\n");

    /*Open files where the output will be written.*/
    printf("\nOpening the output file.\n");
    output_file = fopen(output_file_name,"w");
    if (output_file == NULL)
    {
        fprintf(stderr,
                "The output file %s did not open properly.\n",
                output_file_name);
        exit(EXIT_FAILURE);
    }
    output_file_gnu_plot = fopen(output_file_name_gnu_plot,"w");
    if (output_file_gnu_plot == NULL)
    {
        fprintf(stderr,
                "The output file %s did not open properly.\n",
                output_file_name_gnu_plot);
        exit(EXIT_FAILURE);
    }

    /*Get the array dimension bounds for the comparison.*/
    printf("\nGetting the array dimensions of the RFM and GRT comparison.\n");
    if (rfm_optical_depth.ntime <= grt_optical_depth.ntime)
    {
        time_end = rfm_optical_depth.ntime;
    }
    else
    {
        time_end = grt_optical_depth.ntime;
    }
    if (rfm_optical_depth.nlat <= grt_optical_depth.nlat)
    {
        lat_end = rfm_optical_depth.nlat;
    }
    else
    {
        lat_end = grt_optical_depth.nlat;
    }
    if (rfm_optical_depth.nlon <= grt_optical_depth.nlon)
    {
        lon_end = rfm_optical_depth.nlon;
    }
    else
    {
        lon_end = grt_optical_depth.nlon;
    }
    if (rfm_optical_depth.npfull <= grt_optical_depth.npfull)
    {
        pfull_end = rfm_optical_depth.npfull;
    }
    else
    {
        pfull_end = grt_optical_depth.npfull;
    }
    if (rfm_optical_depth.nwavenumber <= grt_optical_depth.nwavenumber)
    {
        wavenumber_end = rfm_optical_depth.nwavenumber;
        wavenumber_start = 1;
        rfm_shift = 0;
        grt_shift = -1;
    }
    else
    {
        wavenumber_end = grt_optical_depth.nwavenumber;
        wavenumber_start = 0;
        rfm_shift = 1;
        grt_shift = 0;
    }

    /*Write out information about the comparison.*/
    printf("The results will be compared at:\n");
    printf("Times:           %d - %lu\n",0,time_end-1);
    printf("Latitudes:       %d - %lu\n",0,lat_end-1);
    printf("Longitudes:      %d - %lu\n",0,lon_end-1);
    printf("Pressure layers: %d - %lu\n",0,pfull_end-1);
    printf("Wavenumbers:     %d - %lu\n",1,wavenumber_end-1);
    fprintf(output_file,"The results will be compared at:\n");
    fprintf(output_file,"Times:           %d - %lu\n",0,time_end-1);
    fprintf(output_file,"Latitudes:       %d - %lu\n",0,lat_end-1);
    fprintf(output_file,"Longitudes:      %d - %lu\n",0,lon_end-1);
    fprintf(output_file,"Pressure layers: %d - %lu\n",0,pfull_end-1);
    fprintf(output_file,"Wavenumbers:     %d - %lu\n",1,wavenumber_end-1);
    fprintf(output_file,"time, lat, lon, pfull, wavenumber, tau_grt, tau_rfm, diff, %%error\n");
    fprintf(output_file,"\n");

    /*Loop through the optical depth arrays and calcluate the percent
      error as:
      abs((rfm_optical_depth - grt_optical_depth))/rfm_optical_depth.*/
    grt_index = 0;
    rfm_index = 0;
    max_error = -1.0;
    num_rfm_opt_depth_non_zero = 0;
    num_errors_over_one_percent = 0;
    num_errors_over_ten_percent = 0;
    num_errors_over_fifty_percent = 0;
    num_errors_over_one_hundred_percent = 0;
    for (time=0;time<time_end;time++)
    {
        for (lat=0;lat<lat_end;lat++)
        {
            for (lon=0;lon<lon_end;lon++)
            {
                for (pfull=0;pfull<pfull_end;pfull++)
                {
                    for (wavenumber=wavenumber_start;wavenumber<wavenumber_end;wavenumber++)
                    {
                        /*Calculate the absolute and percent differences in the optical depths.*/
                        rfm_index = time*((rfm_optical_depth.nlat)*(rfm_optical_depth.nlon)*(rfm_optical_depth.npfull)*(rfm_optical_depth.nwavenumber)) +
                                    lat*((rfm_optical_depth.nlon)*(rfm_optical_depth.npfull)*(rfm_optical_depth.nwavenumber)) +
                                    lon*((rfm_optical_depth.npfull)*(rfm_optical_depth.nwavenumber)) +
                                    pfull*(rfm_optical_depth.nwavenumber) +
                                    wavenumber + rfm_shift;
                        grt_index = time*((grt_optical_depth.nlat)*(grt_optical_depth.nlon)*(grt_optical_depth.npfull)*(grt_optical_depth.nwavenumber)) +
                                    lat*((grt_optical_depth.nlon)*(grt_optical_depth.npfull)*(grt_optical_depth.nwavenumber)) +
                                    lon*((grt_optical_depth.npfull)*(grt_optical_depth.nwavenumber)) +
                                    pfull*(grt_optical_depth.nwavenumber) +
                                    wavenumber + grt_shift;
                        diff = rfm_optical_depth.OpticalDepth[rfm_index] -
                               grt_optical_depth.OpticalDepth[grt_index];

                        if (diff < 0.0)
                        {
                            diff = diff*(-1.0);
                        }
                        if (rfm_optical_depth.OpticalDepth[rfm_index] != 0.0)
                        {
                            num_rfm_opt_depth_non_zero++;
                            percent_error = (100.)*(diff/(rfm_optical_depth.OpticalDepth[rfm_index]));
                        }
                        else
                        {
                            percent_error = 0.0;
                        }

                        /*Count the number of percent errors that are greater
                          than or equal to 1, 10, 50, and 100 percent.*/
                        if (percent_error >= 1.0)
                        {
                            num_errors_over_one_percent++;
                        }
                        if (percent_error >= 10.0)
                        {
                            num_errors_over_ten_percent++;
                        }
                        if (percent_error >= 50.0)
                        {
                            num_errors_over_fifty_percent++;
                        }
                        if (percent_error >= 100.0)
                        {
                            num_errors_over_one_hundred_percent++;
                        }

                        /*Keep track of the maximum error.*/
                        if (percent_error > 0 && percent_error > max_error)
                        {
                            max_error = percent_error;
                        }

                        /*Write out values/differences to a file.*/
                        fprintf(output_file,
                                "%lu %lu %lu %lu %lu %e %e %e %e\n",
                                time,
                                lat,
                                lon,
                                pfull,
                                wavenumber,
                                grt_optical_depth.OpticalDepth[grt_index],
                                rfm_optical_depth.OpticalDepth[rfm_index],
                                diff,
                                percent_error);
                        if (percent_error >= 0.0)
                        {
                            fprintf(output_file_gnu_plot,
                                    "%lu %lu %lu %lu %lu %e %e %e %e\n",
                                    time,
                                    lat,
                                    lon,
                                    pfull,
                                    wavenumber,
                                    grt_optical_depth.OpticalDepth[grt_index],
                                    rfm_optical_depth.OpticalDepth[rfm_index],
                                    diff,
                                    percent_error);
                        }
/*
                        fprintf(output_file,"Time:       %lu\n",time);
                        fprintf(output_file,"lat:        %lu\n",lat);
                        fprintf(output_file,"lon:        %lu\n",lon);
                        fprintf(output_file,"pfull:      %lu\n",pfull);
                        fprintf(output_file,"wavenumber: %lu\n",wavenumber);
                        fprintf(output_file,"grt_index:  %lu\n",grt_index);
                        fprintf(output_file,"rfm_index:  %lu\n",rfm_index);
                        fprintf(output_file,"tau(grt)    %e\n",grt_optical_depth.OpticalDepth[grt_index]);
                        fprintf(output_file,"tau(rfm)    %e\n",rfm_optical_depth.OpticalDepth[rfm_index]);
                        fprintf(output_file,"diff:       %e\n",diff);
                        fprintf(output_file,"error:      %e\n",percent_error);
                        fprintf(output_file,"\n");
*/
                    }
                }
            }
        }
    }
    printf("Compared %lu elements.  Max error = %f%%.\n",
           (time_end*lat_end*lon_end*pfull_end*wavenumber_end),
           max_error);
    printf("Number of RFM optical depth values > 0      = %lu (%f%%)\n",
           num_rfm_opt_depth_non_zero,
           100.0*((double)(num_rfm_opt_depth_non_zero))/
           ((double)(time_end*lat_end*lon_end*pfull_end*wavenumber_end)));
    printf("Number of optical depth differences >= 1%%   = %lu (%f%%)\n",
           num_errors_over_one_percent,
           100.0*((double)(num_errors_over_one_percent))/
           ((double)(num_rfm_opt_depth_non_zero)));
    printf("Number of optical depth differences >= 10%%  = %lu (%f%%)\n",
           num_errors_over_ten_percent,
           100.0*((double)(num_errors_over_ten_percent))/
           ((double)(num_rfm_opt_depth_non_zero)));
    printf("Number of optical depth differences >= 50%%  = %lu (%f%%)\n",
           num_errors_over_fifty_percent,
           100.0*((double)(num_errors_over_fifty_percent))/
           ((double)(num_rfm_opt_depth_non_zero)));
    printf("Number of optical depth differences >= 100%% = %lu (%f%%)\n",
           num_errors_over_one_hundred_percent,
           100.0*((double)(num_errors_over_one_hundred_percent))/
           ((double)(num_rfm_opt_depth_non_zero)));

    fprintf(output_file,"Compared %lu elements.  Max error = %e.\n",
            (time_end*lat_end*lon_end*pfull_end*wavenumber_end),
            max_error);
    fprintf(output_file,"Number of RFM optical depth values > 0      = %lu (%f%%)\n",
           num_rfm_opt_depth_non_zero,
           100.0*((double)(num_rfm_opt_depth_non_zero))/
           ((double)(time_end*lat_end*lon_end*pfull_end*wavenumber_end)));
    fprintf(output_file,"Number of optical depth differences >= 1%%   = %lu (%f%%)\n",
           num_errors_over_one_percent,
           100.0*((double)(num_errors_over_one_percent))/
           ((double)(num_rfm_opt_depth_non_zero)));
    fprintf(output_file,"Number of optical depth differences >= 10%%  = %lu (%f%%)\n",
           num_errors_over_ten_percent,
           100.0*((double)(num_errors_over_ten_percent))/
           ((double)(num_rfm_opt_depth_non_zero)));
    fprintf(output_file,"Number of optical depth differences >= 50%%  = %lu (%f%%)\n",
           num_errors_over_fifty_percent,
           100.0*((double)(num_errors_over_fifty_percent))/
           ((double)(num_rfm_opt_depth_non_zero)));
    fprintf(output_file,"Number of optical depth differences >= 100%% = %lu (%f%%)\n",
           num_errors_over_one_hundred_percent,
           100.0*((double)(num_errors_over_one_hundred_percent))/
           ((double)(num_rfm_opt_depth_non_zero)));

    printf("Verification complete.\n");

    /*Close the output files.*/
    printf("\nClosing the output file.\n");
    ioerr = fclose(output_file);
    if (ioerr != 0)
    {
        fprintf(stderr,
                "The output file %s did not close properly.\n",
                output_file_name);
        exit(EXIT_FAILURE);
    }
    output_file = NULL;
    ioerr = fclose(output_file_gnu_plot);
    if (ioerr != 0)
    {
        fprintf(stderr,
                "The output file %s did not close properly.\n",
                output_file_name_gnu_plot);
        exit(EXIT_FAILURE);
    }
    output_file_gnu_plot = NULL;

    /*Free and nullify the file names.*/
    free(output_file_name_gnu_plot);
    output_file_gnu_plot = NULL;
    output_file_name = NULL;
    grt_filename = NULL;
    rfm_filename = NULL;

    /*Free all malloc'd data.*/
    printf("\nFreeing remaining memory.\n");
    freeOpticalDepth(&grt_optical_depth);
    freeOpticalDepth(&rfm_optical_depth);
    printf("Free complete.\n\n");

    return EXIT_SUCCESS;
}
