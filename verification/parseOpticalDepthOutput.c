#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "parseOpticalDepthOutput.h"

/*---------------------------------------------------------------------------*/
/*Handle errors by printing an error message and exiting with a
 *non-zero status. */
#define NCERR(e) {fprintf(stderr, "Error: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}

/*---------------------------------------------------------------------------*/
/*Read in the optical depth from the grtcode NetCDF output file.  The smallest
  wavenumber in the file is one, not zero like the RFM reference.*/
int readOpticalDepthFromGrtcodeFile(char fname[],
                                    OpticalDepth_t *GrtOutput)
{
    /*Local variables*/
    int retval;        /*NetCDF error code.*/
    int ncid;          /*NetCDF generated dataset id.*/
    int varid;         /*NetCDF variable id.*/
    int dimid;         /*NetCDF dimension id.*/
    int ndims_in;      /*Number of dimensions defined for the NetCDF dataset.*/
    int nvars_in;      /*Number of variables defined for the NetCDF dataset.*/
    int ngatts_in;     /*Number of global attributes defined for the NetCDF dataset.*/
    int unlimdimid_in; /*Id of the unlimited dimension for the NetCDF dataset.*/

    /*Open the grtcode NetCDF output file.*/
    if ((retval = nc_open(fname,
                          NC_NOWRITE,
                          &ncid)))
    {
        NCERR(retval);
    }

    /*Inquire about the file.  Get the number of dimensions, variables,
      and global attributes, and id of the unlimited dimension.*/
    if ((retval = nc_inq(ncid,
                         &ndims_in,
                         &nvars_in,
                         &ngatts_in,
                         &unlimdimid_in)))
    {
        NCERR(retval);
    }

    /*Get the time bounds from the NetCDF dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "time",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &GrtOutput->ntime)))
    {
        NCERR(retval);
    }

    /*Get the latitude coordinate bounds from the NetCDF dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "lat",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &GrtOutput->nlat)))
    {
        NCERR(retval);
    }

    /*Get the longitude coordinate bounds from the NetCDF dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "lon",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &GrtOutput->nlon)))
    {
        NCERR(retval);
    }

    /*Get the pressure layer midpoint (pfull) bounds from the NetCDF
      dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "pfull",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &GrtOutput->npfull)))
    {
        NCERR(retval);
    }

    /*Get the wavenumber bounds from the NetCDF dimension.*/
    if ((retval = nc_inq_dimid(ncid,
                               "wavenumber",
                               &dimid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_inq_dimlen(ncid,
                                dimid,
                                &GrtOutput->nwavenumber)))
    {
        NCERR(retval);
    }

    /*Null initialize the optical depth array.*/
    GrtOutput->OpticalDepth = NULL;

    /*Get the size of the optical depth array.*/
    const size_t opt_dep_len = (GrtOutput->ntime)*(GrtOutput->nlat)*
                               (GrtOutput->nlon)*(GrtOutput->npfull)*
                               (GrtOutput->nwavenumber);

    /*Make sure that the optical depth array is not too big.*/
    if (opt_dep_len > 5000000)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromGrtcodeFile): the optical"
                " depth array (%lu) is larger than 5 million elements.\n",
                opt_dep_len);
        exit(EXIT_FAILURE);
    }

    /*Malloc space for the optical depth array.*/
    if (GrtOutput->OpticalDepth != NULL)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromGrtcodeFile): the optical"
                " depth array is not null.\n");
        exit(EXIT_FAILURE);
    }
    GrtOutput->OpticalDepth = (REAL_t *)malloc(sizeof(REAL_t)*opt_dep_len);
    if (GrtOutput->OpticalDepth == NULL)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromGrtcodeFile): the optical"
                " depth array malloc failed.\n");
        exit(EXIT_FAILURE);
    }

    /*Get the optical depth variable id and read in the data.*/
    if ((retval = nc_inq_varid(ncid,
                               "OpticalDepth",
                               &varid)))
    {
        NCERR(retval);
    }
    if ((retval = nc_get_var_float(ncid,
                                   varid,
                                   GrtOutput->OpticalDepth)))
    {
        NCERR(retval);
    }

    /*Close the grtcode NetCDF output file.*/
    if ((retval = nc_close(ncid)))
    {
        NCERR(retval);
    }
    ncid = 0;

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Read in the optical depth from the RFM output file. These optical depth
  values include values calculated at wavenumber = 0.*/
int readOpticalDepthFromRfmFile(char fname[],
                                OpticalDepth_t *RfmOutput,
                                size_t num_layers,
                                size_t num_wavenumbers)
{
    /*Local variables*/
    struct stat status;        /*File status.*/
    size_t tmp_file_len;       /*Length of a file.*/
    FILE *tmp_file_ptr = NULL; /*File pointer.*/
    char *buffer = NULL;       /*Character buffer.*/
    int found_number;          /*Found number flag.*/
    size_t cols;               /*Number of columns on a line.*/
    int first_line;            /*First line flag.*/
    size_t tmp_num_lines;      /*Number of lines in a file.*/
    int ioerr;                 /*I/O error code.*/
    size_t i;                  /*Loop variable.*/

    /*Set the number of vertical layers and wavenumber grid points.*/
    RfmOutput->npfull = num_layers;
    RfmOutput->nwavenumber = num_wavenumbers;

    /*The RFM output file is only given at one time, latitude, and longitude
      point.*/
    RfmOutput->ntime = 1;
    RfmOutput->nlat = 1;
    RfmOutput->nlon = 1;

    /*Get the size of the file.*/
    stat(fname,
         &status);
    tmp_file_len = (size_t)(status.st_size) + 1;

    /*Make sure that file is a valid length.*/
    if (tmp_file_len == 0 || tmp_file_len > 10000000)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromRfmFile): the file %s"
                " does not have a vaild size (%lu bytes).\n",
                fname,
                tmp_file_len);
    }

    /*Open the file.*/
    tmp_file_ptr = fopen(fname,"r");
    if (tmp_file_ptr == NULL)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromRfmFile): the file %s"
                " did not open properly.\n",
                fname);
        exit(EXIT_FAILURE);
    }

    /*Read in the file and store it in a buffer. Count the number of lines
      in the file and the number of columns per line.  The file is setup so
      that the every line except the first line contains a separate wavenumber,
      followed by optical depth values for each columns starting at the first
      layer above the surface.*/
    buffer = (char *)malloc(sizeof(char)*tmp_file_len);
    tmp_num_lines = 0;
    found_number = 0;
    cols = 0;
    first_line = 1;
    for (i=0;i<tmp_file_len-1;i++)
    {
        buffer[i] = getc(tmp_file_ptr);

        /*Make sure the read in character is allowed in the file.*/
        if (buffer[i] != '0' &&
            buffer[i] != '1' &&
            buffer[i] != '2' &&
            buffer[i] != '3' &&
            buffer[i] != '4' &&
            buffer[i] != '5' &&
            buffer[i] != '6' &&
            buffer[i] != '7' &&
            buffer[i] != '8' &&
            buffer[i] != '9' &&
            buffer[i] != ' ' &&
            buffer[i] != '.' &&
            buffer[i] != 'E' &&
            buffer[i] != '+' &&
            buffer[i] != '-' &&
            buffer[i] != '\n')
        {
            fprintf(stderr,
                    "Error(readOpticalDepthFromRfmFile): An invalid"
                    " character (%c) was read in on line (%lu) from the"
                    " file %s.\n",
                    buffer[i],
                    tmp_num_lines+1,
                    fname);
            exit(EXIT_FAILURE);
        }

        /*Handle the end of a line.*/
        if (buffer[i] == '\n')
        {
            if (!first_line)
            {
                /*Make sure that the number of columns in the line is equal
                  to the inputted number of layers plus 1.*/
                if (cols != num_layers+1)
                {
                    fprintf(stderr,
                            "Error(readOpticalDepthFromRfmFile): the number"
                            " of columns (%lu) in line (%lu) of the file %s"
                            " does not match the number of inputted layers"
                            " (%lu) plus one\n.",
                            cols,
                            tmp_num_lines+1,
                            fname,
                            num_layers);
                    exit(EXIT_FAILURE);
                }
            }
            first_line = 0;
            tmp_num_lines++;
            cols = 0;
            found_number = 0;
        }
        else if (!first_line)
        {
            if (!found_number)
            {
                /*Look for either a space or a number.*/
                if (buffer[i] != ' ' &&
                    buffer[i] != '0' &&
                    buffer[i] != '1' &&
                    buffer[i] != '2' &&
                    buffer[i] != '3' &&
                    buffer[i] != '4' &&
                    buffer[i] != '5' &&
                    buffer[i] != '6' &&
                    buffer[i] != '7' &&
                    buffer[i] != '8' &&
                    buffer[i] != '9')
                {
                    fprintf(stderr,
                            "Error(readOpticalDepthFromRfmFile): Found (%c)"
                            " on line %lu when expecting either a number or a"
                            " space.\n",
                            buffer[i],
                            tmp_num_lines+1);
                    exit(EXIT_FAILURE);
                }

                /*Found a number.*/
                if (buffer[i] == '0' ||
                    buffer[i] == '1' ||
                    buffer[i] == '2' ||
                    buffer[i] == '3' ||
                    buffer[i] == '4' ||
                    buffer[i] == '5' ||
                    buffer[i] == '6' ||
                    buffer[i] == '7' ||
                    buffer[i] == '8' ||
                    buffer[i] == '9')
                {
                    found_number = 1;
                    cols++;
                }

            }
            else if (buffer[i] == ' ')
            {
                /*Look for another number.*/
                found_number = 0;
            }
        }
    }
    buffer[tmp_file_len-1] = '\0';

    /*Make sure that the number of lines - 1 equals in the inputted number
      of wavenumber grid points.*/
    if (num_wavenumbers != tmp_num_lines - 1)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromRfmFile): the number of lines"
                "minus one (%lu) in the file %s does not match in the"
                "inputted number of wavenumber grid points (%lu).\n",
                (size_t)(tmp_num_lines - 1),
                fname,
                num_wavenumbers);
        exit(EXIT_FAILURE);
    }

    /*Close the file.*/
    ioerr = fclose(tmp_file_ptr);
    if (ioerr != 0)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromRfmFile): the file %s"
                " did not close properly.\n",
                fname);
        exit(EXIT_FAILURE);
    }
    tmp_file_ptr = NULL;

    /*Null initialize the optical depth array.*/
    RfmOutput->OpticalDepth = NULL;

    /*Get the size of the optical depth array.*/
    const size_t opt_dep_len = (RfmOutput->ntime)*(RfmOutput->nlat)*
                               (RfmOutput->nlon)*(RfmOutput->npfull)*
                               (RfmOutput->nwavenumber);

    /*Make sure that the optical depth array is not too big.*/
    if (opt_dep_len > 5000000)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromRfmFile): the optical"
                " depth array (%lu) is larger than 5 million elements.\n",
                opt_dep_len);
        exit(EXIT_FAILURE);
    }

    /*Malloc space for the optical depth array.*/
    if (RfmOutput->OpticalDepth != NULL)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromRfmFile): the optical"
                " depth array is not null.\n");
        exit(EXIT_FAILURE);
    }
    RfmOutput->OpticalDepth = (REAL_t *)malloc(sizeof(REAL_t)*opt_dep_len);
    if (RfmOutput->OpticalDepth == NULL)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromRfmFile): the optical"
                " depth array malloc failed.\n");
        exit(EXIT_FAILURE);
    }

    /*Copy the data from the buffer into the optical depth array.*/
    first_line = 1;
    found_number = 0;
    size_t num_buff_count = 0;
    char number_buffer[64];
    int wave_num_flag = 1;
    float tmp_wave_number;
    unsigned int j;
    size_t layer_counter = 0;
    float tmp_opt_depth;
    int tmp_opt_depth_index;
    for (i=0;i<64;i++)
    {
        number_buffer[i] = '\0';
    }

    size_t tmp_buff_len = strlen(buffer)+1;
    for (i=0;i<tmp_buff_len;i++)
    {
        if (buffer[i] == '\n')
        {
            if (!first_line &&
                found_number &&
                !wave_num_flag)
            {
                tmp_opt_depth = atof(number_buffer);
                tmp_opt_depth_index = 0*(RfmOutput->nlat)*(RfmOutput->nlon)*(RfmOutput->npfull)*(RfmOutput->nwavenumber) +
                                      0*(RfmOutput->nlon)*(RfmOutput->npfull)*(RfmOutput->nwavenumber) +
                                      0*(RfmOutput->npfull)*(RfmOutput->nwavenumber) +
                                      ((num_layers-1)-layer_counter)*(RfmOutput->nwavenumber) +
                                      (int)tmp_wave_number;
                RfmOutput->OpticalDepth[tmp_opt_depth_index] = tmp_opt_depth;
            }

            first_line = 0;
            wave_num_flag = 1;
            found_number = 0;
            layer_counter = 0;

            /*Reset the number buffer and counter.*/
            for (j=0;j<num_buff_count;j++)
            {
                number_buffer[j] = '\0';
            }
            num_buff_count = 0;
        }
        else if (!first_line)
        {
            if (!found_number)
            {
                /*Found a number.*/
                if (buffer[i] == '0' ||
                    buffer[i] == '1' ||
                    buffer[i] == '2' ||
                    buffer[i] == '3' ||
                    buffer[i] == '4' ||
                    buffer[i] == '5' ||
                    buffer[i] == '6' ||
                    buffer[i] == '7' ||
                    buffer[i] == '8' ||
                    buffer[i] == '9')
                {
                    found_number = 1;
                }
            }
            if (found_number)
            {
                if (buffer[i] == ' ')
                {
                    found_number = 0;
                    number_buffer[num_buff_count] = '\0';
                    if (wave_num_flag)
                    {
                        tmp_wave_number = atof(number_buffer);
                        wave_num_flag = 0;
                    }
                    else
                    {
                        tmp_opt_depth = atof(number_buffer);
                        tmp_opt_depth_index = 0*(RfmOutput->nlat)*(RfmOutput->nlon)*(RfmOutput->npfull)*(RfmOutput->nwavenumber) +
                                              0*(RfmOutput->nlon)*(RfmOutput->npfull)*(RfmOutput->nwavenumber) +
                                              0*(RfmOutput->npfull)*(RfmOutput->nwavenumber) +
                                              ((num_layers-1)-layer_counter)*(RfmOutput->nwavenumber) +
                                              (int)tmp_wave_number;

                        RfmOutput->OpticalDepth[tmp_opt_depth_index] = tmp_opt_depth;
                        layer_counter++;
                    }

                    /*Reset the number buffer and counter.*/
                    for (j=0;j<num_buff_count;j++)
                    {
                        number_buffer[j] = '\0';
                    }
                    num_buff_count = 0;
                }
                else
                {
                    number_buffer[num_buff_count] = buffer[i];
                    num_buff_count++;
                }
            }
        }
    }

    /*Free the buffer.*/
    if (buffer == NULL)
    {
        fprintf(stderr,
                "Error(readOpticalDepthFromRfmFile): the buffer is null"
                " so it cannot be freed.\n");
        exit(EXIT_FAILURE);
    }
    free(buffer);
    buffer = NULL;

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*Free the optical depth array.*/
int freeOpticalDepth(OpticalDepth_t *OpticalDepth_struct)
{
    /*Make sure that the optical depth array is not null.*/
    if (OpticalDepth_struct->OpticalDepth == NULL)
    {
        fprintf(stderr,
                "Error(freeOpticalDepth): the input optical"
                " depth array is null.\n");
        exit(EXIT_FAILURE);
    }

    /*Free the optical depth array.*/
    free(OpticalDepth_struct->OpticalDepth);

    /*Nullify the optical depth pointer.*/
    OpticalDepth_struct->OpticalDepth = NULL;

    return EXIT_SUCCESS;
}

/*---------------------------------------------------------------------------*/

