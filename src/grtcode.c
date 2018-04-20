#include <stdio.h>
#include <stdlib.h>
#include "arguments.h"
#include "constants.h"
#include "continuum.h"
#include "debug.h"
#include "input_fields.h"
#include "launch.h"
#include "model_fields.h"
#include "molecules.h"
#include "o3_continuum.h"
#include "parseHITRANfile.h"
#include "solar_flux.h"
#include "TIPS_2011.h"
#include "utils.h"
#include "write_output.h"
#include "omp.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif

int main(int argc,
         char* argv[])
{
    /*Set default command line argument values.*/
    struct arguments arguments;
    arguments.atmosInputFile = NULL;
    arguments.wingBreadth = 25;
    arguments.h2o_ctm = 0;
    arguments.o3_ctm = 0;
    arguments.device = DEFAULT_DEVICE;
    arguments.host = 0;
    arguments.outputFile = "defaultoutfile.nc";
    arguments.res = 1.0;
    arguments.t = MIN_TIME;
    arguments.T = MIN_TIME - 1;
    arguments.w = MIN_WVN;
    arguments.W = DEFAULT_WVN;
    arguments.x = MIN_LON;
    arguments.X = MIN_LON - 1;
    arguments.y = MIN_LAT;
    arguments.Y = MIN_LAT - 1;
    arguments.nHitFiles = 0;
    int i;
    for (i=0;i<NUM_MOL;++i)
    {
        arguments.molConc[i] = MISSING_CONC;
    }

    /*Parse the program's arguments using arg_parse.*/
    parse_options(argc,
                  argv,
                  &arguments);

    /*Make sure that only one target was specified (device or host).  If
      the device flag is given, an integer argument specifies the number of
      the device that will be used for the run (default is device 0).*/
    if (arguments.device != DEFAULT_DEVICE && arguments.host != 0)
    {
        fatal("more than one target specified (device=%d,host=%d).  Please"
                  " use either --host or --device or neither flag to just"
                  " default to device 0.",
              arguments.device,
              arguments.host);
    }
    int const launchType = arguments.host == 1 ? HOST_LAUNCH : DEVICE_LAUNCH;
    if (launchType == DEVICE_LAUNCH)
    {
        using_gpu();
#ifdef __NVCC__
        /*Set the GPU to be the host's current device.*/
        int const deviceNumber = arguments.device;
        HANDLE_ERROR(cudaSetDevice(deviceNumber));
#endif
    }

    /*Check input wavenumber bounds, input wavenumber resolution, and
      calculate the wavenumber grid size.*/
    check(input_bounds_check(arguments.w,
                             MIN_WVN,
                             &arguments.W,
                             MAX_WVN));
    if (arguments.res > RES_MAX || arguments.res < RES_MIN)
    {
        fatal("input wavenumber resolution (%e 1/cm) must be >= %e and"
                  " <= %e.",
              arguments.res,
              RES_MIN,
              RES_MAX);
    }
    unsigned int const nF = (arguments.W-arguments.w)/arguments.res + 1;
    log_mesg("Calculating optical depth spectra in range [%d,%d]"
                 " at resolution %e (1/cm).",
             arguments.w,
             arguments.W,
             arguments.res);

    /*Read in HITRAN line data.*/
    int const nMols = arguments.nHitFiles;
    line_params_t **hitLines = NULL;
    check(malloc_ptr((void **)(&hitLines),
                     sizeof(*hitLines)*nMols));
    line_flags_t flags = {((unsigned int) -1),1,0};
    int mol;
    for (mol=0;mol<nMols;++mol)
    {
        hitLines[mol] = NULL;
        check(parse_hitran_file(&(hitLines[mol]),
                                arguments.hitFiles[mol],
                                flags,
                                arguments.w,
                                arguments.W));
    }
    log_mesg("Submitted %d molecules:",
             nMols);
    for (mol=0;mol<nMols;++mol)
    {
        char mol_name[8];
        check(get_mol_name(hitLines[mol]->mol,
                           mol_name,
                           8));
        log_mesg("\t%s\t[%u lines]",
                 mol_name,
                 hitLines[mol]->nLines);
    }

    /*Make sure that the input molecular concentrations match the input
      hitran files.*/
    for (i=0;i<NUM_MOL;++i)
    {
        int found = 0;
        for (mol=0;mol<nMols;++mol)
        {
            if (i == (hitLines[mol])->mol)
            {
                found = 1;
                break;
            }
        }
        if (!found && arguments.molConc[i] != MISSING_CONC)
        {
            char mol_name[8];
            check(get_mol_name(i,
                               mol_name,
                               8));
            fatal("concentration for molecule %s was given on the command"
                      " line, but its corresponding HITRAN file was not"
                      " passed in.",
                  mol_name);
        }
    }

    /*Read in input data, process and store it in the form required by
      this model.*/
    log_mesg("Reading input data from file %s.",
             arguments.atmosInputFile);
    req_model_fields_t inputData;
    int mol_ids[nMols];
    for (mol=0;mol<nMols;++mol)
    {
        mol_ids[mol] = hitLines[mol]->mol;
    }
    check(get_input_data(&inputData,
                         arguments.atmosInputFile,
                         mol_ids,
                         nMols,
                         arguments.molConc));

    /*Check input time, longitude, and latitude bounds.*/
    check(input_bounds_check(arguments.t,
                             MIN_TIME,
                             &arguments.T,
                             inputData.ntime-1));
    check(input_bounds_check(arguments.x,
                             MIN_LON,
                             &arguments.X,
                             inputData.nlon-1));
    check(input_bounds_check(arguments.y,
                             MIN_LAT,
                             &arguments.Y,
                             inputData.nlat-1));

    /*Read in the solar flux values.*/
    SolarFlux_t solar_flux;
    check(get_solar_flux("INPUT/solar_flux.csv",
                         &solar_flux,
                         nF,
                         arguments.w,
                         arguments.res,
                         (launchType == DEVICE_LAUNCH)));

    /*Read in the water vapor continuum coefficients.*/
    ContinuumCoefs_t h2o_continuum;
    if (arguments.h2o_ctm)
    {
        /*Read in continuum coefficients and optionally put them on the
          device.*/
        check(get_h2o_continuum_coefs(&h2o_continuum,
                                      nF,
                                      arguments.w,
                                      arguments.res,
                                      (launchType == DEVICE_LAUNCH)));
    }

    /*Read in the ozone continuum coefficients.*/
    OzoneContinuumCoefs_t o3_continuum;
    if (arguments.o3_ctm)
    {
        /*Read in the ozone continuum coefficients.*/
        check(get_ozone_continuum_coefs("INPUT/ozone_continuum/"
                                            "ozone_continuum.csv",
                                        &o3_continuum,
                                        nF,
                                        arguments.w,
                                        arguments.res,
                                        (launchType == DEVICE_LAUNCH)));
    }

#ifdef __NVCC__
    /*Initialize TIPS.*/
    check(initTIPS_d());
#endif

#ifdef __NVCC__
    /*Declare CUDA stream parameters.*/
    int nstreams = -1;
    cudaStream_t *streams = NULL;
#endif

#ifdef _OPENMP
    /*Print out the number of OpenMP threads that will be used.*/
    log_mesg("Using %d OpenMP threads.",
             omp_get_max_threads());
#endif

    /*Allocate space for the data that will be output from the run.*/
    OutputFields_t out;
    check(alloc_output_fields(&out,
                              nF,
                              inputData.nlevel,
                              (launchType == DEVICE_LAUNCH)));

    /*Allocate/set pointers to buffers needed by the computation.*/
    WorkVars_t bufs;
    check(alloc_work_vars(&bufs,
                          inputData.nlevel,
                          MAX_NUM_LINES,
                          nF,
                          (launchType == DEVICE_LAUNCH)));

    /*Initialize output file.*/
    int outfile_ncid;
    check(init_output_file(arguments.outputFile,
                           &outfile_ncid,
                           arguments.X-arguments.x+1,
                           arguments.Y-arguments.y+1,
                           inputData.nlevel,
                           nF,
                           1));

    /*Loop over the atmospheric columns.*/
    int time;
    for (time=arguments.t;time<=arguments.T;++time)
    {
        int lon;
        for (lon=arguments.x;lon<=arguments.X;++lon)
        {
            int lat;
            for (lat=arguments.y;lat<=arguments.Y;++lat)
            {
                if (launchType == HOST_LAUNCH)
                {
                    check(launch_host(&bufs,
                                      &inputData,
                                      &solar_flux,
                                      time,
                                      lon,
                                      lat,
                                      nMols,
                                      hitLines,
                                      nF,
                                      (fp_t)arguments.w,
                                      arguments.res,
                                      arguments.wingBreadth,
                                      arguments.h2o_ctm,
                                      &h2o_continuum,
                                      arguments.o3_ctm,
                                      &o3_continuum,
                                      &out));
                }
                else if(launchType == DEVICE_LAUNCH)
                {
                    using_gpu();
#ifdef __NVCC__
                    check(launch_device(&bufs,
                                        &inputData,
                                        time,
                                        lon,
                                        lat,
                                        nMols,
                                        hitLines,
                                        nF,
                                        (fp_t)arguments.w,
                                        arguments.res,
                                        arguments.wingBreadth,
                                        arguments.h2o_ctm,
                                        &h2o_continuum,
                                        &out));
#endif
                }
                else
                {
                    fatal("invalid launch type (%d) requested, must be"
                              "either %d or %d.",
                          launchType,
                          DEVICE_LAUNCH,
                          HOST_LAUNCH);
                }

                /*Write out the column of output data.*/
                check(write_data_column(outfile_ncid,
                                        out.lw_flux_down,
                                        out.lw_flux_up,
                                        out.sw_flux_down,
                                        out.sw_flux_up,
                                        out.tau_gas,
                                        out.tau_scatter,
                                        time,
                                        lon,
                                        lat,
                                        inputData.nlevel,
                                        nF,
                                        1));
            }
        }
    }

    /*Close the output file.*/
    check(close_output_file(outfile_ncid));

    /*Free/nullify pointers to buffers needed by the computation.*/
    check(free_work_vars(&bufs,
                         (launchType == DEVICE_LAUNCH)));

    /*Free memory storing the data that was output from the run.*/
    check(free_output_fields(&out,
                             (launchType == DEVICE_LAUNCH)));

    /*Free memory storing the ozone continuum coefficients.*/
    if (arguments.o3_ctm)
    {
        check(free_ozone_continuum_coefs(&o3_continuum,
                                         (launchType == DEVICE_LAUNCH)));
    }

    /*Free memory storing the continuum coefficients.*/
    if (arguments.h2o_ctm)
    {
        check(free_continuum_coeffs(&h2o_continuum,
                                    (launchType == DEVICE_LAUNCH)));
    }

    /*Free memory storing the input solar flux values.*/
    check(free_solar_flux(&solar_flux,
                          (launchType == DEVICE_LAUNCH)));

    /*Free memory storing HITRAN line parameters.*/
    for (mol=0;mol<nMols;++mol)
    {
        check(free_line_params_host(&(hitLines[mol]),
                                    flags));
    }
    free(hitLines);

    /*Free memory storing input data.*/
    free_req_model_fields(&inputData);
    log_mesg("Run completed successfully, returning code %d.",
             SUCCESS);
    return SUCCESS;
}
