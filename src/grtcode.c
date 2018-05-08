#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "arguments.h"
#include "constants.h"
#include "debug.h"
#include "device_launch.h"
#include "host_launch.h"
#include "input_fields.h"
#include "model_fields.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "parseHITRANfile.h"
#include "solar_flux.h"
#include "TIPS_2011.h"
#include "utils.h"
#include "write_output.h"
#include "water_vapor_continuum.h"

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#include "query_gpu.cuh"
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
    int launchType;
    int num_devices;
    if (arguments.host != 0)
    {
        if (arguments.device != DEFAULT_DEVICE)
        {
            fatal("more than one target specified (device=%d,host=%d)."
                      "  Please use either --host or --device or neither"
                      " flag (defaults to device 0).",
                  arguments.device,
                  arguments.host);
        }
        launchType = HOST_LAUNCH;
        num_devices = 1;
    }
    else
    {
        using_gpu();
        launchType = DEVICE_LAUNCH;
#ifdef _OPENMP
        /*Use all available devices via OpenMP.*/
        check(get_num_gpus(&num_devices,1));
#else
        /*Use only the input device.*/
        num_devices = 1;
        HANDLE_ERROR(cudaSetDevice(arguments.device));
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
    SolarFlux_t solar_flux_h;
    check(get_solar_flux(&solar_flux_h,
                         nF,
                         arguments.w,
                         arguments.res));
    SolarFlux_t *solar_flux = NULL;
    check(malloc_ptr((void **)(&solar_flux),
                     sizeof(*solar_flux)*num_devices));
    for (i=0;i<num_devices;++i)
    {
        if (launchType == DEVICE_LAUNCH)
        {
            HANDLE_ERROR(cudaSetDevice(i));
            check(put_solar_flux_on_device(&solar_flux_h,
                                           &(solar_flux[i])));
        }
        else
        {
            memcpy(&(solar_flux[i]),
                   &solar_flux_h,
                   sizeof(solar_flux_h));
        }
    }

    WaterVaporContinuumCoefs_t h2o_continuum_h;
    WaterVaporContinuumCoefs_t *h2o_continuum = NULL;
    if (arguments.h2o_ctm)
    {
        /*Read in the water vapor continuum coefficients.*/
        check(get_water_vapor_continuum_coefs(&h2o_continuum_h,
                                              nF,
                                              arguments.w,
                                              arguments.res));
        check(malloc_ptr((void **)(&h2o_continuum),
                         sizeof(*h2o_continuum)*num_devices));
        for (i=0;i<num_devices;++i)
        {
            if (launchType == DEVICE_LAUNCH)
            {
                HANDLE_ERROR(cudaSetDevice(i));
                check(put_water_vapor_coefs_on_device(&h2o_continuum_h,
                                                      &(h2o_continuum[i])));
            }
            else
            {
                memcpy(&(h2o_continuum[i]),
                       &h2o_continuum_h,
                       sizeof(h2o_continuum_h));
            }
        }
    }

    OzoneContinuumCoefs_t o3_continuum_h;
    OzoneContinuumCoefs_t *o3_continuum = NULL;
    if (arguments.o3_ctm)
    {
        /*Read in the ozone continuum coefficients.*/
        check(get_ozone_continuum_coefs(&o3_continuum_h,
                                        nF,
                                        arguments.w,
                                        arguments.res));
        check(malloc_ptr((void **)(&o3_continuum),
                         sizeof(*o3_continuum)*num_devices));
        for (i=0;i<num_devices;++i)
        {
            if (launchType == DEVICE_LAUNCH)
            {
                HANDLE_ERROR(cudaSetDevice(i));
                check(put_ozone_coefs_on_device(&o3_continuum_h,
                                                &(o3_continuum[i])));
            }
            else
            {
                memcpy(&(o3_continuum[i]),
                       &o3_continuum_h,
                       sizeof(o3_continuum_h));
            }
        }
    }

    /*Initialize output file.*/
    int outfile_ncid;
    check(init_output_file(arguments.outputFile,
                           &outfile_ncid,
                           arguments.X-arguments.x+1,
                           arguments.Y-arguments.y+1,
                           inputData.nlevel,
                           nF,
                           1));

    /*Allocate space for the data that will be output from the run.*/
    OutputFields_t *out = NULL;
    check(malloc_ptr((void **)(&out),
                     sizeof(*out)*num_devices));
    for (i=0;i<num_devices;++i)
    {
        check(alloc_output_fields(&(out[i]),
                                  nF,
                                  inputData.nlevel,
                                  (launchType == DEVICE_LAUNCH)));
    }

    /*Allocate/set pointers to buffers needed by the computation.*/
    WorkVars_t *bufs = NULL;
    WorkVars_h_t *bufs_h = NULL;
    if (launchType == DEVICE_LAUNCH)
    {
        using_gpu();
        check(malloc_ptr((void **)(&bufs),
                         sizeof(*bufs)*num_devices));
        for (i=0;i<num_devices;++i)
        {
            HANDLE_ERROR(cudaSetDevice(i));
            check(alloc_work_vars(&(bufs[i]),
                                  inputData.nlevel,
                                  MAX_NUM_LINES,
                                  nF));
        }
    }
    else
    {
        check(malloc_ptr((void **)(&bufs_h),
                         sizeof(*bufs_h)*num_devices));
        for (i=0;i<num_devices;++i)
        {
            check(alloc_work_vars_h(&(bufs_h[i]),
                                    inputData.nlevel,
                                    MAX_NUM_LINES,
                                    nF));
        }
    }

#ifdef __NVCC__
    /*Initialize TIPS.*/
    check(initTIPS_d());

    /*Declare CUDA stream parameters.*/
    int nstreams = -1;
    cudaStream_t *streams = NULL;
#endif

#if defined(_OPENMP) && !defined(__NVCC__)
    /*Print out the number of OpenMP threads that will be used.*/
    log_mesg("Using %d OpenMP threads.",
             omp_get_max_threads());
#endif

    /*Loop over the atmospheric columns.*/
    int time;
    int lon;
    int lat;

#if defined(__NVCC__) && defined(_OPENMP)
/*
                         shared(num_devices,arguments,launchType,bufs,bufs_h, \
                                inputData,solar_flux,hitLines,h2o_continuum, \
                                o3_continuum,out,outfile_ncid,stderr) \
*/
#pragma omp parallel for num_threads(num_devices) \
                         collapse(3) \
                         default(shared) \
                         private(time,lon,lat)
#endif
    for (time=arguments.t;time<=arguments.T;++time)
    {
        for (lon=arguments.x;lon<=arguments.X;++lon)
        {
            for (lat=arguments.y;lat<=arguments.Y;++lat)
            {
                int index = num_devices - 1;
                if (launchType == HOST_LAUNCH)
                {
                    launch_h(&(bufs_h[index]),
                             &inputData,
                             &(solar_flux[index]),
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
                             &(h2o_continuum[index]),
                             arguments.o3_ctm,
                             &(o3_continuum[index]),
                             &(out[index]));
                }
                else if (launchType == DEVICE_LAUNCH)
                {
                    using_gpu();
#ifdef _OPENMP
                    index = omp_get_thread_num();
                    HANDLE_ERROR(cudaSetDevice(index));
                    log_mesg("thread %d setting current device to %d.",
                             omp_get_thread_num(),
                             index);
#endif
                    launch(&(bufs[index]),
                           &inputData,
                           &(solar_flux[index]),
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
                           &(h2o_continuum[index]),
                           arguments.o3_ctm,
                           &(o3_continuum[index]),
                           &(out[index]));
                }

                /*Write out the column of output data.*/
#pragma omp critical (output)
                write_data_column(outfile_ncid,
                                  out[index].lw_flux_down,
                                  out[index].lw_flux_up,
                                  out[index].sw_flux_down,
                                  out[index].sw_flux_up,
                                  out[index].tau_gas,
                                  out[index].tau_scatter,
                                  time-arguments.t,
                                  lon-arguments.x,
                                  lat-arguments.y,
                                  inputData.nlevel,
                                  nF,
                                  0);
            }
        }
    }

    /*Close the output file.*/
    check(close_output_file(outfile_ncid));

    /*Free/nullify pointers to buffers needed by the computation.*/
    if (launchType == DEVICE_LAUNCH)
    {
        using_gpu();
#ifdef __NVCC__
        for (i=0;i<num_devices;++i)
        {
            log_mesg("freeing work buffers on device %d.",
                     i);
            HANDLE_ERROR(cudaSetDevice(i));
            check(free_work_vars(&(bufs[i])));
        }
#endif
        free(bufs);
    }
    else
    {
        for (i=0;i<num_devices;++i)
        {
            check(free_work_vars_h(&(bufs_h[i])));
        }
        free(bufs_h);
    }

    /*Free memory storing the data that was output from the run.*/
    for (i=0;i<num_devices;++i)
    {
#ifdef __NVCC__
       log_mesg("freeing output buffers on device %d.",
                i);
        HANDLE_ERROR(cudaSetDevice(i));
#endif
        check(free_output_fields(&(out[i]),
                                 (launchType == DEVICE_LAUNCH)));
    }
    free(out);

    /*Free memory storing the ozone continuum coefficients.*/
    if (arguments.o3_ctm)
    {
        check(free_ozone_continuum_coefs(&o3_continuum_h));
    }

    /*Free memory storing the continuum coefficients.*/
    if (arguments.h2o_ctm)
    {
        check(free_water_vapor_continuum_coeffs(&h2o_continuum_h));
    }

    /*Free memory storing the input solar flux values.*/
    check(free_solar_flux(&solar_flux_h));

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
