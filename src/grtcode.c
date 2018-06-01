#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "water_vapor_continuum.h"
#include "write_output.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#include "query_gpu.cuh"
#endif

#ifdef use_MPI
#include "mpi.h"
#endif


int main(int argc,
         char* argv[])
{
    /*Initialize MPI (if necessary).*/
    int rank = 0;
    int num_ranks = 1;
#ifdef use_MPI
    mpi_check(MPI_Init(NULL,NULL));
    mpi_check(MPI_Comm_size(MPI_COMM_WORLD,&num_ranks));
    mpi_check(MPI_Comm_rank(MPI_COMM_WORLD,&rank));
#endif

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
    arguments.workers = 1;
    arguments.write_spectra = 0;
    int i;
    for (i=0;i<NUM_MOL;++i)
    {
        arguments.molConc[i] = MISSING_CONC;
    }

    /*Parse the program's arguments using arg_parse.*/
    parse_options(argc,
                  argv,
                  &arguments);

    /*Set the number of workers.  More than one requires OpenMP.*/
    int num_workers = arguments.workers;
    if (num_workers < 1)
    {
        fatal("at least one worker (%d) is required.",
              num_workers);
    }
    else if (num_workers > 1)
    {
#ifndef _OPENMP
        fatal("OpenMP is required when the number of workers (%d) > 1.",
              num_workers);
#endif
    }

    /*Set configuration for how the program will run (host vs. device).*/
    int launch_type;
    if (arguments.host)
    {
        launch_type = HOST_LAUNCH;
    }
    else
    {
        using_gpu();
        launch_type = DEVICE_LAUNCH;
        int num_devices;
        check(get_num_gpus(&num_devices,1));
        if (num_devices < 1)
        {
            fatal("the number of devices found (%d) must be > 0.  Use"
                      " --host for to run on only CPUs.",
                  num_devices);
        }
        if (num_workers > num_devices)
        {
            log_mesg("the number of workers requested (%d) exceeds the"
                         " number of devices found (%d).  All devices"
                         " will be used.",
                     num_workers,
                     num_devices);
            num_workers = num_devices;
        }
        if (num_workers == 1)
        {
            if (arguments.device < 0 || arguments.device >= num_devices)
            {
                fatal("the specified device (%d) does not exist on this"
                          " system.  Try a number in the range [0,%d].",
                      arguments.device,
                      num_devices-1);
            }
#ifdef __NVCC__
            HANDLE_ERROR(cudaSetDevice(arguments.device));
#endif
            log_mesg("Using device %d.",
                     arguments.device);
        }
        else
        {
            log_mesg("Using devices 0 - %d.",
                     num_workers-1);
        }
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
    unsigned int const nws = (arguments.W-arguments.w)/arguments.res + 1;
    log_mesg("Calculating optical depth spectra in range [%d,%d]"
                 " at resolution %e (1/cm).",
             arguments.w,
             arguments.W,
             arguments.res);

    /*Read in HITRAN line data.*/
    int const num_mols = arguments.nHitFiles;
    line_params_t **hit_lines = NULL;
    check(malloc_ptr((void **)(&hit_lines),
                     sizeof(*hit_lines)*num_mols));
    line_flags_t flags = {((unsigned int) -1),1,0};
    int mol;
    for (mol=0;mol<num_mols;++mol)
    {
        hit_lines[mol] = NULL;
        check(parse_hitran_file(&(hit_lines[mol]),
                                arguments.hitFiles[mol],
                                flags,
                                arguments.w,
                                arguments.W));
    }
    log_mesg("Submitted %d molecules:",
             num_mols);
    for (mol=0;mol<num_mols;++mol)
    {
        char mol_name[8];
        check(get_mol_name(hit_lines[mol]->mol,
                           mol_name,
                           8));
        log_mesg("\t%s\t[%u lines]",
                 mol_name,
                 hit_lines[mol]->nLines);
    }

    /*Make sure that the input molecular concentrations match the input
      hitran files.*/
    for (i=0;i<NUM_MOL;++i)
    {
        int found = 0;
        for (mol=0;mol<num_mols;++mol)
        {
            if (i == (hit_lines[mol])->mol)
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
    req_model_fields_t input_data;
    int mol_ids[num_mols];
    for (mol=0;mol<num_mols;++mol)
    {
        mol_ids[mol] = hit_lines[mol]->mol;
    }
    check(get_input_data(&input_data,
                         arguments.atmosInputFile,
                         mol_ids,
                         num_mols,
                         arguments.molConc));

    /*Check input time, longitude, and latitude bounds.*/
    check(input_bounds_check(arguments.t,
                             MIN_TIME,
                             &arguments.T,
                             input_data.ntime-1));
    check(input_bounds_check(arguments.x,
                             MIN_LON,
                             &arguments.X,
                             input_data.nlon-1));
    check(input_bounds_check(arguments.y,
                             MIN_LAT,
                             &arguments.Y,
                             input_data.nlat-1));

    /*Determine the bounds for the rank.*/
    double num_times = arguments.T - arguments.t + 1;
    int time_chunk = (int)(ceil(num_times/((double)num_ranks)));
    int t_start = rank*time_chunk + arguments.t;
    int t_end;
    if (t_start > arguments.T)
    {
        /*This rank will not participate in the calculation.*/
        t_end = t_start - 1;
    }
    else
    {
        t_end = t_start + time_chunk - 1;
        if (t_end > arguments.T)
        {
            t_end = arguments.T;
        }
    }
    log_mesg("rank=%d, t=%d, T=%d, workers=%d\n",
             rank,
             t_start,
             t_end,
             num_workers);
    mpi_check(MPI_Barrier(MPI_COMM_WORLD));

    /*Read in the solar flux values.*/
    SolarFlux_t solar_flux_h;
    check(get_solar_flux(&solar_flux_h,
                         nws,
                         arguments.w,
                         arguments.res));
    SolarFlux_t *solar_flux = NULL;
    if (launch_type == HOST_LAUNCH)
    {
        solar_flux = &solar_flux_h;
    }
    else
    {
        check(malloc_ptr((void **)(&solar_flux),
                         sizeof(*solar_flux)*num_workers));
        for (i=0;i<num_workers;++i)
        {
#ifdef __NVCC__
            HANDLE_ERROR(cudaSetDevice(i));
#endif
            check(put_solar_flux_on_device(&solar_flux_h,
                                           &(solar_flux[i])));
        }
    }

    /*Read in the water vapor continuum coefficients.*/
    WaterVaporContinuumCoefs_t h2o_continuum_h;
    WaterVaporContinuumCoefs_t *h2o_continuum = NULL;
    if (arguments.h2o_ctm)
    {
        check(get_water_vapor_continuum_coefs(&h2o_continuum_h,
                                              nws,
                                              arguments.w,
                                              arguments.res));
        if (launch_type == HOST_LAUNCH)
        {
            h2o_continuum = &h2o_continuum_h;
        }
        else
        {
            check(malloc_ptr((void **)(&h2o_continuum),
                             sizeof(*h2o_continuum)*num_workers));
            for (i=0;i<num_workers;++i)
            {
#ifdef __NVCC__
                HANDLE_ERROR(cudaSetDevice(i));
#endif
                check(put_water_vapor_coefs_on_device(&h2o_continuum_h,
                                                      &(h2o_continuum[i])));
            }
        }
    }

    /*Read in the ozone continuum coefficients.*/
    OzoneContinuumCoefs_t o3_continuum_h;
    OzoneContinuumCoefs_t *o3_continuum = NULL;
    if (arguments.o3_ctm)
    {
        check(get_ozone_continuum_coefs(&o3_continuum_h,
                                        nws,
                                        arguments.w,
                                        arguments.res));
        if (launch_type == HOST_LAUNCH)
        {
            o3_continuum = &o3_continuum_h;
        }
        else
        {
            check(malloc_ptr((void **)(&o3_continuum),
                             sizeof(*o3_continuum)*num_workers));
            for (i=0;i<num_workers;++i)
            {
#ifdef __NVCC__
                HANDLE_ERROR(cudaSetDevice(i));
#endif
                check(put_ozone_coefs_on_device(&o3_continuum_h,
                                                &(o3_continuum[i])));
            }
        }
    }

    /*Allocate space for the data that will be output from the run.*/
    OutputFields_t *out = NULL;
    int num_buffers;
    if (launch_type == HOST_LAUNCH)
    {
        num_buffers = 1;
    }
    else
    {
        num_buffers = num_workers;
    }
    check(malloc_ptr((void **)(&out),
                     sizeof(*out)*num_buffers));
    for (i=0;i<num_buffers;++i)
    {
        check(alloc_output_fields(&(out[i]),
                                  nws,
                                  input_data.nlevel,
                                  (launch_type == DEVICE_LAUNCH)));
    }

    /*Allocate/set pointers to buffers needed by the computation.*/
    WorkVars_h_t *bufs_h = NULL;
    WorkVars_t *bufs = NULL;
    if (launch_type == HOST_LAUNCH)
    {
        check(malloc_ptr((void **)(&bufs_h),
                         sizeof(*bufs_h)));
        check(alloc_work_vars_h(bufs_h,
                                input_data.nlevel,
                                MAX_NUM_LINES,
                                nws));
    }
    else
    {
        check(malloc_ptr((void **)(&bufs),
                         sizeof(*bufs)*num_buffers));
        for (i=0;i<num_buffers;++i)
        {
#ifdef __NVCC__
            HANDLE_ERROR(cudaSetDevice(i));
            check(alloc_work_vars(&(bufs[i]),
                                  input_data.nlevel,
                                  MAX_NUM_LINES,
                                  nws));
#endif
        }
    }

    /*Initialize output file.*/
    int const max_file_name_length = 128;
    char outfile_name[2*max_file_name_length];
    if (strlen(arguments.outputFile) > max_file_name_length)
    {
        fatal("name of the output file (%s) must be <= %d characters long.",
              arguments.outputFile,
              max_file_name_length);
    }
    if (num_ranks > 1)
    {
        snprintf(outfile_name,
                 2*max_file_name_length,
                 "%s.%d",
                 arguments.outputFile,
                 rank);
    }
    else
    {
        snprintf(outfile_name,
                 2*max_file_name_length,
                 "%s",
                 arguments.outputFile);
    }
    int outfile_ncid;
    check(init_output_file(outfile_name,
                           &outfile_ncid,
                           arguments.X-arguments.x+1,
                           arguments.Y-arguments.y+1,
                           input_data.nlevel,
                           nws,
                           arguments.write_spectra));

    /*Device specific initialization.*/
    int nstreams = -1;
#ifdef __NVCC__
    cudaStream_t *streams = NULL;
#endif
    if (launch_type == DEVICE_LAUNCH)
    {
        /*Initialize TIPS.*/
        for (i=0;i<num_workers;++i)
        {
#ifdef __NVCC__
            HANDLE_ERROR(cudaSetDevice(i));
            check(initTIPS_d());
#endif
        }
    }

#ifdef _OPENMP
    /*Print out the number of OpenMP threads that will be used.*/
    log_mesg("Using %d OpenMP threads per rank.",
             num_workers);
    omp_set_num_threads(num_workers);
    omp_set_nested(1);
#endif

    /*Loop over the atmospheric columns.*/
    int time;
    int lon;
    int lat;
    int *err;
    check(malloc_ptr((void **)&err,
                     sizeof(*err)*num_workers));
    memset(err,
           0,
           sizeof(*err)*num_workers);

#pragma omp parallel for num_threads(num_buffers) \
                         collapse(3) \
                         default(shared) \
                         private(time,lon,lat)
    for (time=t_start;time<=t_end;++time)
    {
        for (lon=arguments.x;lon<=arguments.X;++lon)
        {
            for (lat=arguments.y;lat<=arguments.Y;++lat)
            {
#ifdef _OPENMP
                int index = omp_get_thread_num();
#else
                int index = 0;
#endif
                int ecode;

                if (launch_type == HOST_LAUNCH)
                {
                    ecode = launch_h(bufs_h,
                                     &input_data,
                                     solar_flux,
                                     time,
                                     lon,
                                     lat,
                                     num_mols,
                                     hit_lines,
                                     nws,
                                     (fp_t)arguments.w,
                                     arguments.res,
                                     arguments.wingBreadth,
                                     arguments.h2o_ctm,
                                     h2o_continuum,
                                     arguments.o3_ctm,
                                     o3_continuum,
                                     out);
                }
                else if (launch_type == DEVICE_LAUNCH)
                {
#ifdef __NVCC__
                    HANDLE_ERROR(cudaSetDevice(index));
                    ecode = launch(&(bufs[index]),
                                   &input_data,
                                   &(solar_flux[index]),
                                   time,
                                   lon,
                                   lat,
                                   num_mols,
                                   hit_lines,
                                   nws,
                                   (fp_t)arguments.w,
                                   arguments.res,
                                   arguments.wingBreadth,
                                   arguments.h2o_ctm,
                                   &(h2o_continuum[index]),
                                   arguments.o3_ctm,
                                   &(o3_continuum[index]),
                                   &(out[index]));
#endif
                }
                err[index] = err[index] | ecode;

                /*Write out the column of output data.*/
#pragma omp critical (output)
                write_data_column(outfile_ncid,
                                  out[index].lw_flux_down,
                                  out[index].lw_flux_up,
                                  out[index].sw_flux_down,
                                  out[index].sw_flux_up,
                                  out[index].tau_gas,
                                  out[index].tau_scatter,
                                  time-t_start,
                                  lon-arguments.x,
                                  lat-arguments.y,
                                  input_data.nlevel,
                                  nws,
                                  arguments.write_spectra);
            }
        }
    }

    /*Check for errors.*/
    for (i=0;i<num_workers;++i)
    {
        if (err[i] != 0)
        {
            fatal("errors (code=%d) detected during column calculations.",
                  err);
        }
    }
    free(err);

    /*Close the output file.*/
    check(close_output_file(outfile_ncid));

    /*Free buffers needed by the computation.*/
    if (launch_type == HOST_LAUNCH)
    {
        check(free_work_vars_h(bufs_h));
        free(bufs_h);
    }
    else
    {
        for (i=0;i<num_buffers;++i)
        {
#ifdef __NVCC__
            HANDLE_ERROR(cudaSetDevice(i));
            check(free_work_vars(&(bufs[i])));
#endif
        }
        free(bufs);
    }

    /*Free memory storing the data that was output from the run.*/
    for (i=0;i<num_buffers;++i)
    {
        check(free_output_fields(&(out[i]),
                                 (launch_type == DEVICE_LAUNCH)));
    }
    free(out);

    /*Free memory storing the ozone continuum coefficients.*/
    if (arguments.o3_ctm)
    {
        if (launch_type == DEVICE_LAUNCH)
        {
            for (i=0;i<num_workers;++i)
            {
#ifdef __NVCC__
                HANDLE_ERROR(cudaSetDevice(i));
#endif
                check(remove_ozone_coefs_from_device(&(o3_continuum[i])));
            }
            free(o3_continuum);
        }
        check(free_ozone_continuum_coefs(&o3_continuum_h));
    }

    /*Free memory storing the continuum coefficients.*/
    if (arguments.h2o_ctm)
    {
        if (launch_type == DEVICE_LAUNCH)
        {
            for (i=0;i<num_workers;++i)
            {
#ifdef __NVCC__
                HANDLE_ERROR(cudaSetDevice(i));
#endif
                check(remove_water_vapor_coefs_from_device(&(h2o_continuum[i])));
            }
            free(h2o_continuum);
        }
        check(free_water_vapor_continuum_coeffs(&h2o_continuum_h));
    }

    /*Free memory storing the input solar flux values.*/
    if (launch_type == DEVICE_LAUNCH)
    {
        for (i=0;i<num_workers;++i)
        {
#ifdef __NVCC__
            HANDLE_ERROR(cudaSetDevice(i));
#endif
            check(remove_solar_flux_from_device(&(solar_flux[i])));
        }
        free(solar_flux);
    }
    check(free_solar_flux(&solar_flux_h));

    /*Free memory storing input data.*/
    free_req_model_fields(&input_data);

    /*Free memory storing HITRAN line parameters.*/
    for (mol=0;mol<num_mols;++mol)
    {
        check(free_line_params_host(&(hit_lines[mol]),
                                    flags));
    }
    free(hit_lines);

/*
#ifdef use_MPI
*/
    mpi_check(MPI_Finalize());
/*
#endif
*/

    log_mesg("Run completed successfully, returning code %d.",
             SUCCESS);
    return SUCCESS;
}
