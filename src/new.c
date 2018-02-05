#include "debug.h"
#include "model_field.h"

int main(int argc,
         char* argv[])
{
    /*Set default command line argument values.*/
    struct arguments arguments;
    arguments.atmos_input_file = NULL;
    arguments.atmos_input_file_format = NULL;
    arguments.device = 0;
    arguments.host = 0;
    arguments.nhitfiles = 0;
    arguments.nmolConc = 0;
    arguments.nmolConcOver = 0;
    arguments.output_file = "defaultoutfile.nc";
    arguments.wingBreadth = 25;
    arguments.ctm = 0;
    arguments.w = 1;
    arguments.W = 3000;
    arguments.t = 0;
    arguments.T = -1;
    arguments.x = 0;
    arguments.X = -1;
    arguments.y = 0;
    arguments.Y = -1;
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

    /*Make sure that only one target was specified (device or host).  If
      mpi is being used, then specify the number of devices (if a GPU
      run is being performed) or number of CPU cores (if a host-only run is
      begin performed) will be used.*/
    if (arguments.device != 0 && arguments.host != 0)
    {
        fatal("more than one target specified (device=%d,host=%d).  Please"
                  " use either --host or --device or neither flag to just"
                  " default to device 0.",
              arguments.device,
              arguments.host);
    }
    int const launchType = arguments.host == 1 ? 0 : 1;
    if (launchType == 1)
    {
#ifdef __NVCC__
        /*Set the GPU to be the host's current device.*/
        const int device_number = arguments.device;
        HANDLE_ERROR(cudaSetDevice(device_number));
#endif
    }

    /*Initialize the fields required by the model.*/
    req_model_fields_t atmos_state;
    check(init_req_model_fields(&atmos_state));

#ifdef FOO
    /*Read in atmospheric data from the inputted netCDF file.*/
    radiationOutputFields_t atmosData;
    getAndSetAtmosFieldsFromFile(arguments.atmos_input_file,
                                 arguments.atmos_input_file_format,
                                 &atmosData);


























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
    printf("\nSubmitted %d molecules.\n",
           nMols);

    /*Initialize the output file.*/
    int ncid;
    int varid[11];
    char *OUTPUT_FNAME = NULL;
    unsigned int compute_lat_beg = 0;
    unsigned int compute_lat_end = atmosData.nlat;
    unsigned int compute_lon_beg = 0;
    unsigned int compute_lon_end = atmosData.nlon;
    unsigned int lat;

    if (arguments.mpi != 0)
    {
        OUTPUT_FNAME = (char *)malloc(strlen(arguments.output_file)+9);
        if (OUTPUT_FNAME == NULL)
        {
            fprintf(stderr,
                    "Error(main): malloc failed for %zu bytes of"
                        " OUTPUT_FNAME.\n",
                    strlen(arguments.output_file) + 9);
            exit(EXIT_FAILURE);
        }

        /*Split up the latitudes amongst ranks.*/
        lat = atmosData.nlat/world_size;
        if (lat*world_size != atmosData.nlat)
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
        if (compute_lat_end > atmosData.nlat)
        {
            compute_lat_end = atmosData.nlat;
        }
        compute_lon_beg = 0;
        compute_lon_end = atmosData.nlon;
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
                           varid,
                           OUTPUT_FNAME,
                           compute_lat_end - compute_lat_beg,
                           compute_lon_end - compute_lon_beg,
                           numLayers,
                           nF);

    /*Check time bounds.*/
    if (arguments.T < 0)
    {
        arguments.T = atmosData.ntime - 1;
    }
    else if ((size_t)arguments.T > atmosData.ntime-1)
    {
        fprintf(stderr,
                "Upper time bound %d excepts the maximum time level (%zu) in"
                    " the input file.\n",
                arguments.T,
                atmosData.ntime-1);
        exit(EXIT_FAILURE);
    }
    if (arguments.t < 0)
    {
        fprintf(stderr,
                "Lower time bound %d must be >= 0.\n",
                arguments.t);
        exit(EXIT_FAILURE);
    }
    else if (arguments.t > arguments.T)
    {
        fprintf(stderr,
                "Lower time bound %d cannot be > upper time bound %d.\n",
                arguments.t,
                arguments.T);
        exit(EXIT_FAILURE);
    }
    int time = arguments.T - arguments.t + 1;

    /*Check longitude bounds.*/
    if (arguments.X < 0)
    {
        arguments.X = atmosData.nlon - 1;
    }
    else if ((size_t)arguments.X > atmosData.nlon-1)
    {
        fprintf(stderr,
                "Upper longitude bound %d excepts the maximum longitude"
                    " index (%zu) in the input file.\n",
                arguments.X,
                atmosData.nlon-1);
        exit(EXIT_FAILURE);
    }
    if (arguments.x < 0)
    {
        fprintf(stderr,
                "Lower longitude bound %d must be >= 0.\n",
                arguments.x);
        exit(EXIT_FAILURE);
    }
    else if (arguments.x > arguments.X)
    {
        fprintf(stderr,
                "Lower longitude bound %d cannot be > upper longitude"
                    " bound %d.\n",
                arguments.x,
                arguments.X);
        exit(EXIT_FAILURE);
    }

    /*Check latitude bounds.*/
    if (arguments.Y < 0)
    {
        arguments.Y = atmosData.nlat - 1;
    }
    else if ((size_t)arguments.Y > atmosData.nlat-1)
    {
        fprintf(stderr,
                "Upper latitude bound %d excepts the maximum latitude"
                    " index (%zu) in the input file.\n",
                arguments.Y,
                atmosData.nlat-1);
        exit(EXIT_FAILURE);
    }
    if (arguments.y < 0)
    {
        fprintf(stderr,
                "Lower latitude bound %d must be >= 0.\n",
                arguments.y);
        exit(EXIT_FAILURE);
    }
    else if (arguments.y > arguments.Y)
    {
        fprintf(stderr,
                "Lower latitude bound %d cannot be > upper latitude"
                    " bound %d.\n",
                arguments.y,
                arguments.Y);
        exit(EXIT_FAILURE);
    }

    /*Setup HITRAN lines.*/
    RefLinePtrs_t HitLines[nMols];
    RefLine_flags_t flags= {((unsigned int) -1),1,0}; /*(host cuda malloc default,
                                                         host=True,
                                                         device=false)*/
    unsigned int mol;
    for (mol=0;mol<nMols;++mol)
    {
        HitLines[mol] = parseHITRANfile(hitFnameList[mol],
                                        flags,
                                        arguments.w,
                                        arguments.W);

        /*Check the molecular configurations.  For all molecules whose
          partial pressure is not taken from the input NetCDF file, calculate
          the partial pressure from the concentrations inputted on the
          command line.*/
/*
        checkMolConfig(&arguments,
                       HitLines[mol].mol,
                       &atmosData,
                       time);
*/
    }

    REAL_t *CS_h = NULL;
    REAL_t *CF_h = NULL;
    REAL_t *T0_h = NULL;
    REAL_t *T0F_h = NULL;
    REAL_t *CS_d = NULL;
    REAL_t *CF_d = NULL;
    REAL_t *T0_d = NULL;
    REAL_t *T0F_d = NULL;
    if (arguments.ctm == 1)
    {
        /*Read in continuum coefficients.*/
        CS_h = (REAL_t *)calloc(nF,sizeof(REAL_t));
        parseCKD("INPUT/continuum/296MTCKD25_S.ccf",
                 CS_h,
                 nF,
                 arguments.w,
                 arguments.res);

        CF_h = (REAL_t *)calloc(nF,sizeof(REAL_t));
        parseCKD("INPUT/continuum/296MTCKD25_F.ccf",
                 CF_h,
                 nF,
                 arguments.w,
                 arguments.res);

        T0_h = (REAL_t *)calloc(nF,sizeof(REAL_t));
        parseCKD("INPUT/continuum/CKDS.ppp",
                 T0_h,
                 nF,
                 arguments.w,
                 arguments.res);

        T0F_h = (REAL_t *)calloc(nF,sizeof(REAL_t));
        parseCKD("INPUT/continuum/CKDF.ppp",
                 T0F_h,
                 nF,
                 arguments.w,
                 arguments.res);

#ifdef __NVCC__
        if (launchType == 1)
        {
            /*Make device copies.*/
            HANDLE_ERROR(cudaMalloc(&CS_d,
                                    nF*sizeof(REAL_t)));
            HANDLE_ERROR(cudaMemcpy(CS_d,
                                    CS_h,
                                    nF*sizeof(REAL_t),
                                    cudaMemcpyHostToDevice));

            HANDLE_ERROR(cudaMalloc(&CF_d,
                                    nF*sizeof(REAL_t)));
            HANDLE_ERROR(cudaMemcpy(CF_d,
                                    CF_h,
                                    nF*sizeof(REAL_t),
                                    cudaMemcpyHostToDevice));

            HANDLE_ERROR(cudaMalloc(&T0_d,
                                    nF*sizeof(REAL_t)));
            HANDLE_ERROR(cudaMemcpy(T0_d,
                                    T0_h,
                                    nF*sizeof(REAL_t),
                                    cudaMemcpyHostToDevice));

            HANDLE_ERROR(cudaMalloc(&T0F_d,
                                    nF*sizeof(REAL_t)));
            HANDLE_ERROR(cudaMemcpy(T0F_d,
                                    T0F_h,
                                    nF*sizeof(REAL_t),
                                    cudaMemcpyHostToDevice));
        }
#endif
    }

    /*Write out dimension data.*/
/*
    writeDimensionData(ncid,
                       varid[0],
                       (size_t)(arguments.T - arguments.t + 1),
                       );

    writeDimensionData(ncid,
                       varid[1],
                       (size_t)(compute_lat_end - compute_lat_beg),
                       );

    writeDimensionData(ncid,
                       varid[2],
                       (size_t)(compute_lon_end - compute_lon_beg),
                       );

    writeDimensionData(ncid,
                       varid[3],
                       (size_t)(numLayers),
                       );
*/

    REAL_t *wvn = (REAL_t *)malloc(sizeof(REAL_t)*nF);
    int i;
    for (i=0;i<nF;i++)
    {
        wvn[i] = arguments.w + i*arguments.res;
    }
    writeDimensionData(ncid,
                       varid[5],
                       (size_t)(nF),
                       wvn);
    free(wvn);

    /*Declare stream parameters.*/
#ifdef __NVCC__
    int nstreams = -1;
    cudaStream_t *streams = NULL;
#endif

#ifdef _OPENMP
    /*Print out the number of OpenMP threads that will be used.*/
    fprintf(stdout,
            "\nUsing %d OpenMP threads.\n",
            omp_get_max_threads());
#endif

    /*Initialize TIPS.*/
#ifdef __NVCC__
    initTIPS_d();
#else
    initTIPS();
#endif

    /*Compute the spectra.*/
    unsigned int lon;
    REAL_t *out = NULL;
    REAL_t *fluxesDown = NULL;
    REAL_t *fluxesUp = NULL;
    REAL_t *fluxesDown_accumulated = NULL;
    REAL_t *fluxesUp_accumulated = NULL;

    for (time=arguments.t;time<=arguments.T;++time)
    {
        for (lat=arguments.y;lat<=arguments.Y;++lat)
        {
            for (lon=arguments.x;lon<=arguments.X;++lon)
            {
                if (launchType == 0)
                {
                    if (out == NULL)
                    {
                        out = (REAL_t*)calloc(nF*numLayers,
                                              sizeof(REAL_t));
                        fluxesDown = (REAL_t*)calloc(nF*(numLayers+1),
                                                     sizeof(REAL_t));
                        fluxesUp = (REAL_t*)calloc(nF*(numLayers+1),
                                                   sizeof(REAL_t));
                        fluxesDown_accumulated = (REAL_t*)calloc((numLayers+1),
                                                                 sizeof(REAL_t));
                        fluxesUp_accumulated = (REAL_t*)calloc((numLayers+1),
                                                               sizeof(REAL_t));
                    }

                    /*Calculate line spectra on the host.*/
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
                                lon,
                                arguments.ctm,
                                CS_h,
                                CF_h,
                                T0_h,
                                T0F_h,
                                fluxesDown,
                                fluxesUp);
                }
                else if (launchType == 1)
                {
#ifdef __NVCC__
                    if (out == NULL)
                    {
                        HANDLE_ERROR(cudaHostAlloc(&out,
                                                   nF*numLayers*sizeof(REAL_t),
                                                   cudaHostAllocDefault));
                        HANDLE_ERROR(cudaHostAlloc(&fluxesDown,
                                                   nF*(numLayers+1)*sizeof(REAL_t),
                                                   cudaHostAllocDefault));
                        HANDLE_ERROR(cudaHostAlloc(&fluxesUp,
                                                   nF*(numLayers+1)*sizeof(REAL_t),
                                                   cudaHostAllocDefault));
                        HANDLE_ERROR(cudaHostAlloc(&fluxesDown_accumulated,
                                                   (numLayers+1)*sizeof(REAL_t),
                                                   cudaHostAllocDefault));
                        HANDLE_ERROR(cudaHostAlloc(&fluxesUp_accumulated,
                                                   (numLayers+1)*sizeof(REAL_t),
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
                                  lon,
                                  arguments.ctm,
                                  CS_d,
                                  CF_d,
                                  T0_d,
                                  T0F_d,
                                  fluxesDown,
                                  fluxesUp);
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
                            "Error(main): unknown launch type (%d)"
                                " requested.\n",
                            launchType);
                    exit(EXIT_FAILURE);
                }

                /*Sum the fluxes. Should this be a kernel?*/
                sum_fluxes(nF,
                           numLayers+1,
                           fluxesDown,
                           fluxesDown_accumulated,
                           arguments.res);
                sum_fluxes(nF,
                           numLayers+1,
                           fluxesUp,
                           fluxesUp_accumulated,
                           arguments.res);

                /*Write out the output file.*/
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
                                                out,
                                                fluxesDown,
                                                fluxesUp,
                                                fluxesDown_accumulated,
                                                fluxesUp_accumulated);

                memset(out,
                       0,
                       nF*numLayers*sizeof(REAL_t));
            }
        }
    }

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
    if (launchType == 1)
    {
#ifdef __NVCC__
        cudaFreeHost(out);
        cudaFreeHost(fluxesDown);
        cudaFreeHost(fluxesUp);
        cudaFreeHost(fluxesDown_accumulated);
        cudaFreeHost(fluxesUp_accumulated);
#endif
    }
    else
    {
        free(out);
        free(fluxesDown);
        free(fluxesUp);
        free(fluxesDown_accumulated);
        free(fluxesUp_accumulated);
    }

    if (arguments.ctm == 1)
    {
        /*Free malloc'd arrays.*/
        free(CS_h);
        free(CF_h);
        free(T0_h);
        free(T0F_h);

#ifdef __NVCC__
        if (launchType == 1)
        {
            /*Free device copies.*/
            HANDLE_ERROR(cudaFree(CS_d));
            HANDLE_ERROR(cudaFree(CF_d));
            HANDLE_ERROR(cudaFree(T0_d));
            HANDLE_ERROR(cudaFree(T0F_d));
        }
#endif
    }


#endif

    return EXIT_SUCCESS;
}
