program test
    use,intrinsic :: iso_fortran_env
    use,intrinsic :: iso_c_binding
#ifdef _OPENMP
    use omp_lib
#else
#error "You must build this code with OpenMP."
#endif
    use molecular_lines_f
    implicit none

!Utility macro for switching precision.  If you want to run in double
!precision, you must include the -DDOUBLE_PRECISION when building the
!library.
#ifdef SINGLE_PRECISION
#define FP c_float
#else
#define FP c_double
#endif

    type(GrtContext_t),dimension(:),allocatable :: context !Library context.
    integer(kind=c_int) :: num_levels !Number of atmospheric levels.
    real(kind=c_double) :: w0 !Lower bound [1/cm] of spectral grid.
    real(kind=c_double) :: wn !Upper bound [1/cm] of spectral grid.
    real(kind=c_double) :: wres !Resolution [1/cm] of spectral grid.
    character(kind=c_char,len=64) :: h2o_ctm_dir !Path to directory
                                                 !containing the water
                                                 !vapor continuum input
                                                 !files.
    character(kind=c_char,len=64) :: o3_ctm_dir !Path to directory containing
                                                !ozone continuum input files.
    integer(kind=c_int) :: num_contexts !Number of library contexts.
    logical :: host_only !Flag telling if the run will occur on only the
                         !host CPU.
    integer(kind=c_int),parameter :: host_id = -1 !Special id that signifies
                                                  !a host-only run.
    integer :: argc !Number of command line arguments.
    character(len=128) :: argv !Command line argument.
    integer(kind=c_int) :: rc !Return code from library calls.
    integer(kind=c_int) :: gid !GPU id for a context.
    real(kind=c_double) :: min_line_center_wavenumber !Lower bound for line
                                                      !center [1/cm].
    real(kind=c_double) :: max_line_center_wavenumber !Upper bound for line
                                                      !center [1/cm].
    integer(kind=c_int64_t) :: num_wpoints !Spectral grid size.
    integer(kind=c_int) :: num_layers !Number of atmospheric layers.  This
                                      !must be equal to the number of
                                      !atmospheric levels - 1.
    integer(kind=c_int) :: num_columns !Number of atmospheric columns.
    real(kind=FP),dimension(:,:),allocatable :: pressure !Atmospheric level
                                                         !pressures [atm].
    real(kind=FP),dimension(:,:),allocatable :: temperature !Atmospheric level
                                                            !temperatures [K].
    real(kind=FP),dimension(:,:),allocatable :: ppmv !Atmospheric level
                                                     !molecular abundances
                                                     ![ppmv].
    real(kind=FP),dimension(:,:,:),allocatable :: optical_depth !Atmospheric
                                                                !layer optical
                                                                !depths.
    integer(kind=c_int) :: i
    integer(kind=c_int) :: j

    !Command line argument controls how many GPUs will be used.
    num_contexts = 1
    host_only = .false.
    argc = command_argument_count()
    if (argc .eq. 1) then
        call get_command_argument(1,argv)
        if (trim(argv) .eq. "--host") then
            host_only = .true.
        else
            read(argv,*) num_contexts
        endif
    elseif (argc .gt. 1) then
        call get_command_argument(0,argv)
        write(error_unit,*) "Usage: "//trim(argv)//" [--host|num_gpus]"
        stop 1
    endif
    allocate(context(num_contexts))

    !Increase the verbosity of the library output.
    call grt_set_verbosity_f(3)

    num_levels = 25
    w0 = 1._c_double
    wn = 500._c_double
    wres = 0.1_c_double
    h2o_ctm_dir = "water_vapor_continuum"
    o3_ctm_dir = "ozone_continuum"

    do i = 1,num_contexts

        if (argc .eq. 0) then
            !Initalize library context pointers.
            rc = grt_context_init_f(context(i), &
                                    num_levels, &
                                    w0, &
                                    wn, &
                                    wres, &
                                    "HITRAN_files/hitran2012.par", &
                                    h2o_ctm_dir=h2o_ctm_dir, &
                                    o3_ctm_dir=o3_ctm_dir)
        else
            !Determine the GPU id for the context.
            if (host_only) then
                gid = host_id
            else
                gid = i - 1
            endif

            !Initalize library context pointers.
            rc = grt_context_init_f(context(i), &
                                    num_levels, &
                                    w0, &
                                    wn, &
                                    wres, &
                                    "HITRAN_files/hitran2012.par", &
                                    h2o_ctm_dir=h2o_ctm_dir, &
                                    o3_ctm_dir=o3_ctm_dir, &
                                    gpu_id=gid)
        endif
        call check_rc(rc)

        !Only water vapor lines with line centers in the range 1 - 1000 [1/cm]
        !will be included in the calculation.
        min_line_center_wavenumber = 1._c_double
        max_line_center_wavenumber = 1000._c_double

        !Add water vapor to the library context.
        rc = grt_add_molecule_f(context(i), &
                                H2O, &
                                min_line_center_wavenumber, &
                                max_line_center_wavenumber)
        call check_rc(rc)

        !Add ozone to the library context.
        rc = grt_add_molecule_f(context(i), &
                                O3)
        call check_rc(rc)
    enddo

    rc = grt_get_spectral_grid_size_f(context(1), &
                                      num_wpoints)
    call check_rc(rc)
    num_layers = num_levels - 1

    !Allocate necessary arrays.
    num_columns = 2*num_contexts
    allocate(pressure(num_levels,num_columns))
    allocate(temperature(num_levels,num_columns))
    allocate(ppmv(num_levels,num_columns))
    allocate(optical_depth(num_wpoints,num_layers,num_columns))

    !Loop over some columns.
!$omp parallel do num_threads(num_contexts) &
!$omp&            default(none) &
!$omp&            shared(context,num_levels,num_columns,num_wpoints, &
!$omp&                   pressure,temperature,ppmv,optical_depth) &
!$omp&            private(i,j,gid,rc)
    do i = 1,num_columns
        gid = omp_get_thread_num() + 1

        !Make up some data for the column.
        do j = 1,num_levels
            pressure(j,i) = 0.1 + 150.*j
            temperature(j,i) = 230. + 2.3*j
            ppmv(j,i) = 300. + 0.2*j
        enddo

        !Set the water vapor abundance for the library context.
        rc = grt_set_molecule_ppmv_f(context(gid), &
                                     H2O, &
                                     ppmv(:,i))
        call check_rc(rc)

        !Make up some more data for the column.
        do j = 1,num_levels
            ppmv(j,i) = 325. - 3.3*j
        enddo

        !Set the water vapor abundance for the library context.
        rc = grt_set_molecule_ppmv_f(context(gid), &
                                     O3, &
                                     ppmv(:,i))
        call check_rc(rc)

        !Calculate the optical depths.
        rc = grt_calculate_optical_depth_f(context(gid), &
                                           pressure(:,i), &
                                           temperature(:,i), &
                                           optical_depth(:,:,i))
        call check_rc(rc)
    enddo

    !Clean up.
    deallocate(pressure)
    deallocate(temperature)
    deallocate(ppmv)
    deallocate(optical_depth)

    !Free memory allocated by the library context.
    do i = 1,num_contexts
        rc = grt_context_free_f(context(i))
        call check_rc(rc)
    enddo
    deallocate(context)


    contains


    !Utility routine that traps errors returned by the library.
    subroutine check_rc(rc)
        integer(kind=c_int),intent(in) :: rc
        character(len=256) :: mesg
        if (rc .ne. 0) then
            call grt_errstr_f(rc, &
                              mesg)
            write(error_unit,*) trim(mesg)
            stop 1
        endif
    end subroutine check_rc


end program test
