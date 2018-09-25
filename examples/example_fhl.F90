program test
    use iso_c_binding
    use molecular_lines_fhl
    implicit none

!Utility macro for switching precision.  If you want to run in double
!precision, you must include the -DDOUBLE_PRECISION when building the
!library.
#ifdef SINGLE_PRECISION
#define FP c_float
#else
#define FP c_double
#endif

    character(kind=c_char,len=64) :: namelist_file !Path to the namelist
                                                   !file.
    integer(kind=c_int) :: num_levels !Number of atmospheric levels.
    integer(kind=c_int64_t) :: num_wpoints !Spectral grid size.
    integer(kind=c_int) :: num_molecules !Number of molecules the library
                                         !is using.  This must be equal to
                                         !the size of the hitran_files
                                         !array.
    integer(kind=c_int) :: num_layers !Number of atmospheric layers.  This
                                      !must be equal to the number of
                                      !atmospheric levels - 1.
    real(kind=FP),dimension(:),allocatable :: pressure !Atmospheric level
                                                       !pressures [atm].
    real(kind=FP),dimension(:),allocatable :: temperature !Atmospheric level
                                                          !temperatures [K].
    real(kind=FP),dimension(:,:),allocatable :: ppmv !Atmospheric level
                                                     !molecular abundances
                                                     ![ppmv].
    real(kind=FP),dimension(:,:),allocatable :: optical_depth !Atmospheric
                                                              !layer optical
                                                              !depths.
    integer(kind=c_int) :: num_columns !Number of atmospheric columns.
    integer(kind=c_int) :: i
    integer(kind=c_int) :: j
    integer(kind=c_int) :: k

    namelist_file = "namelist/example.nml"

    !Set verbosity.
    call grt_set_verbosity_f(2)

    !Initalize the library.
    call grt_context_init_fhl(namelist_filepath=trim(namelist_file))

    !Get the number of atmospheric levels.
    num_levels = grt_get_num_levels_fhl()

    !Get the size of the spectral grid.
    num_wpoints = grt_get_spectral_grid_size_fhl()

    !Allocate data arrays.
    num_molecules = grt_get_num_molecules_fhl()
    num_layers = num_levels - 1
    allocate(pressure(num_levels))
    allocate(temperature(num_levels))
    allocate(ppmv(num_levels,num_molecules))
    allocate(optical_depth(num_wpoints,num_layers))

    !Loop over columns.
    num_columns = 4
    do i = 1,num_columns

        !Make up some data.
        do j = 1,num_levels
            pressure(j) = 0.1 + 150.*real(j,kind=FP)
            temperature(j) = 230. + 2.3*real(j,kind=FP)
        enddo

        !Set molecular abundances.
        do k = 1,num_molecules
            do j = 1,num_levels
                ppmv(j,k) = 300. + 0.2*real(j,kind=FP) + real(k,kind=FP)
            enddo
        enddo

        !Calculate the optical depths.
        call grt_calculate_optical_depth_fhl(pressure, &
                                             temperature, &
                                             optical_depth, &
                                             xh2o=ppmv(:,1), &
                                             xco2=ppmv(:,2), &
                                             xo3=ppmv(:,3), &
                                             xn2o=ppmv(:,4), &
                                             xco=ppmv(:,5), &
                                             xch4=ppmv(:,6), &
                                             xo2=ppmv(:,7))
    enddo

    !Clean up.
    deallocate(pressure)
    deallocate(temperature)
    deallocate(ppmv)
    deallocate(optical_depth)

    !Free memory allocated by the context.
    call grt_context_free_fhl()

end program test
