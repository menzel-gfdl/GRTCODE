program test
    use iso_c_binding
    use molecular_lines_fhl
    implicit none

!Utility macro for switching precision.  If you want to run in double
!precision, you must include the -DDOUBLE_PRECISION when building the
!library.
#ifdef DOUBLE_PRECISION
#define FP c_double
#else
#define FP c_float
#endif

    character(kind=c_char,len=64) :: namelist_file !Path to the namelist
                                                   !file.
    character(kind=c_char,len=64),dimension(2) :: hitran_files !Array of paths
                                                               !to HITRAN
                                                               !database
                                                               !input files.
    character(kind=c_char,len=64) :: h2o_ctm_dir !Path to directory
                                                 !containing the water
                                                 !vapor continuum input
                                                 !files.
    character(kind=c_char,len=64) :: o3_ctm_dir !Path to directory containing
                                                !ozone continuum input files.
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
    hitran_files = (/"HITRAN_files/water_vapor.hitran12.par", &
                     "HITRAN_files/ozone.hitran12.par      "/)
    h2o_ctm_dir = "water_vapor_continuum"
    o3_ctm_dir = "ozone_continuum"

    !Initalize the library.
    call grt_context_init_fhl(hitran_files, &
                              namelist_filepath=trim(namelist_file), &
                              h2o_ctm_dir=trim(h2o_ctm_dir), &
                              o3_ctm_dir=trim(o3_ctm_dir))

    !Get the number of atmospheric levels.
    num_levels = grt_get_num_levels_fhl()

    !Get the size of the spectral grid.
    num_wpoints = grt_get_spectral_grid_size_fhl()

    !Allocate data arrays.
    num_molecules = size(hitran_files)
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

        !Set molecular abundances.  Molecule index corresponds to its
        !place in the input hitran_files array.  For example, in this case
        !the hitran_files array was declare as:
        !hitran_files(1) = water vapor file.
        !hitran_files(2) = ozone input file.
        !so
        !ppmv(:,1) = water vapor abundances.
        !ppmv(:,2) = ozone abundances.
        do k = 1,num_molecules
            do j = 1,num_levels
                ppmv(j,k) = 300. + 0.2*real(j,kind=FP) + real(k,kind=FP)
            enddo
        enddo

        !Calculate the optical depths.
        call grt_calculate_optical_depth_fhl(pressure, &
                                             temperature, &
                                             ppmv, &
                                             optical_depth)
    enddo

    !Clean up.
    deallocate(pressure)
    deallocate(temperature)
    deallocate(ppmv)
    deallocate(optical_depth)

    !Free memory allocated by the context.
    call grt_context_free_fhl()

end program test
