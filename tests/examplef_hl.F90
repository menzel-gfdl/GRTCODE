program test
    use iso_c_binding
    use molecular_lines_fhl
    implicit none

#ifdef DOUBLE_PRECISION
#define FP c_double
#else
#define FP c_float
#endif

    character(kind=c_char,len=64) :: namelist_file
    character(kind=c_char,len=64),dimension(2) :: hitran_files
    character(kind=c_char,len=64) :: h2o_ctm_dir
    character(kind=c_char,len=64) :: o3_ctm_dir
    integer(kind=c_int) :: num_levels
    integer(kind=c_int64_t) :: num_wpoints
    integer(kind=c_int) :: num_molecules
    integer(kind=c_int) :: num_layers
    real(kind=FP),dimension(:),allocatable :: pressure
    real(kind=FP),dimension(:),allocatable :: temperature
    real(kind=FP),dimension(:,:),allocatable :: ppmv
    real(kind=FP),dimension(:,:),allocatable :: optical_depth
    integer(kind=c_int) :: num_columns
    integer(kind=c_int) :: i
    integer(kind=c_int) :: j
    integer(kind=c_int) :: k

    !Initalize the library.
    namelist_file = "namelist/example.nml"
    hitran_files = (/"HITRAN_files/01_hit12.par","HITRAN_files/03_hit12.par"/)
    h2o_ctm_dir = "water_vapor_continuum"
    o3_ctm_dir = "ozone_continuum"
    call grt_context_init_fhl(hitran_files, &
                              namelist_filepath=trim(namelist_file), &
                              h2o_ctm_dir=trim(h2o_ctm_dir), &
                              o3_ctm_dir=trim(o3_ctm_dir))

    !Allocate data arrays.
    num_levels = grt_get_num_levels()
    num_wpoints = grt_get_spectral_grid_size()
    num_molecules = size(hitran_files)
    num_layers = num_levels - 1
    allocate(pressure(num_levels))
    allocate(temperature(num_levels))
    allocate(ppmv(num_levels,num_molecules))
    allocate(optical_depth(num_wpoints,num_layers))

    !Mimic looping over columns.
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
        call calculate_optical_depth_fhl(pressure, &
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
