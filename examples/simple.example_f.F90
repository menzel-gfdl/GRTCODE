program test
    use iso_c_binding
    use molecular_lines_f
    implicit none

!Utility macro for switching precision.  If you want to run in double
!precision, you must include the -DDOUBLE_PRECISION when building the
!library.
#ifdef DOUBLE_PRECISION
#define FP c_double
#else
#define FP c_float
#endif

    type(GrtContext_t) :: context !Library context.
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
    integer(kind=c_int) :: rc !Return code from library calls.
    integer(kind=c_int64_t) :: num_wpoints !Spectral grid size.
    integer(kind=c_int) :: num_layers !Number of atmospheric layers.  This
                                      !must be equal to the number of
                                      !atmospheric levels - 1.
    character(kind=c_char,len=64) :: hitran_path !Path to HITRAN database
                                                 !input file.
    integer(kind=c_int) :: h2o !Water vapor molecule id.  This is set by the
                               !library.
    integer(kind=c_int) :: o3 !Ozone molecule id.  This is set by the
                              !library.
    real(kind=FP),dimension(:),allocatable :: pressure !Atmospheric level
                                                       !pressures [atm].
    real(kind=FP),dimension(:),allocatable :: temperature !Atmospheric level
                                                          !temperatures [K].
    real(kind=FP),dimension(:),allocatable :: ppmv !Atmospheric level
                                                   !molecular abundances
                                                   ![ppmv].
    real(kind=FP),dimension(:,:),allocatable :: optical_depth !Atmospheric
                                                              !layer optical
                                                              !depths.
    integer(kind=c_int) :: num_columns !Number of atmospheric columns.
    integer(kind=c_int) :: i
    integer(kind=c_int) :: j

    num_levels = 25
    w0 = 1._c_double
    wn = 500._c_double
    wres = 0.1_c_double
    h2o_ctm_dir = "water_vapor_continuum"
    o3_ctm_dir = "ozone_continuum"

    !Initalize a library context context.
    rc = grt_context_init_f(context, &
                            num_levels, &
                            w0, &
                            wn, &
                            wres, &
                            h2o_ctm_dir=trim(h2o_ctm_dir), &
                            o3_ctm_dir=trim(o3_ctm_dir))
    call check_rc(rc)

    !Add water vapor to the context.
    hitran_path = ""
    hitran_path = "HITRAN_files/water_vapor.hitran12.par"
    rc = grt_add_molecule_f(context, &
                            trim(hitran_path), &
                            h2o, &
                            min_line_center_wavenumber=1._c_double, &
                            max_line_center_wavenumber=300._c_double)
    call check_rc(rc)

    !Add ozone to the context.
    hitran_path = ""
    hitran_path = "HITRAN_files/ozone.hitran12.par"
    rc = grt_add_molecule_f(context, &
                            trim(hitran_path), &
                            o3)
    call check_rc(rc)

    !Get the size of the spectral grid.
    rc = grt_get_spectral_grid_size_f(context, &
                                      num_wpoints)
    call check_rc(rc)
    num_layers = num_levels - 1

    !Allocate necessary arrays.
    allocate(pressure(num_levels))
    allocate(temperature(num_levels))
    allocate(ppmv(num_levels))
    allocate(optical_depth(num_wpoints,num_layers))

    !Loop over columns.
    num_columns = 3
    do i = 1,num_columns

        !Make up some data.
        do j = 1,num_levels
            pressure(j) = 0.1 + 150.*real(j,kind=FP)
            temperature(j) = 230. + 2.3*real(j,kind=FP)
            ppmv(j) = 300. + 0.2*real(j,kind=FP)
        enddo

        !Set water vapor abundance for the context.
        rc = grt_set_molecule_ppmv_f(context, &
                                     h2o, &
                                     ppmv)
        call check_rc(rc)

        !Make up some more data.
        do j = 1,num_levels
            ppmv(j) = 600. + 0.3*real(j,kind=FP)
        enddo

        !Set ozone abundance for the context.
        rc = grt_set_molecule_ppmv_f(context, &
                                     o3, &
                                     ppmv)
        call check_rc(rc)

        !Calculate the optical depths.
        rc = grt_calculate_optical_depth_f(context, &
                                           pressure, &
                                           temperature, &
                                           optical_depth)
        call check_rc(rc)
    enddo

    !Clean up.
    deallocate(pressure)
    deallocate(temperature)
    deallocate(ppmv)
    deallocate(optical_depth)

    !Free memory allocated by the context.
    rc = grt_context_free_f(context)
    call check_rc(rc)


    contains


    !Utility routine that traps errors returned by the library.
    subroutine check_rc(rc)
        use iso_fortran_env
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
