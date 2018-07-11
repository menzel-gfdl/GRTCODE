program test
    use iso_c_binding
    use molecular_lines
    implicit none

#ifdef DOUBLE_PRECISION
#define FP c_double
#else
#define FP c_float
#endif

    type(GrtContext_t) :: context
    integer(kind=c_int) :: num_levels
    real(kind=c_double) :: w0
    real(kind=c_double) :: wn
    real(kind=c_double) :: wres
    integer(kind=c_int64_t) :: num_wpoints
    integer(kind=c_int) :: rc
    integer(kind=c_int) :: num_columns
    character(kind=c_char,len=64) :: h2o_hitran
    integer(kind=c_int) :: h2o
    character(kind=c_char,len=64) :: h2o_ctm_dir
    character(kind=c_char,len=64) :: o3_hitran
    integer(kind=c_int) :: o3
    character(kind=c_char,len=64) :: o3_ctm_dir
    real(kind=FP),dimension(:),allocatable :: pressure
    real(kind=FP),dimension(:),allocatable :: temperature
    real(kind=FP),dimension(:),allocatable :: ppmv
    real(kind=FP),dimension(:,:),allocatable :: optical_depth
    integer(kind=c_int) :: i
    integer(kind=c_int) :: j

    !Initialize parameters.
    num_levels = 25
    w0 = 1._c_double
    wn = 500._c_double
    wres = 0.1_c_double
    num_columns = 4
    h2o_hitran = "HITRAN_FILES/01_hit12.par"
    h2o_ctm_dir = "water_vapor_continuum"
    o3_hitran = "HITRAN_FILES/03_hit12.par"
    o3_ctm_dir = "ozone_continuum"

    !Initalize a context.
    rc = grt_context_init_f(context, &
                            num_levels, &
                            w0, &
                            wn, &
                            wres, &
                            num_wpoints, &
                            h2o_ctm_dir=trim(h2o_ctm_dir), &
                            o3_ctm_dir=trim(o3_ctm_dir))
    call check_rc(rc)

    !Add water vapor to the context.
    rc = add_molecule_f(context, &
                        trim(h2o_hitran), &
                        h2o, &
                        min_line_center_wavenumber=1._c_double, &
                        max_line_center_wavenumber=300._c_double)
    call check_rc(rc)

    !Add ozone to the context.
    rc = add_molecule_f(context, &
                        trim(o3_hitran), &
                        o3, &
                        min_line_center_wavenumber=200._c_double, &
                        max_line_center_wavenumber=500._c_double)
    call check_rc(rc)

    !Mimic getting input data.
    allocate(pressure(num_levels))
    allocate(temperature(num_levels))
    allocate(ppmv(num_levels))
    allocate(optical_depth(num_wpoints,(num_levels-1)))

    !Mimic looping over columns.
    do i = 1,num_columns

        !Make up some data.
        do j = 1,num_levels
            pressure(j) = 0.1 + 150.*real(j,kind=FP)
            temperature(j) = 230. + 2.3*real(j,kind=FP)
            ppmv(j) = 300. + 0.2*real(j,kind=FP)
        enddo

        !Set water vapor ppmv.
        rc = set_molecule_ppmv_f(context, &
                                 h2o, &
                                 ppmv)
        call check_rc(rc)

        !Set ozone ppmv.
        do j = 1,num_levels
            ppmv(j) = 600. + 0.3*real(j,kind=FP)
        enddo
        rc = set_molecule_ppmv_f(context, &
                                 o3, &
                                 ppmv)
        call check_rc(rc)

        !Calculate the optical depths.
        rc = calculate_optical_depth_f(context, &
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
