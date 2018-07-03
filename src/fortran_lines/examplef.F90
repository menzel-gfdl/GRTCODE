program test
    use iso_c_binding
    use molecular_lines
    implicit none

#ifdef DOUBLE_PRECISION
#define FP c_double
#else
#define FP c_float
#endif

    type(c_ptr) :: context
    integer(kind=c_int) :: num_levels
    real(kind=c_double) :: w0
    real(kind=c_double) :: wn
    real(kind=c_double) :: wres
    integer(kind=c_int64_t) :: num_wpoints
    integer(kind=c_int) :: rc
    integer(kind=c_int) :: num_columns
    character(kind=c_char,len=64) :: h2o_hitran
    integer(kind=c_int) :: h2o
    real(kind=FP),dimension(:),allocatable :: pressure
    real(kind=FP),dimension(:),allocatable :: temperature
    real(kind=FP),dimension(:),allocatable :: ppmv
    real(kind=FP),dimension(:,:),allocatable :: optical_depth
    integer(kind=c_int) :: i
    integer(kind=c_int) :: j

    !Initialize parameters.
    num_levels = 25
    w0 = 2._c_double
    wn = 100._c_double
    wres = 0.1_c_double
    num_columns = 4
    h2o_hitran = "h2o_hit12.par"//c_null_char

    !Initalize library.
    rc = initialize_grt(context, &
                        num_levels, &
                        w0, &
                        wn, &
                        wres, &
                        num_wpoints, &
                        use_gpu=1, &
                        use_h2o_ctm=1)
    call check_rc(rc)

    !Add water vapor.
    rc = add_molecule(context, &
                      h2o_hitran, &
                      h2o, &
                      min_line_center_wavenumber=4._c_double, &
                      max_line_center_wavenumber=8._c_double)
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
        rc = set_molecule_ppmv(context, &
                               h2o, &
                               ppmv)
        call check_rc(rc)

        !Calculate the optical depths.
        rc = calculate_optical_depth(context, &
                                     pressure, &
                                     temperature, &
                                     optical_depth)
        call check_rc(rc)
!       write(*,*) optical_depth

    enddo

    !Clean up.
    deallocate(pressure)
    deallocate(temperature)
    deallocate(ppmv)
    deallocate(optical_depth)

    !Finalize library.
    rc = finalize_grt(context)
    call check_rc(rc)


    contains


    subroutine check_rc(rc)
        integer(kind=c_int),intent(in) :: rc
        if (rc .ne. 0) then
            stop 1
        endif
    end subroutine check_rc


end program test
