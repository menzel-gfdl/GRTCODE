!> @file
module molecular_lines
    use iso_c_binding
    implicit none
    private


    !> \defgroup lowlevelfortranapi Low Level Fortran API
    !! \section Overview
    !!     The low level fortran API simply provides bindings directly
    !!     to the functions described in @ref capi.
    !! \section Example
    !! \include examplef.F90


    public :: initialize_grt_f
    public :: finalize_grt_f
    public :: add_molecule_f
    public :: set_molecule_ppmv_f
    public :: calculate_optical_depth_f
    public :: grt_errstr_f


#ifdef DOUBLE_PRECISION
#define FP c_double
#else
#define FP c_float
#endif


    interface
        !> \brief function 1
        function initialize_grt_f(context, &
                                  num_levels, &
                                  w0, &
                                  wn, &
                                  wres, &
                                  num_wpoints, &
                                  wcutoff, &
                                  use_gpu, &
                                  num_threads, &
                                  use_h2o_ctm, &
                                  use_o3_ctm) &
            result(return_code) &
            bind(c,name="initialize_grt")
            use iso_c_binding
            implicit none
            type(c_ptr),intent(inout) :: context
            integer(kind=c_int),value,intent(in) :: num_levels
            real(kind=c_double),value,intent(in) :: w0
            real(kind=c_double),value,intent(in) :: wn
            real(kind=c_double),value,intent(in) :: wres
            integer(kind=c_int64_t),intent(inout) :: num_wpoints
            real(kind=c_double),intent(in),optional :: wcutoff
            integer(kind=c_int),intent(in),optional :: use_gpu
            integer(kind=c_int),intent(in),optional :: num_threads
            integer(kind=c_int),intent(in),optional :: use_h2o_ctm
            integer(kind=c_int),intent(in),optional :: use_o3_ctm
            integer(kind=c_int) :: return_code
        end function initialize_grt_f
    end interface


    !> \brief function 1
    interface
        function finalize_grt_f(context) &
            result(return_code) &
            bind(c,name="finalize_grt")
            use iso_c_binding
            implicit none
            type(c_ptr),intent(inout) :: context
            integer(kind=c_int) :: return_code
        end function finalize_grt_f
    end interface


    interface
        function add_molecule(context, &
                              hitran_filepath, &
                              molecule_id, &
                              min_line_center_wavenumber, &
                              max_line_center_wavenumber) &
            result(return_code) &
            bind(c)
            use iso_c_binding
            implicit none
            type(c_ptr),value,intent(in) :: context
            character(kind=c_char,len=1),dimension(*),intent(in) :: hitran_filepath
            integer(kind=c_int),intent(inout) :: molecule_id
            real(kind=c_double),intent(in),optional :: min_line_center_wavenumber
            real(kind=c_double),intent(in),optional :: max_line_center_wavenumber
            integer(kind=c_int) :: return_code
        end function add_molecule
    end interface


    interface
        function set_molecule_ppmv_f(context, &
                                     molecule_id, &
                                     ppmv) &
            result(return_code) &
            bind(c,name="set_molecule_ppmv")
            use iso_c_binding
            implicit none
            type(c_ptr),value,intent(in) :: context
            integer(kind=c_int),value,intent(in) :: molecule_id
            real(kind=FP),dimension(*),intent(in) :: ppmv
            integer(kind=c_int) :: return_code
        end function set_molecule_ppmv_f
    end interface


    interface
        function calculate_optical_depth_f(context, &
                                           pressure, &
                                           temperature, &
                                           optical_depth) &
            result(return_code) &
            bind(c,name="calculate_optical_depth")
            use iso_c_binding
            implicit none
            type(c_ptr),value,intent(in) :: context
            real(kind=FP),dimension(*),intent(in) :: pressure
            real(kind=FP),dimension(*),intent(in) :: temperature
            real(kind=FP),dimension(*),intent(inout) :: optical_depth
            integer(kind=c_int) :: return_code
        end function calculate_optical_depth_f
    end interface


    interface
        function grt_errstr(code, &
                            buf, &
                            buf_size) &
            result(return_code) &
            bind(c)
            use iso_c_binding
            implicit none
            integer(kind=c_int),value,intent(in) :: code
            character(kind=c_char,len=1),dimension(*),intent(inout) :: buf
            integer(kind=c_int),value,intent(in) :: buf_size
            integer(kind=c_int) :: return_code
        end function grt_errstr
    end interface


    contains


        function add_molecule_f(context, &
                                hitran_filepath, &
                                molecule_id, &
                                min_line_center_wavenumber, &
                                max_line_center_wavenumber) &
            result(return_code)
            use iso_c_binding
            implicit none

            type(c_ptr),value,intent(in) :: context
            character(len=*),intent(in) :: hitran_filepath
            integer(kind=c_int),intent(inout) :: molecule_id
            real(kind=c_double),intent(in),optional :: min_line_center_wavenumber
            real(kind=c_double),intent(in),optional :: max_line_center_wavenumber
            integer(kind=c_int) :: return_code

            character(kind=c_char,len=1024) :: buf

            if (len_trim(hitran_filepath) .ge. len(buf)) then
                stop 1
            endif
            buf = ""
            buf = trim(hitran_filepath)//c_null_char
            return_code = add_molecule(context, &
                                       buf, &
                                       molecule_id, &
                                       min_line_center_wavenumber, &
                                       max_line_center_wavenumber)
        end function add_molecule_f


        subroutine grt_errstr_f(code, &
                                buf)
            use iso_fortran_env
            use iso_c_binding
            implicit none

            integer(kind=c_int),intent(in) :: code
            character(len=*),intent(inout) :: buf

            integer(kind=c_int) :: return_code

            buf = ""
            return_code = grt_errstr(code, &
                                     buf, &
                                     len(buf,kind=c_int))
            if (return_code .ne. 0) then
                write(error_unit,*) "GRT: error while getting error string."
                stop 1
            endif
        end subroutine grt_errstr_f


#ifdef FOO
    !> \defgroup highlevelfortranapi High Level Fortran API
    !! \section Overview
    !!     The high level fortran API provides a more simplified interface,
    !!     at the expense of more fine-grained control.


    !Namelist variables.
    integer(kind=c_int) :: num_levels = 60
    real(kind=c_double) :: w0 = 1._c_double
    real(kind=c_double) :: wn = 50000._c_double
    real(kind=c_double) :: wres = 0.1_c_double
    namelist /grt_nml/ num_levels, &
                       w0, &
                       wn, &
                       wres


    !Private variables.
    type(c_ptr) :: context


    subroutine check_rc(rc)
        integer(kind=c_int),intent(in) :: rc
        if (rc .ne. 0) then
            stop 1
        endif
    end subroutine check_rc


    subroutine grt_init(hitran_filepaths, &
                        molecule_ids, &
                        num_wpoints)
        character(len=*),dimension(:),intent(in) :: hitran_filepaths
        integer(kind=c_int),dimension(:),intent(in) :: molecule_ids
        integer(kind=c_int64_t),intent(inout) :: num_wpoints
        character(kind=c_char,len=1024) :: buf
        integer(kind=c_int) :: rc
        integer(kind=c_int) :: i
        if (size(hitran_filepaths) .ne. size(molecule_ids)) then
            stop 1
        endif
        do i = 10,99
            inquire()
            open(unit=i,file="grt_nml",action="read")
            read(i,grt_nml)
            close(i)
            exit
        enddo
        if (i .eq. 100) then
            stop 1
        endif


        rc = initialize_grt(context, &
                            num_levels, &
                            w0, &
                            wn, &
                            wres, &
                            num_wpoints)
        call check_rc(rc)
        do i = 1,size(hitran_filepaths)
            buf = ""
            buf = trim(hitran_filepaths(i))//c_null_char
            rc = add_molecule(context, &
                              buf, &
                              molecule_ids(i))
            call check_rc(rc)
        enddo

    end subroutine grt_init


    subroutine grt_end()
        integer(kind=c_int) :: rc
        rc = finalize_grt(context)
        call check_rc(rc)
    end subroutine grt_end

#endif


end module molecular_lines
