!> @file
module molecular_lines
    use iso_c_binding
    implicit none
    private


    !> \defgroup lowlevelfortranapi Low Level Fortran API
    !! \section Overview
    !!     The low level fortran API simply provides bindings directly
    !!     to the functions described in @ref capi.  See the @ref capi
    !!     documentation for further details.
    !! \section Example
    !! \include examplef.F90
    public :: GrtContext_t
    public :: grt_context_init_f
    public :: grt_context_free_f
    public :: add_molecule_f
    public :: set_molecule_ppmv_f
    public :: calculate_optical_depth_f
    public :: grt_errstr_f


#ifdef DOUBLE_PRECISION
#define FP c_double
#else
#define FP c_float
#endif


    !> @ingroup lowlevelfortranapi
    !! @brief Library context.
    type GrtContext_t
        private
        type(c_ptr) :: p !< c pointer containing the address of a struct
                         !! where all data required by the context is
                         !! stored.
    end type GrtContext_t


    interface
        function grt_context_init(context, &
                                  num_levels, &
                                  w0, &
                                  wn, &
                                  wres, &
                                  num_wpoints, &
                                  wcutoff, &
                                  gpu_id, &
                                  num_threads, &
                                  h2o_ctm_dir, &
                                  o3_ctm_dir) &
            result(return_code) &
            bind(c)
            use iso_c_binding
            implicit none
            type(c_ptr),intent(inout) :: context
            integer(kind=c_int),value,intent(in) :: num_levels
            real(kind=c_double),value,intent(in) :: w0
            real(kind=c_double),value,intent(in) :: wn
            real(kind=c_double),value,intent(in) :: wres
            integer(kind=c_int64_t),intent(inout) :: num_wpoints
            real(kind=c_double),intent(in),optional :: wcutoff
            integer(kind=c_int),intent(in),optional :: gpu_id
            integer(kind=c_int),intent(in),optional :: num_threads
            character(kind=c_char,len=1),dimension(*),intent(in),optional :: h2o_ctm_dir
            character(kind=c_char,len=1),dimension(*),intent(in),optional :: o3_ctm_dir
            integer(kind=c_int) :: return_code
        end function grt_context_init
    end interface


    interface
        function grt_context_free(context) &
            result(return_code) &
            bind(c)
            use iso_c_binding
            implicit none
            type(c_ptr),intent(inout) :: context
            integer(kind=c_int) :: return_code
        end function grt_context_free
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
        function set_molecule_ppmv(context, &
                                   molecule_id, &
                                   ppmv) &
            result(return_code) &
            bind(c)
            use iso_c_binding
            implicit none
            type(c_ptr),value,intent(in) :: context
            integer(kind=c_int),value,intent(in) :: molecule_id
            real(kind=FP),dimension(*),intent(in) :: ppmv
            integer(kind=c_int) :: return_code
        end function set_molecule_ppmv
    end interface


    interface
        function calculate_optical_depth(context, &
                                         pressure, &
                                         temperature, &
                                         optical_depth) &
            result(return_code) &
            bind(c)
            use iso_c_binding
            implicit none
            type(c_ptr),value,intent(in) :: context
            real(kind=FP),dimension(*),intent(in) :: pressure
            real(kind=FP),dimension(*),intent(in) :: temperature
            real(kind=FP),dimension(*),intent(inout) :: optical_depth
            integer(kind=c_int) :: return_code
        end function calculate_optical_depth
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


        !> @ingroup lowlevelfortranapi
        !! @brief Initialize a context.
        !! @return 0 if completed successfully, or else an error code.
        function grt_context_init_f(context, &
                                    num_levels, &
                                    w0, &
                                    wn, &
                                    wres, &
                                    num_wpoints, &
                                    wcutoff, &
                                    gpu_id, &
                                    num_threads, &
                                    h2o_ctm_dir, &
                                    o3_ctm_dir) &
            result(return_code)
            type(GrtContext_t),intent(inout) :: context !< Library context.
            integer(kind=c_int),intent(in) :: num_levels !< Number of atmospheric levels.
            real(kind=c_double),intent(in) :: w0 !< Lowest wavenumber [1/cm] on spectral grid.
            real(kind=c_double),intent(in) :: wn !< Highest wavenumber [1/cm] on spectral grid.
            real(kind=c_double),intent(in) :: wres !< Spectral grid resolution [1/cm].
            integer(kind=c_int64_t),intent(inout) :: num_wpoints !< Number of points on
                                                                 !!spectral grid.
            real(kind=c_double),intent(in),optional :: wcutoff !< Cutoff [1/cm] from spectral
                                                               !! line center.  If NULL, this
                                                               !! defaults to 25 [1/cm].*/
            integer(kind=c_int),intent(in),optional :: gpu_id !< Id of the GPU that will be associated
                                                              !! with this context.  If NULL, then
                                                              !! use GPU 0 if at least one GPU
                                                              !! exists on the system, or else
                                                              !! set to -1 (corresponding to
                                                              !! a host only run.
            integer(kind=c_int),intent(in),optional :: num_threads !< If running on the host CPU,
                                                                   !! determines the maximum number
                                                                   !! of OpenMP threads that will
                                                                   !! be used.  Defaults to
                                                                   !! omp_get_max_threads (or one
                                                                   !! if not build with OpenMP).
            character(len=*),intent(in),optional :: h2o_ctm_dir !< Directory containing the
                                                                !! provided water vapor continuum
                                                                !! input files.  If NULL, then
                                                                !! the water vapor continuum
                                                                !! is not included in the optical
                                                                !! depth calculation.
            character(len=*),intent(in),optional :: o3_ctm_dir !< Directory containing the
                                                               !! provided ozone continuum
                                                               !! input files.  If NULL, then
                                                               !! the ozone continuum is not
                                                               !! included in the optical depth
                                                               !!  calculation.
            integer(kind=c_int) :: return_code

            return_code = grt_context_init(context%p, &
                                           num_levels, &
                                           w0, &
                                           wn, &
                                           wres, &
                                           num_wpoints, &
                                           wcutoff, &
                                           gpu_id, &
                                           num_threads, &
                                           trim(h2o_ctm_dir)//c_null_char, &
                                           trim(o3_ctm_dir)//c_null_char)
        end function grt_context_init_f


        !> @ingroup lowlevelfortranapi
        !! @brief Release memory allocated by the context.
        !! @return 0 if completed successfully, or else an error code.
        function grt_context_free_f(context) &
            result(return_code)
            type(GrtContext_t),intent(inout) :: context !< Library context.
            integer(kind=c_int) :: return_code

            return_code = grt_context_free(context%p)
        end function grt_context_free_f


        !> @ingroup lowlevelfortranapi
        !! @brief Add a molecule.  The optical depths of all added molecules
        !!        will be computed and summed to give the total optical optical
        !!        depth of each atmospheric layer at each spectral grid point.
        !! @return 0 if completed successfully, or else an error code.
        function add_molecule_f(context, &
                                hitran_filepath, &
                                molecule_id, &
                                min_line_center_wavenumber, &
                                max_line_center_wavenumber) &
            result(return_code)
            type(GrtContext_t),intent(in) :: context !< Library context.
            character(len=*),intent(in) :: hitran_filepath !< Path to HITRAN ascii file containing
                                                           !! molecular line parameters.
            integer(kind=c_int),intent(inout) :: molecule_id !< Id that is associated with the molecule.
            real(kind=c_double),intent(in),optional :: min_line_center_wavenumber !< Lower bound [1/cm] of spectral range.
                                                                                  !! Only lines with line center wavenumbers
                                                                                  !! greater than or eqaul to this will be
                                                                                  !! computed.  Defaults to 1 [1/cm].
            real(kind=c_double),intent(in),optional :: max_line_center_wavenumber !< Upper bound [1/cm] of spectral range.
                                                                                  !! Only lines with line center wavenumbers
                                                                                  !! less than or equal to this will be
                                                                                  !! computed.  Defaults to 50,000 [1/cm].
            integer(kind=c_int) :: return_code

            return_code = add_molecule(context%p, &
                                       trim(hitran_filepath)//c_null_char, &
                                       molecule_id, &
                                       min_line_center_wavenumber, &
                                       max_line_center_wavenumber)
        end function add_molecule_f


        !> @ingroup lowlevelfortranapi
        !! @brief Update the abundances [ppmv] for a molecule.
        !! @return 0 if completed successfully, or else an error code.
        function set_molecule_ppmv_f(context, &
                                     molecule_id, &
                                     ppmv) &
            result(return_code)
            type(GrtContext_t),intent(in) :: context !< Library context.
            integer(kind=c_int),intent(in) :: molecule_id !< Molecule id returned by
                                                          !! @ref add_molecule_f.
            real(kind=FP),dimension(*),intent(in) :: ppmv !< Array of molecular abundances [ppmv].
                                                          !! The size of this array must be
                                                          !! eqaul to the number of atmospheric
                                                          !! levels.
            integer(kind=c_int) :: return_code

            return_code = set_molecule_ppmv(context%p, &
                                            molecule_id, &
                                            ppmv)
        end function set_molecule_ppmv_f


        !> @ingroup lowlevelfortranapi
        !! @brief Calculate the total optical depth of each atmospheric
        !!        layer at each spectral grid point.
        !! @return 0 if completed successfully, or else an error code.
        function calculate_optical_depth_f(context, &
                                           pressure, &
                                           temperature, &
                                           optical_depth) &
            result(return_code)
            type(GrtContext_t),value,intent(in) :: context !< Library context.
            real(kind=FP),dimension(*),intent(in) :: pressure !< Array of atmospheric pressures [atm].
                                                              !! The size of this array must be
                                                              !! eqaul to the number of atmospheric
                                                              !! levels.
            real(kind=FP),dimension(*),intent(in) :: temperature !< Array of atmospheric temperatures [K].
                                                                 !! The size of this array must be
                                                                 !! eqaul to the number of atmospheric
                                                                 !! levels.
            real(kind=FP),dimension(*),intent(inout) :: optical_depth !< Array of atmospheric optical depths.
                                                                      !! The size of this array must be equal
                                                                      !! to the number of atmospheric layers
                                                                      !! times the number of spectral grid
                                                                      !! points.  Memory is layed out as
                                                                      !! (wavenumber,layer) (i.e., the
                                                                      !! fastest changing dimension is the
                                                                      !! one corresponding to the spectral
                                                                      !! grid.)
            integer(kind=c_int) :: return_code

            return_code = calculate_optical_depth(context%p, &
                                                  pressure, &
                                                  temperature, &
                                                  optical_depth)
        end function calculate_optical_depth_f


        !> @ingroup lowlevelfortranapi
        !! @brief Return a message describing the input return code.
        subroutine grt_errstr_f(code, &
                                buf)
            use iso_fortran_env

            integer(kind=c_int),intent(in) :: code !< Code returned from one of the
                                                   !! GRT functions.
            character(len=*),intent(inout) :: buf !< Buffer where message will be stored.

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
