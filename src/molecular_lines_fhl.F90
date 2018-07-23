!> @file
module molecular_lines_fhl
    use iso_c_binding
#ifdef _OPENMP
    use omp_lib
#endif
    use molecular_lines_f
    implicit none
    private


    !> \defgroup highlevelfortranapi High Level Fortran API
    !! \brief This module provides Fortran functions that calculate optical
    !!     depths from molecular lines in each layer of an atmospheric column.
    !!      To use this module, include the line:\n\n
    !!      ```use molecular_lines_fhl```\n\n
    !!      in your fortran application.
    !! \section Overview
    !!     The high level fortran API provides simplified interfaces, at the
    !!     cost of thread-safety and more fine-grained control.  The goal
    !!     of this API is to make the library more accessible to users
    !!     familiar with Fortran paradigms currently used in GCMs.
    !! \section Namelist
    !!     Runtime arguments can be supplied via a Fortran namelist titled
    !!     <c>molecular_lines_nml</c>.  A path to the file containing
    !!     the namelist may be passed as an optional argument to the
    !!     @ref grt_context_init_fhl routine.  An example namelist is
    !!     included in a directory named namelist in the base directory of
    !!     this repository.  If no path is provided, default
    !!     values are used.
    !!     \param num_levels <b> Integer(kind=c_int) </b> number of
    !!                       atmospheric levels.  Defaults to 60.
    !!     \param w0 <b> Real(kind=c_double) </b> lower bound [1/cm] of
    !!               spectral grid.  Defaults to 1.
    !!     \param wn <b> Real(kind=c_double) </b> upper bound [1/cm] of
    !!               spectral grid.  Defaults to 3250.
    !!     \param wres <b> Real(kind=c_double) </b> resolution [1/cm] of
    !!                 spectral grid.  Defaults to 0.1.
    !! \section Example
    !! \include examplef_hl.F90


    !Public routines
    public :: grt_context_init_fhl
    public :: grt_context_free_fhl
    public :: grt_calculate_optical_depth_fhl
    public :: grt_get_num_levels_fhl
    public :: grt_get_spectral_grid_size_fhl
    public :: grt_set_verbosity_f
    public :: grt_get_verbosity_f

#ifdef DOUBLE_PRECISION
#define FP c_double
#else
#define FP c_float
#endif


    !Private module variables
    type(GrtContext_t) :: context


    !Namelist variables.
    integer(kind=c_int) :: num_levels = 60 !Number of atmospheric levels.
    real(kind=c_double) :: w0 = 1._c_double !Lower bound of spectral grid.
    real(kind=c_double) :: wn = 3250._c_double !Upper bound of spectral grid.
    real(kind=c_double) :: wres = 0.1_c_double !Resolution of spectral grid.
    namelist /molecular_lines_nml/ num_levels, &
                                   w0, &
                                   wn, &
                                   wres


    contains


        !Utility routine that writes an error message to stderr.
        subroutine error(mesg)
            use iso_fortran_env

            !Inputs/outputs
            character(len=*),intent(in) :: mesg

            write(error_unit,*) "Error: "//trim(mesg)
        end subroutine error


        !Utility routine that catches when this library is called in an
        !OpenMP threaded region.
        subroutine omp_thread_trap()
#ifdef _OPENMP
            if (omp_get_level() .gt. 0) then
                if (omp_get_thread_num() .eq. 0) then
                    call error("This routine is not thread-safe.  You" &
                               //" cannot call it in an OpenMP region.")
                endif
!$omp barrier
                stop 1
            endif
#endif
        end subroutine


        !Utility routine that crashes the program if an error occurs in one
        !of the GRT routines.
        subroutine check_rc(rc)

            !Inputs/outputs
            integer(kind=c_int),intent(in) :: rc

            !Local variables
            character(len=1024) :: err_mesg

            if (rc .ne. 0) then
                call grt_errstr_f(rc, &
                                  err_mesg)
                call error(trim(err_mesg))
                stop 1
            endif
        end subroutine check_rc


        !> @ingroup highlevelfortranapi
        !! @brief Return the number of atmospheric levels.
        function grt_get_num_levels_fhl() result(n)

            !Inputs/outputs
            integer(kind=c_int) :: n

            n = num_levels
        end function grt_get_num_levels_fhl


        !> @ingroup highlevelfortranapi
        !! @brief Return the size of the spectral grid.
        function grt_get_spectral_grid_size_fhl() result(n)

            !Inputs/outputs
            integer(kind=c_int64_t) :: n

            !Local variables
            integer(kind=c_int) :: return_code

            return_code = grt_get_spectral_grid_size_f(context, &
                                                       n)
            call check_rc(return_code)
        end function grt_get_spectral_grid_size_fhl


        !> @ingroup highlevelfortranapi
        !! @brief Initialize the library with the molecules that correspond
        !!        to the input HITRAN file paths.
        subroutine grt_context_init_fhl(hitran_filepaths, &
                                        namelist_filepath, &
                                        wcutoff, &
                                        gpu_id, &
                                        num_threads, &
                                        h2o_ctm_dir, &
                                        o3_ctm_dir)

            !Inputs/outputs
            character(len=*),dimension(:),intent(in) :: hitran_filepaths !< Array of paths to HITRAN
                                                                         !! database files.
            character(len=*),intent(in),optional :: namelist_filepath !< Path to namelist file.
            real(kind=c_double),intent(in),optional :: wcutoff !< Cutoff [1/cm] from spectral
                                                               !! line center.  If NULL, this
                                                               !! defaults to 25 [1/cm].
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

            !Local variables
            integer(kind=c_int) :: return_code
            logical :: nml_exists
            logical :: in_use
            integer(kind=c_int) :: io_status
            integer(kind=c_int) :: mol_id
            integer(kind=c_int) :: i

            call omp_thread_trap()

            if (present(namelist_filepath)) then
                !Read the namelist.
                inquire(file=trim(namelist_filepath), &
                        exist=nml_exists)
                if (nml_exists) then
                    do i = 10,99
                        inquire(unit=i, &
                                opened=in_use)
                        if (.not. in_use) then
                            open(unit=i, &
                                 file=trim(namelist_filepath), &
                                 action="read")
                            read(i, &
                                 molecular_lines_nml, &
                                 iostat=io_status)
                            close(i)
                            if (io_status .gt. 0) then
                                call error("Error while reading " &
                                           //trim(namelist_filepath))
                                stop 1
                            endif
                            exit
                        endif
                    enddo
                    if (i .eq. 100) then
                        call error("Could not find available file unit for " &
                                   //trim(namelist_filepath))
                        stop 1
                    endif
                else
                    call error(trim(namelist_filepath)//" does not exist.")
                    stop 1
                endif
            endif

            !Initialize the library context.
            return_code = grt_context_init_f(context, &
                                             num_levels, &
                                             w0, &
                                             wn, &
                                             wres, &
                                             wcutoff, &
                                             gpu_id, &
                                             num_threads, &
                                             trim(h2o_ctm_dir), &
                                             trim(o3_ctm_dir))
            call check_rc(return_code)

            !Add the molecules associated with the input HITRAN files to
            !the library context.
            do i = 1,size(hitran_filepaths)
                return_code = grt_add_molecule_f(context, &
                                                 hitran_filepaths(i), &
                                                 mol_id)
                call check_rc(return_code)
            enddo
        end subroutine grt_context_init_fhl


        !> @ingroup highlevelfortranapi
        !! @brief Release memory allocated by library.
        subroutine grt_context_free_fhl()

            !Local variables
            integer(kind=c_int) :: return_code

            call omp_thread_trap()
            return_code = grt_context_free_f(context)
            call check_rc(return_code)
        end subroutine grt_context_free_fhl


        !> @ingroup highlevelfortranapi
        !! @brief Calculate the total optical depth of each atmospheric
        !!        layer at each spectral grid point.
        !! @return 0 if completed successfully, or else an error code.
        subroutine grt_calculate_optical_depth_fhl(pressure, &
                                                   temperature, &
                                                   ppmv, &
                                                   optical_depth)

            !Inputs/outputs
            real(kind=FP),dimension(:),intent(in) :: pressure !< Array of atmospheric pressures [atm].
                                                              !! The size of this array must be
                                                              !! equal to the number of atmospheric
                                                              !! levels.
            real(kind=FP),dimension(:),intent(in) :: temperature !< Array of atmospheric temperatures [K].
                                                                 !! The size of this array must be
                                                                 !! equal to the number of atmospheric
                                                                 !! levels.
            real(kind=FP),dimension(:,:),intent(in) :: ppmv !< Array of molecular abundances [ppmv].
                                                            !! This array must be of size (num_layers,molecules).
            real(kind=FP),dimension(:,:),intent(inout) :: optical_depth !< Array of atmospheric optical depths.
                                                                        !! The size of this array must be equal
                                                                        !! to the number of atmospheric layers
                                                                        !! times the number of spectral grid
                                                                        !! points.  Memory is layed out as
                                                                        !! (wavenumber,layer) (i.e., the
                                                                        !! fastest changing dimension is the
                                                                        !! one corresponding to the spectral
                                                                        !! grid.)
            !Local variables
            integer(kind=c_int) :: return_code
            integer(kind=c_int) :: i
            integer(kind=c_int) :: num_mols

            call omp_thread_trap()
            return_code = grt_get_num_molecules_f(context, &
                                                  num_mols)
            call check_rc(return_code)

            if (size(ppmv,1) .ne. num_levels .or. size(ppmv,2) .ne. &
                num_mols) then
                call error("input ppmv array must be of size" &
                           //"(num_levels,num_molecules).")
                stop 1
            endif
            do i = 1,num_mols
                return_code = grt_set_molecule_ppmv_f(context, &
                                                      i-1, &
                                                      ppmv(:,i))
                call check_rc(return_code)
            enddo
            if (size(pressure) .ne. num_levels .or. size(temperature) .ne. &
                num_levels) then
                call error("input pressure and temperature must" &
                           //" be of size num_levels.")
                stop 1
            endif
            return_code = grt_calculate_optical_depth_f(context, &
                                                        pressure, &
                                                        temperature, &
                                                        optical_depth)
            call check_rc(return_code)
        end subroutine grt_calculate_optical_depth_fhl


end module molecular_lines_fhl
