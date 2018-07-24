!> @file
module molecular_lines_fhl
    use iso_c_binding
#ifdef _OPENMP
    use omp_lib
#endif
    use molecular_lines_f
    implicit none
    private


    !> \defgroup highlevelfortranapi High-level Fortran API
    !!@section Overview
    !!    Given an atmospheric column made up of at least one layer,
    !!    this code calculates the total optical depth of each layer
    !!    at each point on an input spectral grid.  Unlike the
    !!    @ref lowlevelfortranapi, this API provides simplified interfaces,
    !!    at the cost of thread-safety and more fine-grained control with the
    !!    goal of making this library more accessible to users
    !!    familiar with Fortran paradigms currently used in GCMs.
    !!    To use this API,
    !!
    !!        use molecular_lines_fhl
    !!
    !!    and follow these steps:
    !!        -# Initialize the library by calling the
    !!           @ref grt_context_init_fhl function.  Here you must provide
    !!           an array of paths HITRAN database files, where each path
    !!           in the array corresponds to a molecule that will be used
    !!           in the optical depth calculation.
    !!           @attention You must call this function before calling any
    !!               other function included in this library.  Failure to do
    !!               so will result in undefined behavior.
    !!
    !!           In addition, the following optional parameters may be set (not
    !!           passing these arguments in implies using the default values):
    !!           - A path to a file containg a Fortran namelist.  As described
    !!             below, the number of levels per atmospheric column and
    !!             parameters (lower bound, upper bound, and resolution) for
    !!             the spectral grid on which the optical depths will be
    !!             calculated are namelist controlled.  If this argument is
    !!             not passed in, the default values described in the section
    !!             below will be used.
    !!           - A cutoff value for the molecular lines.  This
    !!             value denotes how far (in terms of wavenumber) from the line
    !!             center each molecular line is calculated out to.  By default
    !!             this value is set to 25 [1/cm], so if for example a line
    !!             center lies at the wavenumber 150 [1/cm], then it contributes
    !!             to the optical depth at all spectral grid points in the range
    !!             125 <= w <= 175 [1/cm].
    !!           - The id of the GPU device you'd like the library to use
    !!             (see the readme for a simple way to determine the ids of
    !!             any GPUs you may have on your system).  If you do not specify a
    !!             GPU id when initializing the library,
    !!             the library will query the system for available
    !!             GPUs.  If any are found, then the first device (device 0) will be used.
    !!             If you wish to force the code to run on your host CPU, pass in a value
    !!             of -1.
    !!           - If you wish to run with either the water vapor
    !!             or ozone continua, you must provide a path to a directory
    !!             containing the necessary input files.  The required input
    !!             files are included with this library in directories named
    !!             "water_vapor_continuum" and "ozone_continuum" respectively,
    !!             and are located in the base directory of this repository.
    !!        -# Calculate the optical depth for each layer in the column at each
    !!           spectral grid point by calling the @ref grt_calculate_optical_depth_fhl
    !!           function.
    !!           @attention All input arrays must be contiguous.  In addition,
    !!               the number of elements in the input pressure [atm] and
    !!               temperature [K] arrays must be equal to the number of
    !!               atmospheric levels.  The input abundance array [ppmv]
    !!               must be two-dimensional, and layed out in memory as
    !!               (level,molecule), where the molecules must be in the
    !!               same order as they were in the array of HITRAN file paths
    !!               that was passed into the @ref grt_context_init_fhl
    !!               routine.  The input optical depth array must also be
    !!               two-dimensional, and layed out in memory as (wavenumber,
    !!               layers), where the number of atmospheric layers must
    !!               equal the number of atmospheric levels minus one and the
    !!               size of the wavenumber dimension must equal the number
    !!               of spectral grid points (returned by the
    !!               @ref grt_get_spectral_grid_size_fhl function).  If any of
    !!               these arrays are not contiguous or have an incorrect size,
    !!               the behavior is undefined.
    !!
    !!        -# Release the memory allocated by the libray by calling the
    !!           @ref grt_context_free_f function.
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
    !!     @note Each atmospheric level corresponds to an interface
    !!         between adjacent atmospheric layers or the lower/upper
    !!         edge of the atmosphere.  Thus, the number of atmospheric
    !!         levels = the number of atmospheric layers plus one.
    !!         Since at least one atmospheric layer is required,
    !!         the number of atmospheric levels must be greater than or
    !!
    !!@section Limitations
    !!    This API aims to simply the function interfaces by managing
    !!    the library context internally.  Because of this, only a single
    !!    library context is used, and thus:
    !!    - The library is no longer thread-safe.  Calling these functions
    !!      in OpenMP threaded regions will lead to undefined behavior.
    !!    - The library is restricted to only using a single GPU.
    !!
    !!    If greater control over parallelism is required, please use the
    !!    @ref lowlevelfortranapi.
    !!@section Example
    !!Here is a simple example demonstrating how to use this library.
    !!@include example_fhl.F90
    !!In order to build this code, copy this code into a file and
    !!(assuming you have gcc installed), run:
    !!
    !!    $ gfortran <file> -o example.x -I<path to library include directory> \
    !!          -L<path to library lib directory> -lmolecular_lines
    !!
    !!To run this example on your GPU, make sure that you have compiled
    !!the library using the NVCC compiler (i.e., by using the provided
    !!Makefile.nvcc).


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
                                                               !! line center.  Defaults to 25 [1/cm].
            integer(kind=c_int),intent(in),optional :: gpu_id !< Id of the GPU that will be associated
                                                              !! with this context.  If not passed in,
                                                              !! then use GPU 0 if at least one GPU
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
                                                                !! input files.  If not passed in,
                                                                !! then the water vapor continuum
                                                                !! is not included in the optical
                                                                !! depth calculation.
            character(len=*),intent(in),optional :: o3_ctm_dir !< Directory containing the
                                                               !! provided ozone continuum
                                                               !! input files.  If not passed in,
                                                               !! then the ozone continuum is not
                                                               !! included in the optical depth
                                                               !! calculation.

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
