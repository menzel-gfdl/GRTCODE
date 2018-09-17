!> @file
module molecular_lines_f
    use iso_c_binding
    implicit none
    private


    !> \defgroup lowlevelfortranapi Low-level Fortran API
    !!@section Overview
    !!    Given an atmospheric column made up of at least one layer,
    !!    this code calculates the total optical depth of each layer
    !!    at each point on an input spectral grid.  To use this API,
    !!
    !!        use molecular_lines_f
    !!
    !!    and follow these steps:
    !!        -# Declare library context(s).  These variables must be
    !!           of type GrtContext_t, as in:
    !!
    !!               type(GrtContext_t) :: x
    !!
    !!           Each context is a derived type that holds
    !!           the address of a struct that contains data that is
    !!           required by the library.
    !!           @note Each context can be associated with a single GPU.  If
    !!               you wish to use multiple GPUs, you must create a context
    !!               for each GPU you wish to use and pass
    !!               in the appropriate device id when initializing the context
    !!               (i.e., by calling the @ref grt_context_init_f
    !!               function).
    !!
    !!        -# Initialize the context(s) by calling the
    !!           @ref grt_context_init_f function.  Here you must provide the
    !!           number of levels per atmospheric column and parameters
    !!           (lower bound, upper bound, and resolution) for the
    !!           spectral grid on which the optical depths will be calculated.
    !!           @attention You must call this function before calling any
    !!               other function included in this library.  Failure to do
    !!               so will result in undefined behavior.
    !!           @note Each atmospheric level corresponds to an interface
    !!               between adjacent atmospheric layers or the lower/upper
    !!               edge of the atmosphere.  Thus, the number of atmospheric
    !!               levels = the number of atmospheric layers plus one.
    !!               Since at least one atmospheric layer is required,
    !!               the number of atmospheric levels must be greater than or
    !!               equal to two.
    !!
    !!           In addition, the following optional parameters may be set (not
    !!           passing these arguments in implies using the default values):
    !!           - A cutoff value for the molecular lines.  This
    !!             value denotes how far (in terms of wavenumber) from the line
    !!             center each molecular line is calculated out to.  By default
    !!             this value is set to 25 [1/cm], so if for example a line
    !!             center lies at the wavenumber 150 [1/cm], then it contributes
    !!             to the optical depth at all spectral grid points in the range
    !!             125 <= w <= 175 [1/cm].
    !!           - The id of the GPU device you wish to associate with this
    !!             context (see the readme for a simple way to determine the ids of
    !!             any GPUs you may have on your system).  If you do not specify a
    !!             GPU id when initializing the context,
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
    !!        -# Add each molecule that you want included in the optical
    !!           depth calculation by calling the @ref grt_add_molecule_f function.
    !!           Each added molecule requires an ascii [HITRAN](http://hitran.org)
    !!           database file containing the necessary molecular line parameters.  The
    !!           format of these files must match that described in Table 1 of
    !!           [Rothman et al. 2013, Journal of Quantitative Spectroscopy & Radiative
    !!           Transfer, 130](http://dx.doi.org/10.1016/j.jqsrt.2013.07.002).
    !!           Example HITRAN database files for a select set of molecules are included
    !!           with this library in a directory labeled HITRAN_files in the base of this
    !!           repository.
    !!           @attention Ozone and water vapor continua will only be
    !!               included in the optical depth calculation if the ozone
    !!               and water vapor molecules are added.
    !!
    !!        -# Set the abundance [ppmv] of each added molecule by calling
    !!           the @ref grt_set_molecule_ppmv_f function.
    !!           @attention The input abundance array must be contiguous and
    !!               its size (number of elements) must be equal to the number
    !!               of atmospheric levels passed into the @ref grt_context_init_f
    !!               function, or else the behavior is undefined.
    !!
    !!        -# Calculate the optical depth for each layer in the column at each
    !!           spectral grid point by calling the @ref grt_calculate_optical_depth_f
    !!           function.
    !!           @attention All input arrays must be contiguous.  In addition,
    !!               the number of elements in the input pressure [mb] and
    !!               temperature [K] arrays must be equal to the number of
    !!               atmospheric levels.  The number of elements in the input
    !!               optical depth array must be equal to the number of
    !!               atmospheric layers (number of atmospheric levels minus one)
    !!               times the number of spectral grid points (returned by
    !!               the @ref grt_get_spectral_grid_size_f function).  If any of
    !!               these arrays are not contiguous or have an incorrect size,
    !!               the behavior is undefined.
    !!
    !!        -# Release the memory allocated by the context(s) by calling the
    !!           @ref grt_context_free_f function.
    !!@section Example
    !!Here is a simple example demonstrating how to use this library.
    !!@include example_f.F90
    !!In order to build this code, copy this code into a file, modify
    !!the paths to input files (as needed), and (assuming you have gcc
    !!installed) run:
    !!
    !!    $ gfortran <file> -fopenmp -o example.x -I<path to library include directory> \
    !!          -L<path to library lib directory> -lmolecular_lines \
    !!          -Wl,-rpath=<path to library lib directory>
    !!
    !!To run this example on your GPU, make sure that you have compiled
    !!the library using the NVCC compiler (i.e., by using the provided
    !!Makefile.nvcc).


    public :: grt_context_init_f
    public :: grt_context_free_f
    public :: grt_add_molecule_f
    public :: grt_set_molecule_ppmv_f
    public :: grt_calculate_optical_depth_f
    public :: grt_get_num_molecules_f
    public :: grt_get_spectral_grid_size_f
    public :: grt_errstr_f
    public :: grt_set_verbosity_f
    public :: grt_get_verbosity_f


#ifdef SINGLE_PRECISION
#define FP c_float
#else
#define FP c_double
#endif


    !> Library context
    !! @ingroup lowlevelfortranapi
    type,public :: GrtContext_t
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
            real(kind=c_double),intent(in),optional :: wcutoff
            integer(kind=c_int),intent(in),optional :: gpu_id
            integer(kind=c_int),intent(in),optional :: num_threads
            character(kind=c_char,len=1),dimension(*),intent(in) :: h2o_ctm_dir
            character(kind=c_char,len=1),dimension(*),intent(in) :: o3_ctm_dir
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
        function grt_add_molecule(context, &
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
        end function grt_add_molecule
    end interface


    interface
        function grt_set_molecule_ppmv(context, &
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
        end function grt_set_molecule_ppmv
    end interface


    interface
        function grt_calculate_optical_depth(context, &
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
        end function grt_calculate_optical_depth
    end interface


    interface
        function grt_get_num_molecules(context, &
                                       n) &
            result(return_code) &
            bind(c)
            use iso_c_binding
            implicit none
            type(c_ptr),value,intent(in) :: context
            integer(kind=c_int),intent(inout) :: n
            integer(kind=c_int) :: return_code
        end function grt_get_num_molecules
    end interface


    interface
        function grt_get_spectral_grid_size(context, &
                                            n) &
            result(return_code) &
            bind(c)
            use iso_c_binding
            implicit none
            type(c_ptr),value,intent(in) :: context
            integer(kind=c_int64_t),intent(inout) :: n
            integer(kind=c_int) :: return_code
        end function grt_get_spectral_grid_size
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


    !> @ingroup lowlevelfortranapi
    !! @brief Directly bound to the @ref grt_set_verbosity function defined
    !!        in the @ref capi
    interface
        subroutine grt_set_verbosity_f(level) &
            bind(c,name="grt_set_verbosity")
            use iso_c_binding
            implicit none
            integer(kind=c_int),value,intent(in) :: level
        end subroutine grt_set_verbosity_f
    end interface


    !> @ingroup lowlevelfortranapi
    !! @brief Directly bound to the @ref grt_set_verbosity function defined
    !!        in the @ref capi
    interface
        function grt_get_verbosity_f() &
            result(level) &
            bind(c,name="grt_get_verbosity")
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: level
        end function grt_get_verbosity_f
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
            real(kind=c_double),intent(in),optional :: wcutoff !< Cutoff [1/cm] from spectral
                                                               !! line center.  This
                                                               !! defaults to 25 [1/cm].*/
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
            integer(kind=c_int) :: return_code
            character(len=1024) :: hbuf
            character(len=1024) :: obuf
            if (present(h2o_ctm_dir)) then
                hbuf = trim(h2o_ctm_dir)//c_null_char
            else
                hbuf = "none"//c_null_char
            endif
            if (present(o3_ctm_dir)) then
                obuf = trim(o3_ctm_dir)//c_null_char
            else
                obuf = "none"//c_null_char
            endif
            return_code = grt_context_init(context%p, &
                                           num_levels, &
                                           w0, &
                                           wn, &
                                           wres, &
                                           wcutoff, &
                                           gpu_id, &
                                           num_threads, &
                                           trim(hbuf), &
                                           trim(obuf))
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
        function grt_add_molecule_f(context, &
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
                                                                                  !! computed.  Defaults to 3250 [1/cm].
            integer(kind=c_int) :: return_code

            return_code = grt_add_molecule(context%p, &
                                           trim(hitran_filepath)//c_null_char, &
                                           molecule_id, &
                                           min_line_center_wavenumber, &
                                           max_line_center_wavenumber)
        end function grt_add_molecule_f


        !> @ingroup lowlevelfortranapi
        !! @brief Update the abundances [ppmv] for a molecule.
        !! @return 0 if completed successfully, or else an error code.
        function grt_set_molecule_ppmv_f(context, &
                                         molecule_id, &
                                         ppmv) &
            result(return_code)
            type(GrtContext_t),intent(in) :: context !< Library context.
            integer(kind=c_int),intent(in) :: molecule_id !< Molecule id returned by
                                                          !! @ref grt_add_molecule_f.
            real(kind=FP),dimension(*),intent(in) :: ppmv !< Array of molecular abundances [ppmv].
                                                          !! The size of this array must be
                                                          !! eqaul to the number of atmospheric
                                                          !! levels.
            integer(kind=c_int) :: return_code

            return_code = grt_set_molecule_ppmv(context%p, &
                                                molecule_id, &
                                                ppmv)
        end function grt_set_molecule_ppmv_f


        !> @ingroup lowlevelfortranapi
        !! @brief Calculate the total optical depth of each atmospheric
        !!        layer at each spectral grid point.
        !! @return 0 if completed successfully, or else an error code.
        function grt_calculate_optical_depth_f(context, &
                                               pressure, &
                                               temperature, &
                                               optical_depth) &
            result(return_code)
            type(GrtContext_t),value,intent(in) :: context !< Library context.
            real(kind=FP),dimension(*),intent(in) :: pressure !< Array of atmospheric pressures [mb].
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

            return_code = grt_calculate_optical_depth(context%p, &
                                                      pressure, &
                                                      temperature, &
                                                      optical_depth)
        end function grt_calculate_optical_depth_f


        !> @ingroup lowlevelfortranapi
        !! @brief Get the number of molecules that have been added to
        !!        the context.
        function grt_get_num_molecules_f(context, &
                                         n) &
            result(return_code)
            type(GrtContext_t),intent(in) :: context !< Library context.
            integer(kind=c_int),intent(inout) :: n !< Number of molecules.
            integer(kind=c_int) :: return_code

            return_code = grt_get_num_molecules(context%p, &
                                                n)
        end function grt_get_num_molecules_f


        !> @ingroup lowlevelfortranapi
        !! @brief Get the spectral grid size.
        function grt_get_spectral_grid_size_f(context, &
                                              n) &
            result(return_code)
            type(GrtContext_t),intent(in) :: context !< Library context.
            integer(kind=c_int64_t),intent(inout) :: n !< Spectral grid size.
            integer(kind=c_int) :: return_code

            return_code = grt_get_spectral_grid_size(context%p, &
                                                     n)
        end function grt_get_spectral_grid_size_f


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


end module molecular_lines_f
