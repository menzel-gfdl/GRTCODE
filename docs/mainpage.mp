/** @file/
    @mainpage GRTcode Molecular Lines

    # Overview
        This library calculates optical depth values for atmospheric columns
        on GPU and CPU architectures by explicitly treating the rotational-
        vibrational lines of atmospheric molecules.  Molecular line parameters
        are passed into the model as ascii HITRAN database files, and optical
        depth values in each layer of the atmospheric column are calculated
        on each point of a user-defined input spectral grid.

    # APIs
        A @ref capi, @ref lowlevelfortranapi, and @ref highlevelfortranapi
        are provided.  The @ref lowlevelfortranapi simply provides direct
        bindings to the @ref capi, while the @ref highlevelfortranapi
        provides a more simplified interface (at the expense of less control).
        Click the links to view more detailed information about each API.

    # Requirements
        This library requires a c and fortran compiler, such as the freely
        available gcc and gfortran.  The c compiler must support the c99
        standard, and the fortran compiler must support this [Fortran 2012
        Technical Specification] (https://www.iso.org/standard/45136.html).
        In order to run on a NVIDIA GPU,
        CUDA must also be installed and the library must be built using the
        included NVCC compiler.  This library also requires make and install.

    # Code
        The source code currently resides in this [Gitlab repository]
        (https://gitlab.gfdl.noaa.gov/Raymond.Menzel/grtcodev2), on branch
        modular_lines.  To obtain the code, run:\n\n
        ```$ git clone
             https://gitlab.gfdl.noaa.gov/Raymond.Menzel/grtcodev2.git .```\n
        ```$ git checkout modular_lines```

    # Building
    ### CPU-only
        To build using default settings (assuming you have gcc and gfortran
        installed), run:\n\n
        ```$ make```\n
        ```$ make install```\n\n  It is recommended that you also run:\n\n
        ```$ make test```\n\n
        to make sure that everything built and runs properly.  The make
        variable PREFIX can be used to specify where the library will be
        installed;  if it is not specified the library is installed in the
        current directory.  The default compilers and compiler flags can be
        overridden by specifying CC, FC, CFLAGS, and FFLAGS when running
        make in the usual fashion.  For example,\n\n
        ```$ make CC=icc CFLAGS=-O3 FC=ifort FFLAGS=-O3```\n\n
        This library takes advantage of thread-level parallelism through the
        use of OpenMP.  If your compiler supports OpenMP, make sure you
        activate it by setting the corresponding compiler option.
    ### With GPUs
        If you have CUDA and the NVCC compiler installed, the library
        can built to run on NVIDIA GPUs by running:\n\n
        ```$ make -f Makefile.nvcc```\n
        ```$ make -f Makefile.nvcc install```\n\n  Once again it is
        recommended that you also run:\n\n
        ```$ make -f Makefile.nvcc test```\n\n
        Once again the make variable PREFIX can be used to specify where the
        library will be installed (default is the current directory).
*/
