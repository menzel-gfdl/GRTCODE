!> @file
!! @brief Fortran bindings for utilities.
module rs_utils
use, intrinsic :: iso_c_binding, only: c_char, c_double, c_float, c_int, c_int64_t, &
                                       c_null_char, c_null_ptr, c_ptr
implicit none
private


#ifdef SINGLE_PRECISION
integer, parameter :: fp = c_float
#else
integer, parameter :: fp = c_double
#endif
integer, parameter, public :: grtcode_success = 0
integer, parameter :: grid_struct = 0
integer, parameter :: optics_struct = 1
integer, parameter :: solar_flux_struct = 3
public :: append_null_char


!> @brief Device object.
type, public :: Device_t
  integer(kind=c_int) :: device !< Device identifier.
end type Device_t


interface create_device
  function c_create_device(device, id) &
    result(error_code) &
    bind(c, name="create_device")
    import c_int
    integer(kind=c_int), intent(inout) :: device
    integer(kind=c_int), intent(in), optional :: id
    integer(kind=c_int) :: error_code
  end function c_create_device
  module procedure f_create_device
end interface create_device
public :: create_device


!> @brief One dimensional grid.
type, public :: Grid_t
  type(c_ptr) :: grid !< Pointer to grid object.
end type Grid_t


interface create_spectral_grid
  function c_create_spectral_grid(grid, w0, wn, dw) &
    result(error_code) &
    bind(c, name="create_spectral_grid")
    import c_double, c_int, c_ptr
    type(c_ptr), value :: grid
    real(kind=c_double), intent(in), value :: w0
    real(kind=c_double), intent(in), value :: wn
    real(kind=c_double), intent(in), value :: dw
    integer(kind=c_int) :: error_code
  end function c_create_spectral_grid
  module procedure f_create_spectral_grid
end interface create_spectral_grid
public :: create_spectral_grid
public :: destroy_spectral_grid


!> @brief Gas optical properties.
type, public :: Optics_t
  type(c_ptr) :: optics !< Pointer to optics object.
end type Optics_t


interface create_optics
  function c_create_optics(optics, num_layers, grid, device) &
    result(error_code) &
    bind(c, name="create_optics")
    import c_int, c_ptr
    type(c_ptr), value :: optics
    integer(kind=c_int), intent(in), value :: num_layers
    type(c_ptr), intent(in), value :: grid
    integer(kind=c_int), intent(in) :: device
    integer(kind=c_int) :: error_code
  end function c_create_optics
  module procedure f_create_optics
end interface create_optics
public :: create_optics


interface destroy_optics
  function c_destroy_optics(optics) &
    result(error_code) &
    bind(c, name="destroy_optics")
    import c_int, c_ptr
    type(c_ptr), value :: optics
    integer(kind=c_int) :: error_code
  end function c_destroy_optics
  module procedure f_destroy_optics
end interface destroy_optics
public :: destroy_optics


interface
  function malloc_struct(p, type_) &
    result(error_code) &
    bind(c)
    import c_int, c_ptr
    type(c_ptr) :: p
    integer(kind=c_int), intent(in), value :: type_
    integer(kind=c_int) :: error_code
  end function malloc_struct
end interface
public :: malloc_struct


interface
  function free_struct(p) &
    result(error_code) &
    bind(c)
    import c_int, c_ptr
    type(c_ptr) :: p
    integer(kind=c_int) :: error_code
  end function free_struct
end interface
public :: free_struct


interface
  subroutine rs_set_verbosity(level) &
    bind(c)
    import c_int
    integer(kind=c_int), intent(in), value :: level
  end subroutine rs_set_verbosity
end interface
public :: rs_set_verbosity


interface optical_properties
  function c_optical_properties(optics, tau, omega, g) &
    result(error_code) &
    bind(c, name="optical_properties")
    import c_int, c_ptr, fp
    type(c_ptr), intent(in), value :: optics
    real(kind=fp), dimension(*), intent(inout), optional :: tau
    real(kind=fp), dimension(*), intent(inout), optional :: omega
    real(kind=fp), dimension(*), intent(inout), optional :: g
    integer(kind=c_int) :: error_code
  end function c_optical_properties
  module procedure f_optical_properties
end interface optical_properties
public :: optical_properties


interface spectral_grid_properties
  !> @brief Get the spectral grid properties.
  !! @return RS_SUCCESS or an error code.
  function c_spectral_grid_properties(grid, w0, n, dw) &
    result(return_code) &
    bind(c, name="spectral_grid_properties")
    import c_double, c_int, c_int64_t, c_ptr
    type(c_ptr), intent(in), value :: grid !< Spectral grid.
    real(kind=c_double), intent(out), optional :: w0 !< Grid lower bound.
    integer(kind=c_int64_t), intent(out), optional :: n !< Grid size.
    real(kind=c_double), intent(out), optional :: dw !< Grid spacing.
    integer(kind=c_int) :: return_code
  end function c_spectral_grid_properties
  module procedure f_spectral_grid_properties
end interface spectral_grid_properties
public :: spectral_grid_properties


interface add_optics
  !> @brief Add optical properties together.
  !! @return RS_SUCCESS or an error code.
  function c_add_optics(optics, num_optics, res) &
    result(return_code) &
    bind(c, name="add_optics")
    import c_int, c_ptr
    type(c_ptr), dimension(*) :: optics
    integer(kind=c_int), intent(in), value :: num_optics
    type(c_ptr), value :: res
    integer(kind=c_int) :: return_code
  end function c_add_optics
  module procedure f_add_optics
end interface add_optics
public :: add_optics


type, public :: SolarFlux_t
  type(c_ptr) :: solar_flux !< Pointer to solar flux object.
end type SolarFlux_t


interface create_solar_flux
  !> @brief Read in data for the solar flux.
  !! @return RS_SUCCESS or an error code.
  function c_create_solar_flux(solar_flux, grid, path) &
    result(return_code) &
    bind(c, name="create_solar_flux")
    import c_char, c_int, c_ptr
    type(c_ptr), value :: solar_flux !< Solar flux object.
    type(c_ptr), value :: grid !< Spectral grid.
    character(kind=c_char, len=1), dimension(*), intent(in) :: path !< Solar flux csv file.
    integer(kind=c_int) :: return_code
  end function c_create_solar_flux
  module procedure f_create_solar_flux
end interface create_solar_flux
public :: create_solar_flux


interface destroy_solar_flux
  !> @brief Free memory for the solar flux.
  !! @return RS_SUCCESS or an error code.
  function c_destroy_solar_flux(solar_flux) &
    result(return_code) &
    bind(c, name="destroy_solar_flux")
    import c_int, c_ptr
    type(c_ptr), value :: solar_flux !< Solar flux object.
    integer(kind=c_int) :: return_code
  end function c_destroy_solar_flux
  module procedure f_destroy_solar_flux
end interface destroy_solar_flux
public :: destroy_solar_flux


interface solar_flux_properties
  !> @brief Get the solar flux properties.
  !! @return RS_SUCCESS or an error code.
  function c_solar_flux_properties(solar, flux) &
    result(return_code) &
    bind(c, name="solar_flux_properties")
    import c_int, c_ptr, fp
    type(c_ptr), intent(in), value :: solar !< Solar flux.
    real(kind=fp), dimension(*), intent(inout) :: flux !< Flux.
    integer(kind=c_int) :: return_code
  end function c_solar_flux_properties
  module procedure f_solar_flux_properties
end interface solar_flux_properties
public :: solar_flux_properties


contains


subroutine append_null_char(str_in, array_out)
  character(kind=c_char, len=*), intent(in) :: str_in
  character(kind=c_char, len=1), dimension(:), allocatable, intent(inout) :: array_out
  integer :: i
  integer :: s
  if (allocated(array_out)) then
    deallocate(array_out)
  endif
  s = len_trim(str_in)
  allocate(array_out(s+1))
  do i = 1, s
    array_out(i) = str_in(i:i)
  enddo
  array_out(i) = c_null_char
end subroutine append_null_char


function f_create_device(device, id) &
  result(error_code)
  type(Device_t), intent(inout) :: device
  integer(kind=c_int), intent(in), optional :: id
  integer(kind=c_int) :: error_code
  error_code = c_create_device(device%device, id)
end function f_create_device


function f_create_spectral_grid(grid, w0, wn, dw) &
  result(error_code)
  type(Grid_t), intent(inout) :: grid
  real(kind=c_double), intent(in) :: w0
  real(kind=c_double), intent(in) :: wn
  real(kind=c_double), intent(in) :: dw
  integer(kind=c_int) :: error_code
  grid%grid = c_null_ptr
  error_code = malloc_struct(grid%grid, grid_struct)
  if (error_code .ne. grtcode_success) then
    return
  endif
  error_code = c_create_spectral_grid(grid%grid, w0, wn, dw)
end function f_create_spectral_grid


function destroy_spectral_grid(grid) &
  result(error_code)
  type(Grid_t), intent(inout) :: grid
  integer(kind=c_int) :: error_code
  error_code = free_struct(grid%grid)
end function destroy_spectral_grid


function f_create_optics(optics, num_layers, grid, device) &
  result(error_code)
  type(Optics_t), intent(inout) :: optics
  integer(kind=c_int), intent(in) :: num_layers
  type(Grid_t), intent(in) :: grid
  type(Device_t), intent(in) :: device
  integer(kind=c_int) :: error_code
  optics%optics = c_null_ptr
  error_code = malloc_struct(optics%optics, optics_struct)
  if (error_code .ne. grtcode_success) then
    return
  endif
  error_code = c_create_optics(optics%optics, num_layers, grid%grid, device%device)
end function f_create_optics


function f_destroy_optics(optics) &
  result(error_code)
  type(Optics_t), intent(inout) :: optics
  integer(kind=c_int) :: error_code
  error_code = c_destroy_optics(optics%optics)
  if (error_code .ne. grtcode_success) then
    return
  endif
  error_code = free_struct(optics%optics)
end function f_destroy_optics


function f_optical_properties(optics, tau, omega, g) &
  result(error_code)
  type(Optics_t), intent(in) :: optics
  real(kind=fp), dimension(:,:), intent(inout), optional :: tau
  real(kind=fp), dimension(:,:), intent(inout), optional :: omega
  real(kind=fp), dimension(:,:), intent(inout), optional :: g
  integer(kind=c_int) :: error_code
  error_code = c_optical_properties(optics%optics, tau, omega, g)
end function f_optical_properties


function f_spectral_grid_properties(grid, w0, n, dw) &
  result(return_code)
  type(Grid_t), intent(in) :: grid !< Molecular lines object.
  real(kind=c_double), intent(out), optional :: w0 !< Grid lower bound.
  integer(kind=c_int64_t), intent(out), optional :: n !< Grid size.
  real(kind=c_double), intent(out), optional :: dw !< Grid spacing.
  integer(kind=c_int) :: return_code
  return_code = c_spectral_grid_properties(grid%grid, w0, n, dw)
end function f_spectral_grid_properties


function f_add_optics(optics, res) &
  result(return_code)
  type(Optics_t), dimension(:), intent(in) :: optics
  type(Optics_t), intent(inout) :: res
  integer(kind=c_int) :: return_code
  type(c_ptr), dimension(:), allocatable :: p
  integer(kind=c_int) :: num_optics
  integer :: i
  num_optics = size(optics)
  allocate(p(num_optics))
  do i = 1, num_optics
    p(i) = optics(i)%optics
  enddo
  return_code = c_add_optics(p, num_optics, res%optics)
  deallocate(p)
end function f_add_optics


function f_create_solar_flux(solar_flux, grid, path) &
  result(return_code)
  type(SolarFlux_t), intent(inout) :: solar_flux !< Solar flux object.
  type(Grid_t), intent(in) :: grid !< Spectral grid.
  character(kind=c_char, len=*), intent(in) :: path !< Solar flux csv file.
  integer(kind=c_int) :: return_code
  character(kind=c_char, len=1), dimension(:), allocatable :: buf
  call append_null_char(path, buf)
  solar_flux%solar_flux = c_null_ptr
  return_code = malloc_struct(solar_flux%solar_flux, solar_flux_struct)
  if (return_code .ne. grtcode_success) then
    return
  endif
  return_code = c_create_solar_flux(solar_flux%solar_flux, grid%grid, buf)
  deallocate(buf)
end function f_create_solar_flux


function f_destroy_solar_flux(solar_flux) &
  result(return_code)
  type(SolarFlux_t), intent(inout) :: solar_flux !< Solar flux object.
  integer(kind=c_int) :: return_code
  return_code = c_destroy_solar_flux(solar_flux%solar_flux)
  if (return_code .ne. grtcode_success) then
    return
  endif
  return_code = free_struct(solar_flux%solar_flux)
end function f_destroy_solar_flux


function f_solar_flux_properties(solar, flux) &
  result(return_code)
  type(SolarFlux_t), intent(in) :: solar !< Solar flux.
  real(kind=fp), dimension(:), intent(inout) :: flux !< Flux.
  integer(kind=c_int) :: return_code
  return_code = c_solar_flux_properties(solar%solar_flux, flux)
end function f_solar_flux_properties


end module rs_utils
