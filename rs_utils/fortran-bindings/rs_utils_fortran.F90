!> @file
!! @brief Fortran bindings for utilities.
module rs_utils
use, intrinsic :: iso_c_binding, only: c_double, c_int, c_null_ptr, c_ptr
implicit none
private


integer, parameter, public :: grtcode_success = 0
integer, parameter :: grid_struct = 0
integer, parameter :: optics_struct = 1


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


contains


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


end module rs_utils
