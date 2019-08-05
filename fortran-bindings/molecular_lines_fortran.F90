module molecular_lines
use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_null_char, &
                                       c_null_ptr, c_ptr
use rs_utils


#ifdef SINGLE_PRECISION
integer, parameter :: fp = c_float
#else
integer, parameter :: fp = c_double
#endif
integer, parameter :: molecular_lines_struct = 2
integer(kind=c_int), parameter, public :: H2O = 1
integer(kind=c_int), parameter, public :: CO2 = 2
integer(kind=c_int), parameter, public :: O3 = 3
integer(kind=c_int), parameter, public :: N2O = 4
integer(kind=c_int), parameter, public :: CO = 5
integer(kind=c_int), parameter, public :: CH4 = 6
integer(kind=c_int), parameter, public :: O2 = 7
integer(kind=c_int), parameter, public :: NO = 8
integer(kind=c_int), parameter, public :: SO2 = 9
integer(kind=c_int), parameter, public :: NO2 = 10
integer(kind=c_int), parameter, public :: NH3 = 11
integer(kind=c_int), parameter, public :: HNO3 = 12
integer(kind=c_int), parameter, public :: OH = 13
integer(kind=c_int), parameter, public :: HF = 14
integer(kind=c_int), parameter, public :: HCl = 15
integer(kind=c_int), parameter, public :: HBr = 16
integer(kind=c_int), parameter, public :: HI = 17
integer(kind=c_int), parameter, public :: ClO = 18
integer(kind=c_int), parameter, public :: OCS = 19
integer(kind=c_int), parameter, public :: H2CO = 20
integer(kind=c_int), parameter, public :: HOCl = 21
integer(kind=c_int), parameter, public :: N2 = 22
integer(kind=c_int), parameter, public :: HCN = 23
integer(kind=c_int), parameter, public :: CH3Cl = 24
integer(kind=c_int), parameter, public :: H2O2 = 25
integer(kind=c_int), parameter, public :: C2H2 = 26
integer(kind=c_int), parameter, public :: C2H6 = 27
integer(kind=c_int), parameter, public :: PH3 = 28
integer(kind=c_int), parameter, public :: COF2 = 29
integer(kind=c_int), parameter, public :: SF6_MOL = 30
integer(kind=c_int), parameter, public :: H2S = 31
integer(kind=c_int), parameter, public :: HCOOH = 32
integer(kind=c_int), parameter, public :: HO2 = 33
integer(kind=c_int), parameter, public :: O = 34
integer(kind=c_int), parameter, public :: ClONO2 = 35
integer(kind=c_int), parameter, public :: NOp = 36
integer(kind=c_int), parameter, public :: HOBr = 37
integer(kind=c_int), parameter, public :: C2H4 = 38
integer(kind=c_int), parameter, public :: CH3OH = 39
integer(kind=c_int), parameter, public :: CH3Br = 40
integer(kind=c_int), parameter, public :: CH3CN = 41
integer(kind=c_int), parameter, public :: CF4_MOL = 42
integer(kind=c_int), parameter, public :: C4H2 = 43
integer(kind=c_int), parameter, public :: HC3N = 44
integer(kind=c_int), parameter, public :: H2 = 45
integer(kind=c_int), parameter, public :: CS = 46
integer(kind=c_int), parameter, public :: SO3 = 47
integer(kind=c_int), parameter, public :: C2N2 = 48
integer(kind=c_int), parameter, public :: COCl2 = 49
integer(kind=c_int), parameter, public :: SO = 50
integer(kind=c_int), parameter, public :: C3H4 = 51
integer(kind=c_int), parameter, public :: CH3 = 52
integer(kind=c_int), parameter, public :: CS2 = 53
integer(kind=c_int), parameter, public :: MAX_NUM_MOLECULES = 53
integer(kind=c_int), parameter, public :: CFC11 = 0
integer(kind=c_int), parameter, public :: CFC12 = 1
integer(kind=c_int), parameter, public :: CFC113 = 2
integer(kind=c_int), parameter, public :: CFC114 = 3
integer(kind=c_int), parameter, public :: CFC115 = 4
integer(kind=c_int), parameter, public :: HCFC22 = 5
integer(kind=c_int), parameter, public :: HCFC141b = 6
integer(kind=c_int), parameter, public :: HCFC142b = 7
integer(kind=c_int), parameter, public :: HFC23 = 8
integer(kind=c_int), parameter, public :: HFC125 = 9
integer(kind=c_int), parameter, public :: HFC134a = 10
integer(kind=c_int), parameter, public :: HFC143a = 11
integer(kind=c_int), parameter, public :: HFC152a = 12
integer(kind=c_int), parameter, public :: HFC227ea = 13
integer(kind=c_int), parameter, public :: HFC245fa = 14
integer(kind=c_int), parameter, public :: CCl4 = 15
integer(kind=c_int), parameter, public :: C2F6 = 16
integer(kind=c_int), parameter, public :: CF4 = 17
integer(kind=c_int), parameter, public :: CH2Cl2 = 18
integer(kind=c_int), parameter, public :: NF3 = 19
integer(kind=c_int), parameter, public :: SF6 = 20
integer(kind=c_int), parameter, public :: MAX_NUM_CFCS = 21
integer(kind=c_int), parameter, public :: CIA_N2 = 0
integer(kind=c_int), parameter, public :: CIA_O2 = 1
integer(kind=c_int), parameter, public :: MAX_NUM_CIAS = 2


type, public :: MolecularLines_t
  type(c_ptr) :: ml !< Pointer to molecular lines object.
end type MolecularLines_t


interface create_molecular_lines
  !> @brief Reserve memory for molecular lines.
  !! @return RS_SUCCESS or an error code.
  function c_create_molecular_lines(ml, num_levels, grid, device, hitran_path, h2o_ctm_dir, &
                                    o3_ctm_dir, wcutoff, optical_depth_method) &
    result(return_code) &
    bind(c, name="create_molecular_lines")
    import c_char, c_double, c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: num_levels !< Number of atmospheric levels.
    type(c_ptr), intent(in), value :: grid !< Spectral grid.
    integer(kind=c_int), intent(in) :: device !< Device.
    character(kind=c_char, len=1), dimension(*), intent(in) :: hitran_path !< Path to HITRAN database file.
    character(kind=c_char, len=1), dimension(*), intent(in) :: h2o_ctm_dir !< Path to water vapor continuum directory.
    character(kind=c_char, len=1), dimension(*), intent(in) :: o3_ctm_dir !< Path to ozone continuum directory.
    real(kind=c_double), intent(in), optional :: wcutoff !< Cutoff from line center [1/cm].
    integer(kind=c_int), intent(in), optional :: optical_depth_method !< Method used to calculate the optical depths.
    integer(kind=c_int) :: return_code
  end function c_create_molecular_lines
  module procedure f_create_molecular_lines
end interface create_molecular_lines
public :: create_molecular_lines


interface destroy_molecular_lines
  !> @brief Free memory for the molecular lines.
  !! @return RS_SUCCESS or an error code.
  function c_destroy_molecular_lines(ml) &
    result(return_code) &
    bind(c, name="destroy_molecular_lines")
    import c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int) :: return_code
  end function c_destroy_molecular_lines
  module procedure f_destroy_molecular_lines
end interface destroy_molecular_lines
public :: destroy_molecular_lines


interface add_molecule
  !> @brief Add a molecule.
  !! @return RS_SUCCESS or an error code.
  function grt_add_molecule(ml, molecule_id, min_line_center, max_line_center) &
    result(return_code) &
    bind(c)
    import c_double, c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: molecule_id !< Molecule id.
    real(kind=c_double), intent(in), optional :: min_line_center !< Lower bound [1/cm] for spectral line centers.
    real(kind=c_double), intent(in), optional :: max_line_center !< Upper bound [1/cm] for spectral line centers.
    integer(kind=c_int) :: return_code
  end function grt_add_molecule
  module procedure f_add_molecule
end interface add_molecule
public :: add_molecule


interface set_molecule_ppmv
  !> @brief Update a molecule's ppmv.
  !! @return RS_SUCCESS or an error code.
  function grt_set_molecule_ppmv(ml, molecule_id, ppmv) &
    result(return_code) &
    bind(c)
    import c_int, c_ptr, fp
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: molecule_id  !< Molecule id.
    real(kind=fp), dimension(*), intent(in) :: ppmv !< Abundance [ppmv] (level).
    integer(kind=c_int) :: return_code
  end function grt_set_molecule_ppmv
  module procedure f_set_molecule_ppmv
end interface set_molecule_ppmv
public :: set_molecule_ppmv


interface add_cfc
  !> @brief Add a CFC.
  !! @return RS_SUCCESS or an error code.
  function grt_add_cfc(ml, cfc_id, filepath) &
    result(return_code) &
    bind(c)
    import c_char, c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: cfc_id !< CFC id.
    character(kind=c_char, len=1), dimension(*), intent(in) :: filepath !< Path to CFC cross section csv file.
    integer(kind=c_int) :: return_code
  end function grt_add_cfc
  module procedure f_add_cfc
end interface add_cfc
public :: add_cfc


interface set_cfc_ppmv
  !> @brief Update a CFC's ppmv.
  !! @return RS_SUCCESS or an error code.
  function grt_set_cfc_ppmv(ml, cfc_id, ppmv) &
    result(return_code) &
    bind(c)
    import c_int, c_ptr, fp
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: cfc_id !< CFC id.
    real(kind=fp), dimension(*), intent(in) :: ppmv !< Abundance [ppmv] (level).
    integer(kind=c_int) :: return_code
  end function grt_set_cfc_ppmv
  module procedure f_set_cfc_ppmv
end interface set_cfc_ppmv
public :: set_cfc_ppmv


interface add_cia
  !> @brief Activate collision-induced absorption between two species.
  !! @return RS_SUCCESS or an error code.
  function grt_add_cia(ml, species1, species2, filepath) &
    result(return_code) &
    bind(c)
    import c_char, c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: species1 !< Id of species.
    integer(kind=c_int), intent(in), value :: species2 !< Id of species.
    character(kind=c_char, len=1), dimension(*), intent(in) :: filepath !< Path to cross section csv file.
    integer(kind=c_int) :: return_code
  end function grt_add_cia
  module procedure f_add_cia
end interface add_cia
public :: add_cia


interface set_cia_ppmv
  !> @brief Update a CIA species' ppmv.
  !! @return RS_SUCCESS or an error code.
  function grt_set_cia_ppmv(ml, cia_id, ppmv) &
    result(return_code) &
    bind(c)
    import c_int, c_ptr, fp
    type(c_ptr), value :: ml !< Molecularlines object.
    integer(kind=c_int), intent(in), value :: cia_id !< CIA species id.
    real(kind=fp), dimension(*), intent(in) :: ppmv !< Abundance [ppmv] (level).
    integer(kind=c_int) :: return_code
  end function grt_set_cia_ppmv
  module procedure f_set_cia_ppmv
end interface set_cia_ppmv
public :: set_cia_ppmv


interface calculate_optics
  !> @brief Calcluate the total optical depth in each layer at each spectral grid point.
  !! @return RS_SUCCESS or an error code.
  function grt_calculate_optical_depth(ml, pressure, temperature, optics) &
    result(return_code) &
    bind(c)
    import c_int, c_ptr, fp
    type(c_ptr), value :: ml !< Molecular lines object.
    real(kind=fp), dimension(*), intent(in) :: pressure !< Pressure [mb] (level).
    real(kind=fp), dimension(*), intent(in) :: temperature !< Temperature [K] (level).
    type(c_ptr), value :: optics !< Optics object.
    integer(kind=c_int) :: return_code
  end function grt_calculate_optical_depth
  module procedure f_calculate_optics
end interface calculate_optics
public :: calculate_optics


interface num_molecules
  !> @brief Get the number of molecules.
  !! @return RS_SUCCESS or an error code.
  function grt_get_num_molecules(ml, n) &
    result(return_code) &
    bind(c)
    import c_int, c_ptr
    type(c_ptr), intent(in), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(out) :: n !< Number of molecules.
    integer(kind=c_int) :: return_code
  end function grt_get_num_molecules
  module procedure f_num_molecules
end interface num_molecules
public :: num_molecules


interface grt_errstr
  !> @brief Return a message for an input return code.
  !! @return RS_SUCCESS or an error code.
  function c_grt_errstr(code, buf, buf_size) &
    result(return_code) &
    bind(c, name="grt_errstr")
    import c_char, c_int
    integer(kind=c_int), intent(in), value :: code !< Error code.
    character(kind=c_char, len=1), dimension(*) :: buf !< Buffer to hold error message.
    integer(kind=c_int), intent(in), value :: buf_size !< Size of input buffer.
    integer(kind=c_int) :: return_code
  end function c_grt_errstr
  module procedure f_grt_errstr
end interface grt_errstr
public :: grt_errstr


interface rayleigh_scattering
  !> @brief Calculate the optical properties due to Rayleigh scattering.
  !! @return RS_SUCCESS or an error code.
  function c_rayleigh_scattering(optics, pressure) &
    result(return_code) &
    bind(c, name="rayleigh_scattering")
    import c_int, c_ptr, fp
    type(c_ptr), value :: optics !< Optics object.
    real(kind=fp), dimension(*), intent(in) :: pressure !< Pressure [mb] (level).
    integer(kind=c_int) :: return_code
  end function c_rayleigh_scattering
  module procedure f_rayleigh_scattering
end interface rayleigh_scattering
public :: rayleigh_scattering


contains


function f_create_molecular_lines(ml, num_levels, grid, device, hitran_path, h2o_ctm_dir, &
                                  o3_ctm_dir, wcutoff, optical_depth_method) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: num_levels !< Number of atmospheric levels.
  type(Grid_t), intent(in) :: grid !< Spectral grid.
  type(Device_t), intent(in) :: device !< Device.
  character(kind=c_char, len=*), intent(in) :: hitran_path !< Path to HITRAN database file.
  character(kind=c_char, len=*), intent(in), optional :: h2o_ctm_dir !< Path to water vapor continuum directory.
  character(kind=c_char, len=*), intent(in), optional :: o3_ctm_dir !< Path to ozone continuum directory.
  real(kind=c_double), intent(in), optional :: wcutoff !< Cutoff from line center [1/cm].
  integer(kind=c_int), intent(in), optional :: optical_depth_method !< Method used to calculate the optical depths.
  integer(kind=c_int) :: return_code
  character(kind=c_char, len=1), dimension(:), allocatable :: hitran
  character(kind=c_char, len=1), dimension(:), allocatable :: h2ob
  character(kind=c_char, len=1), dimension(:), allocatable :: o3b
  call append_null_char(hitran_path, hitran)
  if (present(h2o_ctm_dir)) then
    call append_null_char(h2o_ctm_dir, h2ob)
  else
    call append_null_char("none", h2ob)
  endif
  if (present(o3_ctm_dir)) then
    call append_null_char(o3_ctm_dir, o3b)
  else
    call append_null_char("none", o3b)
  endif
  ml%ml = c_null_ptr
  return_code = malloc_struct(ml%ml, molecular_lines_struct)
  if (return_code .ne. grtcode_success) then
    return
  endif
  return_code = c_create_molecular_lines(ml%ml, num_levels, grid%grid, device%device, &
                                         hitran, h2ob, o3b, wcutoff, optical_depth_method)
  deallocate(hitran)
  deallocate(h2ob)
  deallocate(o3b)
end function f_create_molecular_lines


function f_destroy_molecular_lines(ml) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int) :: return_code
  return_code = c_destroy_molecular_lines(ml%ml)
  if (return_code .ne. grtcode_success) then
    return
  endif
  return_code = free_struct(ml%ml)
end function f_destroy_molecular_lines


function f_add_molecule(ml, molecule_id, min_line_center, max_line_center) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: molecule_id !< Molecule id.
  real(kind=c_double), intent(in), optional :: min_line_center !< Lower bound [1/cm] for spectral line centers.
  real(kind=c_double), intent(in), optional :: max_line_center !< Upper bound [1/cm] for spectral line centers.
  integer(kind=c_int) :: return_code
  return_code = grt_add_molecule(ml%ml, molecule_id, min_line_center, max_line_center)
end function f_add_molecule


function f_set_molecule_ppmv(ml, molecule_id, ppmv) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: molecule_id  !< Molecule id.
  real(kind=fp), dimension(:), intent(in) :: ppmv !< Abundance [ppmv] (level).
  integer(kind=c_int) :: return_code
  return_code = grt_set_molecule_ppmv(ml%ml, molecule_id, ppmv)
end function f_set_molecule_ppmv


function f_add_cfc(ml, cfc_id, filepath) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: cfc_id !< CFC id.
  character(kind=c_char, len=*), intent(in) :: filepath !< Path to CFC cross section csv file.
  integer(kind=c_int) :: return_code
  character(kind=c_char, len=1), dimension(:), allocatable :: buf
  call append_null_char(filepath, buf)
  return_code = grt_add_cfc(ml%ml, cfc_id, buf)
  deallocate(buf)
end function f_add_cfc


function f_set_cfc_ppmv(ml, cfc_id, ppmv) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: cfc_id !< CFC id.
  real(kind=fp), dimension(:), intent(in) :: ppmv !< Abundance [ppmv] (level).
  integer(kind=c_int) :: return_code
  return_code = grt_set_cfc_ppmv(ml%ml, cfc_id, ppmv)
end function f_set_cfc_ppmv


function f_add_cia(ml, species1, species2, filepath) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: species1 !< Id of species.
  integer(kind=c_int), intent(in) :: species2 !< Id of species.
  character(kind=c_char, len=*), intent(in) :: filepath !< Path to cross section csv file.
  integer(kind=c_int) :: return_code
  character(kind=c_char, len=1), dimension(:), allocatable :: buf
  call append_null_char(filepath, buf)
  return_code = grt_add_cia(ml%ml, species1, species2, buf)
  deallocate(buf)
end function f_add_cia


function f_set_cia_ppmv(ml, cia_id, ppmv) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecularlines object.
  integer(kind=c_int), intent(in) :: cia_id !< CIA species id.
  real(kind=fp), dimension(:), intent(in) :: ppmv !< Abundance [ppmv] (level).
  integer(kind=c_int) :: return_code
  return_code = grt_set_cia_ppmv(ml%ml, cia_id, ppmv)
end function f_set_cia_ppmv


function f_calculate_optics(ml, pressure, temperature, optics) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecular lines object.
  real(kind=fp), dimension(:), intent(in) :: pressure !< Pressure [mb] (level).
  real(kind=fp), dimension(:), intent(in) :: temperature !< Temperature [K] (level).
  type(Optics_t), intent(inout) :: optics !< Optics object.
  integer(kind=c_int) :: return_code
  return_code = grt_calculate_optical_depth(ml%ml, pressure, temperature, optics%optics)
end function f_calculate_optics


function f_num_molecules(ml, n) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(out) :: n !< Number of molecules.
  integer(kind=c_int) :: return_code
  return_code = grt_get_num_molecules(ml%ml, n)
end function f_num_molecules


function f_grt_errstr(code, buf) &
  result(return_code)
  integer(kind=c_int), intent(in) :: code !< Error code.
  character(kind=c_char, len=*), intent(inout) :: buf !< Buffer to hold error message.
  integer(kind=c_int) :: return_code
  integer :: i
  logical :: eos
  return_code = c_grt_errstr(code, buf, len(buf))
  eos = .false.
  do i = 1, len(buf)
    eos = eos .or. buf(i:i) .eq. c_null_char
    if (eos) then
      buf(i:i) = " "
    endif
  enddo
end function f_grt_errstr


function f_rayleigh_scattering(optics, pressure) &
  result(return_code)
  type(Optics_t), intent(inout) :: optics !< Optics object.
  real(kind=fp), dimension(:), intent(in) :: pressure !< Pressure [mb] (level).
  integer(kind=c_int) :: return_code
  return_code = c_rayleigh_scattering(optics%optics, pressure)
end function f_rayleigh_scattering


end module molecular_lines
