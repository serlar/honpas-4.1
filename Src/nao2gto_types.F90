! *** Module: nao2gto_types ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Management of NAO2GTO data
!!
!! This module defines data structures to handle NAO2GTO options and
!! prescreening tolerances of ERIs in HONPAS. The specifications are read
!! from the FDF input file.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 03.2016 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_types

  use precision, only: dp

  implicit none

  private

  !> \brief Data type to store a chained list of 4-center integral parameters
  type, public :: eri_link_type

    real(dp) :: gto_eri
    integer  :: four_index(4)
    type(eri_link_type), pointer :: next

  end type eri_link_type

  !> \brief Data type to store relevant cell parameters
  type, public :: hfx_cell_type
    integer :: nsc(3)
    real(dp) :: cell(3)
    real(dp) :: cell_r(3)
  end type hfx_cell_type


  !> \brief Data type to store Hartree-Fock exchange options that can be
  !!        read from SIESTA input files
  !!
  !! \par Default values from FDF:
  !!      - potential_type = 1 (1/r/erfc(wr)/r ...)
  !!      - omega          = 0.11
  !!      - cutoff_radius  = ???
  !!      - eps_far        = 1.0e-6  (far-field screening tolerance)
  !!      - eps_pairlist   = 1.0e-6  (build shell pair-list tolerance)
  !!      - eps_schwarz    = 1.0e-6  (Schwarz tolerance)
  !!      - eps_stored     = 1.0e-6  (stored ERIs tolerance)
  !!      - DM_trunc       = .true.  (use sparse DM to screen ERIs)
  !!      - dump_fit_data  = .false. (dump fit data to nao2gto_fit.yml)
  !!      - eri_on_disk    = .false. (ERIs on memory or Disk)
  !!      - far_field      = .true.  (far-near field screening)
  !!      - on_the_fly     = .false. (evaluate ERIs at each SCF step)
  type, public :: hfx_options_type

    integer   ::  potential_type = -1
    real(dp)  ::  omega = -1.0_dp
    real(dp)  ::  cutoff_radius = -1.0_dp
    real(dp)  ::  eps_far = -1.0_dp
    real(dp)  ::  eps_pairlist = -1.0_dp
    real(dp)  ::  eps_schwarz = -1.0_dp
    real(dp)  ::  eps_stored = -1.0_dp
    logical   ::  DM_trunc = .false.
    logical   ::  dump_fit_data = .false.
    logical   ::  eri_on_disk = .false.
    logical   ::  far_field = .false.
    logical   ::  on_the_fly = .false.

  end type hfx_options_type

  !> \brief Data type to store information about orbital pairs
  type, public :: pair_list_element_type

    integer  :: pair(2)
    integer  :: nl_index
    real(dp) :: r1(3), r2(3)
    real(dp) :: dist2

  end type

  !> \brief Data type to store a list of orbital-pair information
  type, public :: pair_list_type

    type(pair_list_element_type), dimension(:), pointer :: element
    integer :: nelement

  end type pair_list_type

  !> Data type to store information about screening coefficients
  type, public :: hfx_screen_coeff_type

    real(dp) :: x(2)

  end type hfx_screen_coeff_type

  !> \brief Data type to point NAO2GTO routines to relevant SIESTA data
  type, public :: hfx_system_type

    real(dp) :: cell(3,3)
    real(dp) :: cell_r(3,3)

    integer , pointer :: maxnh
    integer , pointer :: na
    integer , pointer :: norb
    integer , pointer :: nspin
    integer , pointer :: nua
    integer , pointer :: nuo
    integer , pointer :: nuotot
    integer , pointer :: iaorb(:)
    integer , pointer :: indxua(:)
    integer , pointer :: iphorb(:)
    integer , pointer :: isa(:)
    integer , pointer :: listh(:)
    integer , pointer :: listhptr(:)
    integer , pointer :: nsc(:)
    integer , pointer :: numh(:)
    real(dp), pointer :: xa(:,:)

  end type hfx_system_type

end module nao2gto_types
