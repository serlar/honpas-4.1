! *** Module: gaufre_orbital ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Management of the information on Natural Atomic Orbitals
!!
!! This module stores data and parameters about Natural Atomic Orbitals
!! to be fitted by Gaussians.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 09.2018 Created [Yann Pouillon]
! *****************************************************************************
module gaufre_orbital

  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

  use gaufre_common

  implicit none

  private

  public :: &
&   gaufre_orbital_t

  ! Orbital quantum numbers and parameters
  type :: gaufre_orbital_t
    character(len=3) :: species = "   "
    integer :: qn_n = -1
    integer :: qn_l = -1
    integer :: zeta = -1
    logical :: polarized = .false.
    real(gfdp) :: population = -1.0_gfdp
  end type

end module gaufre_orbital
