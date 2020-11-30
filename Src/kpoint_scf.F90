! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module kpoint_scf_m
  !
  ! Contains data structures and routines to deal with the kpoint-grid
  ! for the self-consistent calculation
  ! Other uses (bands, optical, polarization) have their own structures.
  !
  use precision, only : dp
  
  ! The k-point-type
  use kpoint_t_m

  implicit none

  public :: setup_kpoint_scf
  public :: reset_kpoint_scf
  public :: kpoint_scf
  public :: gamma_scf
  
  private

  logical, save :: gamma_scf
  type(kpoint_t), save :: kpoint_scf

contains

  subroutine setup_kpoint_scf( ucell )
    use parallel, only: Node
    use siesta_options, only: writek
    use m_spin, only: TrSym

    real(dp), intent(in) :: ucell(3,3)

    call kpoint_read(kpoint_scf, '', ucell, TrSym)

    gamma_scf = (kpoint_scf%N == 1 .and. &
        dot_product(kpoint_scf%k(:,1),kpoint_scf%k(:,1)) < 1.0e-20_dp)

    ! Quick-return if non-IO
    if ( Node /= 0 ) return

    call kpoint_write_stdout(kpoint_scf, all=writek)
    call kpoint_write_xml(kpoint_scf)
    call kpoint_write_file(kpoint_scf, 'KP')

  end subroutine setup_kpoint_scf

  subroutine reset_kpoint_scf()

    call kpoint_delete(kpoint_scf)
    gamma_scf = .true.
    
  end subroutine reset_kpoint_scf

end module kpoint_scf_m
