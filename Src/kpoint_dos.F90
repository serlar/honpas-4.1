! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module kpoint_dos_m
  !
  ! Contains data structures and routines to deal with the kpoint-grid
  ! for all DOS related quantities, PDOS, LDOS and EIG-DOS.
  ! Other uses (bands, optical, polarization) have their own structures.
  !
  use precision, only : dp
  
  ! The k-point-type
  use kpoint_t_m

  implicit none

  public :: setup_kpoint_dos
  public :: reset_kpoint_dos
  public :: kpoints_dos
  public :: gamma_dos

  private

  logical, save :: gamma_dos
  type(kpoint_t), save :: kpoints_dos

contains

  subroutine setup_kpoint_dos( ucell )
    use parallel, only: Node
    use siesta_options, only: writek
    use m_spin, only: TrSym

    real(dp), intent(in) :: ucell(3,3)

    ! First try and read the k-points
    call kpoint_read(kpoints_dos, 'DOS', ucell, TrSym)

    if ( kpoints_dos%method == K_METHOD_NONE ) then

      ! The user hasn't specified a specific DOS k-point sampling.
      call kpoint_delete(kpoints_dos)
      call kpoint_read(kpoints_dos, '', ucell, TrSym)

    end if

    gamma_dos = (kpoints_dos%N == 1 .and. &
        dot_product(kpoints_dos%k(:,1),kpoints_dos%k(:,1)) < 1.0e-20_dp)

  end subroutine setup_kpoint_dos

  subroutine reset_kpoint_dos()

    call kpoint_delete(kpoints_dos)
    gamma_DOS = .true.
    
  end subroutine reset_kpoint_dos
  
end module kpoint_dos_m
