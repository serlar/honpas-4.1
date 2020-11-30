! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module kpoint_ldos_m
  !
  ! Contains data structures and routines to deal with the kpoint-grid
  ! only for the LDOS kgrid.
  ! Other uses (bands, optical, polarization) have their own structures.
  !
  use precision, only : dp
  
  ! The k-point-type
  use kpoint_t_m

  implicit none

  public :: setup_kpoint_ldos
  public :: reset_kpoint_ldos
  public :: kpoint_ldos
  public :: gamma_ldos

  private

  logical, save :: gamma_ldos
  type(kpoint_t), save :: kpoint_ldos

contains

  subroutine setup_kpoint_ldos( ucell )
    use parallel, only: Node
    use siesta_options, only: writek
    use m_spin, only: TrSym
    use kpoint_dos_m, only: setup_kpoint_dos, kpoints_dos, gamma_dos

    real(dp), intent(in) :: ucell(3,3)

    ! First try and read the k-points
    call kpoint_read(kpoint_ldos, 'LDOS', ucell, TrSym)

    if ( kpoint_ldos%method == K_METHOD_NONE ) then
      ! Default to the DOS related quantity

      if ( kpoints_DOS%N == 0 ) &
          call setup_kpoint_dos( ucell )
      
      call kpoint_associate(kpoint_ldos, kpoints_dos)
      gamma_ldos = gamma_dos

    else

      gamma_ldos = (kpoint_ldos%N == 1 .and. &
          dot_product(kpoint_ldos%k(:,1),kpoint_ldos%k(:,1)) < 1.0e-20_dp)

    end if
    
    ! Quick-return if non-IO
    if ( Node /= 0 ) return
    
    ! Write to XML file
    call kpoint_write_stdout(kpoint_ldos, writek, 'LDOS')
    call kpoint_write_xml(kpoint_ldos, 'LDOS')
    call kpoint_write_file(kpoint_ldos, 'LDOS.KP')
      
  end subroutine setup_kpoint_ldos

  subroutine reset_kpoint_ldos()
    use kpoint_dos_m, only: kpoints_dos

    if ( kpoint_associated(kpoint_ldos, kpoints_DOS) ) then
      call kpoint_nullify(kpoint_ldos)
    else
      call kpoint_delete(kpoint_ldos)
    end if
    gamma_LDOS = .true.
    
  end subroutine reset_kpoint_ldos

  
end module kpoint_ldos_m
