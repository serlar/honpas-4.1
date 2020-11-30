! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module kpoint_pdos_m
  !
  ! Contains data structures and routines to deal with the kpoint-grid
  ! only for the PDOS kgrid.
  ! Other uses (bands, optical, polarization) have their own structures.
  !
  use precision, only : dp
  
  ! The k-point-type
  use kpoint_t_m

  implicit none

  public :: setup_kpoint_pdos
  public :: reset_kpoint_pdos
  public :: kpoint_pdos
  public :: gamma_pdos

  private

  logical, save :: gamma_pdos
  type(kpoint_t), save :: kpoint_pdos

contains

  subroutine setup_kpoint_pdos( ucell )
    use parallel, only: Node
    use siesta_options, only: writek
    use m_spin, only: TrSym
    use kpoint_dos_m, only: setup_kpoint_dos, kpoints_dos, gamma_dos

    real(dp), intent(in) :: ucell(3,3)

    ! First try and read the k-points
    call kpoint_read(kpoint_pdos, 'PDOS', ucell, TrSym)

    if ( kpoint_pdos%method == K_METHOD_NONE ) then
      ! Default to the DOS related quantity

      if ( kpoints_DOS%N == 0 ) &
          call setup_kpoint_dos( ucell )

      call kpoint_associate(kpoint_pdos, kpoints_dos)
      gamma_pdos = gamma_dos

    else

      gamma_pdos = (kpoint_pdos%N == 1 .and. &
          dot_product(kpoint_pdos%k(:,1),kpoint_pdos%k(:,1)) < 1.0e-20_dp)

    end if
    
    ! Quick-return if non-IO
    if ( Node /= 0 ) return
    
    ! Write to XML file
    call kpoint_write_stdout(kpoint_pdos, writek, 'PDOS')
    call kpoint_write_xml(kpoint_pdos, 'PDOS')
    call kpoint_write_file(kpoint_pdos, 'PDOS.KP')
      
  end subroutine setup_kpoint_pdos

    subroutine reset_kpoint_pdos()
    use kpoint_dos_m, only: kpoints_dos

    if ( kpoint_associated(kpoint_pdos, kpoints_DOS) ) then
      call kpoint_nullify(kpoint_pdos)
    else
      call kpoint_delete(kpoint_pdos)
    end if
    gamma_PDOS = .true.
    
  end subroutine reset_kpoint_pdos
  
end module kpoint_pdos_m
