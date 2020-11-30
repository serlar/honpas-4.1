! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_siesta_tddft

  implicit none
  private

  public :: siesta_tddft

contains

  subroutine siesta_tddft( istep )
#ifdef MPI
    use mpi_siesta
#endif
    use sys, only : die, bye
    use precision, only: dp
    use parallel,     only : IOnode, SIESTA_worker

    use siesta_cml
    use siesta_options

    use m_state_init
    use m_setup_hamiltonian
    use m_setup_H0
    use m_steps,         only: fincoor
    use sparse_matrices, only:H, Dscf, Escf, maxnh, numh, listhptr
    use m_eo,          only: eo
    use m_energies,    only: Etot           ! Total energy
    use atomlist,      only: no_s, no_l, no_u, indxuo
    use m_spin,        only: nspin
    use kpoint_scf_m,    only: kpoint_scf

    use m_compute_energies, only: compute_energies
    use m_state_analysis, only: state_analysis

    use alloc
    use m_initwf,         only: initwf
    use wavefunctions,    only: wavef_ms, iowavef, compute_tddm, compute_tdEdm
    use m_evolve,         only: evolve
    use m_iotddft,        only: write_tddft
    use m_overfsm,        only: overfsm
    use m_final_H_f_stress,    only: final_H_f_stress
    use m_sankey_change_basis, only: sankey_change_basis
    use m_mpi_utils, only: broadcast
    use fdf

    integer, intent(inout)  :: istep
    
    integer :: ik, ispin
    real(dp) :: G2max
#ifdef DEBUG
    call write_debug( '    PRE siesta_tddft' )
#endif

#ifdef SIESTA__PEXSI
    ! Broadcast relevant things for program logic
    ! These were set in read_options, called only by "SIESTA_workers".
    call broadcast(nscf,comm=true_MPI_Comm_World)
#endif

    if ( SIESTA_worker )  then
       ! Initialization tasks for a given geometry
       call state_init( istep )
    end if

#ifdef SIESTA__PEXSI
    if (ionode) call memory_snapshot("after state_init")
#endif

    if ( fdf_get("Sonly",.false.) ) then
       if ( SIESTA_worker ) then
          call timer( 'all', 2 )
          call timer( 'all', 3 )
       end if
       call bye("S only")
    end if

    ! The current structure of the loop tries to reproduce the
    ! historical Siesta usage. It should be made more clear.
    ! Two changes:
    !
    ! -- The number of scf iterations performed is exactly
    !    equal to the number specified (i.e., the "forces"
    !    phase is not counted as a final scf step)
    !
    ! -- At the change to a TranSiesta GF run the variable "first"
    !    is implicitly reset to "true".

    ! This call computes the non-scf part of H and initializes the
    ! real-space grid structures.  It might be better to split the two,
    ! putting the grid initialization into state_init and moving the
    ! calculation of H_0 to the body of the loop, done if first=.true.  This
    ! would suit "analysis" runs in which nscf = 0
    if ( SIESTA_worker ) call setup_H0( G2max )

#ifdef SIESTA__PEXSI
    if (ionode) call memory_snapshot("after setup_H0")
#endif

    ! The first call to change basis only calculates S^+1/2
    call sankey_change_basis ( istep )

    ! Read the initial wavefunctions.
    ! In future reading wavefunctions can be moved before
    ! changebasis and then the first DM can be computed within
    ! changesbasis. Moving up iowavef would require a call to 
    ! ms_scalapack_setup. Which is now implicit in changebasis. 

    if ( istep == 1 ) then
       allocate(wavef_ms(kpoint_scf%N,nspin))
       call iowavef('read',wavef_ms,no_u,kpoint_scf%N,nspin)
       IF (IONode) THEN
       write(6,'(a)') 'Computing DM from initial KS wavefunctions'
       END IF
         call compute_tddm(Dscf)
    end if

    ! Save the final wavefunctions for restart or analysis. 
    if(tdsavewf .and. istep .ge. fincoor) then
      ! The wavefunctions are saved after transforming into the current basis
      ! but before evolving them to future.This keeps the wavefunctions
      ! concurrent with atomic position.
      if(fincoor .gt. 1) call iowavef('write',wavef_ms,no_u,kpoint_scf%N,nspin)
    end if

    do itded = 1 , ntded ! TDED loop
       
       if(IONode) then
          write(*,'(/a)') &
               '                     ************************       '
          write(*,*)'                 TDED Step    =', itded
          write(*,'(a)') &
               '                     ************************       '
       end if

       if (tdsavewf) then
         if (fincoor .eq. 1 .and. itded .eq. ntded) then
           call iowavef('write',wavef_ms,no_u,kpoint_scf%N,nspin)
         endif
       end if

       call setup_hamiltonian( itded )
       
       call evolve(  td_dt )
       
       ! The total simulation time mainly for plotting
       totime = (istep*dt - dt) + (itded*td_dt)
       
       call compute_energies (itded)
       call write_tddft(totime, istep, itded, ntded, rstart_time, &
            etot, eo, no_u,nspin,kpoint_scf%N)
       
    end do ! TDED loop
    call compute_tdEdm (Escf)
    call final_H_f_stress(istep, 1, .false.)
    call state_analysis( istep )
#ifdef DEBUG
    call write_debug( '    POS siesta_tddft' )
#endif

  end subroutine siesta_tddft
  
end module m_siesta_tddft
