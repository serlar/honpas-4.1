! ---
! Copyright (C) 1996-2016 The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! --
MODULE cranknic_evolk

  USE PRECISION,              ONLY: dp
  USE parallel,               ONLY: Node, Nodes, IOnode
  USE MatrixSwitch
  USE siesta_options,         ONLY: eigen_time, extrapol_H_tdks, ntded_sub
  USE units,                  ONLY: eV, Ryd_time

  IMPLICIT NONE

  PUBLIC :: cn_evolk
  PRIVATE

CONTAINS

  SUBROUTINE cn_evolk ( delt)
! ****************************************************************************
! Subroutine to evolve TDKS wavefunctions using Crank-Nicolson method and 
! calculated the new Density Matrix from evolved states. K-point version.
! Roughly based on D. Sanchez-Portal's evolk from the serial version of TDDFT.
! Written by Rafi Ullah, April 2017 
! ****************************************************************************
  
  USE parallelsubs,          ONLY: LocalToGlobalOrb
  USE sparse_matrices,       ONLY: H, S, numh, listh, listhptr, xijo, Dscf
  USE atomlist,              ONLY: no_u, no_l, indxuo
  USE m_spin,                ONLY: nspin
  USE kpoint_scf_m,          ONLY: kpoint_scf
  USE wavefunctions,         ONLY: compute_tddm, wavef_ms, complx_0, complx_1
  USE m_energies,            ONLY: etot 
  USE m_eo,                  ONLY: qo, eo

  IMPLICIT NONE

  REAL(dp)                 :: delt

  TYPE(matrix)             :: Hauxms, Sauxms, wfaux1, wfaux2

#ifdef MPI
  character(len=5), parameter :: m_storage = 'pzdbc'
  character(len=3), parameter :: m_operation = 'lap'
#else
  character(len=5), parameter :: m_storage = 'szden'
  character(len=5), parameter :: m_operation = 'lap'
#endif

  COMPLEX(dp)              :: cvar1, cvar2

  REAL(dp)                 :: kxij, ckxij, skxij
  INTEGER                  :: ik, ispin, i, j, io, jo, ind, juo, nocc
  LOGICAL, SAVE            :: firstime = .true.

#ifdef DEBUG
  call write_debug( '    PRE cn_evolk' )
#endif
  !
  call timer( 'cn_evolk', 1)
  ! 
  !
  DO ik = 1,kpoint_scf%N
    DO ispin =1,nspin
      call m_allocate ( Hauxms, no_u, no_u, m_storage)
      call m_allocate ( Sauxms, no_u, no_u, m_storage)
      nocc = wavef_ms(ik,ispin)%dim2
      ! H, S from sparse to dense matrix switch distribution
      DO i = 1, no_l
        call LocalToGlobalOrb(i,Node,Nodes,io)
        DO j   = 1, numh(i)
          ind  = listhptr(i) + j
          juo  = listh(ind)
          jo   = indxuo (juo)
          kxij = kpoint_scf%k(1,ik)*xijo(1,ind) + kpoint_scf%k(2,ik)*xijo(2,ind) + &
                 kpoint_scf%k(3,ik)*xijo(3,ind)
         ckxij =  cos(kxij)
         skxij = -sin(kxij)
         cvar1 = cmplx(H(ind,ispin)*ckxij,H(ind,ispin)*skxij,dp)
         cvar2 = cmplx(S(ind)*ckxij,S(ind)*skxij,dp)
         call m_set_element(Hauxms, jo, io, cvar1, complx_1, m_operation)
         call m_set_element(Sauxms, jo, io, cvar2, complx_1, m_operation)
        END DO
      END DO
      !
      call evol2new ( Hauxms, Sauxms, ik, ispin, delt )
      !
      IF (eigen_time) THEN
        call m_allocate(wfaux1,no_u, nocc, m_storage)
        call m_allocate(wfaux2,nocc,nocc,m_storage)
        call mm_multiply (Hauxms,'n',wavef_ms(ik,ispin),'n',wfaux1,        &
                          complx_1,complx_0,m_operation)
        call mm_multiply (wavef_ms(ik,ispin),'c',wfaux1,'n',wfaux2,        &
                          complx_1,complx_0,m_operation)
        DO io =1, nocc
          eo(io,ispin,ik) = real(wfaux2%zval(io,io)) + aimag(wfaux2%zval(io,io))
        END DO
        call m_deallocate(wfaux1)
        call m_deallocate(wfaux2)
      END IF
      !
      call m_deallocate(Hauxms)
      call m_deallocate(Sauxms)
    END DO
  END DO
  !
  IF (firstime) THEN
    firstime = .false.
  ELSE
    IF(IOnode) WRITE(6,"(/a,f14.4)")  'siesta: E_KS(eV) =        ', Etot/eV
  END IF
  ! Calculating new density matrix
  IF(IOnode) WRITE(*,'(/a)') 'cn_evolk: Computing DM after evolving TDKS wavefunctions'  
  call compute_tddm (Dscf)

  call timer( 'cn_evolk', 2)
  END SUBROUTINE cn_evolk
  
!--------------------------------------------------------------------------------!
 SUBROUTINE evol2new (Hauxms, Sauxms, ik, ispin, delt)
!********************************************************************************!
! Subroutine that calcualtes new wavefunctions, given old wavefunctions by using !
! first-order Crank-Nicolson method.                                             !
! Roughly based on D. Sanchez-Portal's serial version.                           !
! Parallel version with MatrixSwitch written by Rafi Ullah, April 2017           !
!********************************************************************************!


 USE wavefunctions
 USE m_spin,                ONLY: nspin
 USE kpoint_scf_m,          ONLY: kpoint_scf
 USE cranknic_evolg,        ONLY: Uphi
 USE atomlist,              ONLY: no_u
 USE MatrixSwitch

 IMPLICIT NONE

 TYPE(matrix)                  :: Hauxms, Sauxms
 INTEGER                       :: i, j, l, ik, ispin
 REAL(dp)                      :: delt, deltat, rvar1
#ifdef MPI
 character(len=5), parameter :: m_storage = 'pzdbc'
 character(len=3), parameter :: m_operation = 'lap'
#else
 character(len=5), parameter :: m_storage = 'szden'
 character(len=5), parameter :: m_operation = 'lap'
#endif

 COMPLEX(dp)                   :: alpha

 LOGICAL, SAVE                                  :: firsttime  = .true.
 LOGICAL, DIMENSION (:,:), ALLOCATABLE, SAVE    :: firstimeK

 IF(firsttime) THEN 
   deltat = (delt * Ryd_time) / dble(ntded_sub)
   IF (IOnode) THEN
     WRITE(6,"(/a,f16.6)") 'cn_evolg: TDED time step (fs)       =', delt
     IF (extrapol_H_tdks) THEN
       WRITE(6,"(a,f16.6)")'cn_evolg: TDED time sub-step (fs)   =', delt/ntded_sub
     END IF
   ENDIF
 
   IF(extrapol_H_tdks) THEN
     ALLOCATE(firstimeK(kpoint_scf%N, nspin))
     ALLOCATE(Hsave(kpoint_scf%N, nspin))
     DO i=1,kpoint_scf%N
       DO j=1,nspin
         call m_allocate (Hsave(i,j),no_u, no_u, m_storage) 
       END DO
     END DO
   firstimeK(1:kpoint_scf%N,1:nspin) = .true.
   END IF
 
   firsttime = .false.
 END IF
 !
 IF (extrapol_H_tdks) THEN
   DO l=1,ntded_sub
     IF (firstimeK(ik,ispin)) THEN
       call Uphi (Hauxms, Sauxms, wavef_ms(ik,ispin), no_u, deltat)
     ELSE
       rvar1 = (l -0.5_dp)/dble(ntded_sub)
       call m_add(Hauxms,'n',Hsave(ik,ispin),complx_1,cmplx(-1.0,0.0,dp), m_operation)
       call m_add(Hauxms,'n',Hsave(ik,ispin),complx_1,cmplx(rvar1,0.0,dp),m_operation)
       call Uphi(Hsave(ik,ispin), Sauxms, wavef_ms(ik,ispin), no_u, deltat)
     END IF
   END DO       
    ! Storing Hamiltonian for extrapolation and later correction.
    call m_add(Hauxms,'n',Hsave(ik,ispin),complx_1,complx_0,m_operation)
    firstimeK(ik,ispin) = .false.
 ELSE
   call Uphi (Hauxms, Sauxms, wavef_ms(ik,ispin), no_u, delt*Ryd_time)
 END IF
 !
 END SUBROUTINE evol2new
 !
END MODULE cranknic_evolk
