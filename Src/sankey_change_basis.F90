!
! Copyright (C) 1996-2016 The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
MODULE m_sankey_change_basis 

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sankey_change_basis

  CONTAINS
! *******************************************************************

     SUBROUTINE sankey_change_basis ( istpmove )

!******************************************************************************
! This subroutine transforms the TDKS wavefunctions into new-basis during
! Ehrenfest dynamics within TDDFT after atomic positions are changed. The 
! TDKS are transformed according to the transformation proposed by Tomfohr and 
! Sankey in Ref. phys. stat. sol. (b) 226 No.1, 115-123 (2001). 
! 
! After transforming TDKS states it calculates the density matrix in new basis set. 
! We use MatrixSwitch to manipulate matrices.
!
! This subroutine is roughly based on D. Sanchez-Portal's chgbasis subroutine from, 
! now absolete, serial version of tddft-siesta.
!
! Re-written for parallelization by Adiran Garaizar and Rafi Ullah, October 2015
! Modified by Rafi Ullah, February 2017
!********************************************************************************

  use precision
  use parallel,            only : Node, Nodes,BlockSize, IONode 
  use parallelsubs,        only : GlobalToLocalOrb, GetNodeOrbs,               &
                                  LocalToGlobalOrb
  use m_diag_option,       only : ParallelOverK
  use fdf
  use alloc
  use sys, only: die
#ifdef MPI
  use mpi_siesta,          only : mpi_bcast, mpi_comm_world, mpi_logical
#endif
  use m_spin,              only: nspin
  use kpoint_scf_m,        only: kpoint_scf, gamma_scf
  use atomlist,            only: no_u, indxuo
  use wavefunctions
  use sparse_matrices,     only : numh, listhptr, listh, S, xijo, Dscf
  use MatrixSwitch
  use matdiagon,           only: geteigen 
  !
  implicit none
  !
  integer, intent(in)     :: istpmove
  !
#ifdef MPI
  integer                 :: MPIerror,desch(9)
  external                :: diagkp
#endif
  !
  logical, save           :: frstme = .true.
  integer                 :: io, iuo, juo, nuo, jo, ind, ispin
  integer                 :: ik, j,ierror
  real(dp)                :: skxij,ckxij, kxij, qe 
  complex(dp)             :: cvar1

  type(matrix)                     :: Maux,invsqS,phi
  type(matrix)                     :: Sauxms
  type(matrix),allocatable,save    :: sqrtS(:)
  character(3)                     :: m_operation
  character(5)                     :: m_storage
  !
#ifdef DEBUG
  call write_debug( '    PRE sankey_change_basis' )
#endif
  !
#ifdef MPI
  call GetNodeOrbs(no_u,Node,Nodes,nuo)
  if (frstme) then
    if(ParallelOverK) then
      call die ( "chgbasis: tddft-siesta not parallelized over k-points." )
    endif
      call ms_scalapack_setup(mpi_comm_world,1,'c',BlockSize)
  endif
#else
  Node = 0
  Nodes = 1
  nuo = no_u
#endif
  call timer( 'sankey_change_basis', 1 )
  !
#ifdef MPI
  m_storage='pzdbc'
  m_operation='lap'
#else
  m_storage='szden'
  m_operation='lap'
#endif
  !
  IF (nspin .eq. 4) THEN
    call die ('chgbasis: ERROR: EID not yet prepared for non-collinear spin')
  END IF
    ! Allocate local arrays
  if(frstme) then
    allocate(sqrtS(kpoint_scf%N))
    do ik=1,kpoint_scf%N
      call m_allocate(sqrtS(ik),no_u,no_u,m_storage)
    end do
    frstme=.false.
  endif
  call m_allocate(Maux,no_u,no_u,m_storage)
  call m_allocate(invsqS,no_u,no_u,m_storage)
  !
  do ik = 1,kpoint_scf%N
  call m_allocate(Sauxms,no_u,no_u,m_storage)
    do iuo = 1,nuo
      call LocalToGlobalOrb(iuo, Node, Nodes, io)
      do j = 1,numh(iuo)
        ind = listhptr(iuo) + j
        juo = listh(ind)
        jo  = indxuo (juo)
        if(.not.gamma_scf) then 
          kxij = kpoint_scf%k(1,ik) * xijo(1,ind) +&
          kpoint_scf%k(2,ik) * xijo(2,ind) +&
          kpoint_scf%k(3,ik) * xijo(3,ind)
          ckxij = cos(kxij)
          skxij = -sin(kxij)
        else 
          ckxij=1.0_dp
          skxij=0.0_dp
        endif
        cvar1 =  cmplx(S(ind)*ckxij,S(ind)*skxij,dp)
        call m_set_element(Sauxms, jo, io, cvar1, complx_1, m_operation)
      enddo
    enddo
    !
    if(istpmove.eq.1) then   ! istpmove 
      ! If first step calculate S0^1/2 and save for next step. 
      call calculatesqrtS(Sauxms,invsqS,sqrtS(ik),nuo,m_storage,m_operation)
    elseif(istpmove.gt.1) then 
      ! Calculate both Sn^1/2 and Sn^-1/2 where Sn^1/2 is used in n+1 step. 
      call calculatesqrtS(Sauxms,invsqS,Maux,nuo,m_storage,m_operation)
      !Saux= Sn-1^1/2*Sn^-1/2
      call mm_multiply(invsqS,'n',sqrtS(ik),'n',&
      Sauxms,cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
      ! Passing Sn^1/2 from Maux to sqrtS for next step usage.
      call m_add ( Maux,'n',sqrtS(ik),cmplx(1.0_dp,0.0_dp,dp),              &
                  cmplx(0.0_dp,0.0_dp,dp),m_operation )
      ! C1=S0^1/2*S1^1/2*C0
      qe=2.0d0*kpoint_scf%w(ik)/dble(nspin)
      do ispin=1,nspin
        ! Cn = Saux*Cn-1 where Saux= Sn-1^1/2*Sn^-1/2
        call m_allocate ( phi,wavef_ms(ik,ispin)%dim1,                      &
                          wavef_ms(ik,ispin)%dim2,m_storage )
        call m_add( wavef_ms(ik,ispin),'n',phi,cmplx(1.0,0.0,dp),           &
                   cmplx(0.0_dp,0.0_dp,dp),m_operation )
        call mm_multiply(Sauxms,'n',phi,'n',               &
             wavef_ms(ik,ispin),cmplx(1.0_dp,0.0_dp,dp),     &
             cmplx(0.0_dp,0.0_dp,dp),m_operation)
        call m_deallocate(phi)
      enddo  
    endif   !istpmove 
  call m_deallocate(Sauxms)
  enddo          ! ik 
  !
  IF(istpmove.gt.1) THEN   ! istpmove 
    IF (IONode) THEN
      WRITE(*,'(a)') 'chgbasis: Computing DM in new basis'
    END IF
      call compute_tddm (Dscf)
  END IF
  call m_deallocate(Maux)
  call m_deallocate(invsqS)

  call timer('sankey_change_basis',2)
#ifdef DEBUG
  call write_debug( '    POS sankey_change_basis' )
#endif
  END SUBROUTINE sankey_change_basis 

  SUBROUTINE calculatesqrtS(S,invsqS,sqrtS,nu,m_storage,m_operation)
  
    use precision 
    use matdiagon
    use MatrixSwitch
    use parallelsubs,          only: LocalToGlobalOrb
    use parallel,              only: Node, Nodes
    use wavefunctions,         only: complx_0
    ! 
    implicit none
    ! 
    character(5), intent(in)                      :: m_storage
    character(3), intent(in)                      :: m_operation
    type(matrix), intent(inout)                   :: S,invsqS,sqrtS
    type(matrix)                                  :: SD01, SD02
    complex(dp)                                   :: varaux
    real(dp)                                      :: eig01, eig02
    real(dp), allocatable                         :: eigen(:)
    integer                                       :: no,nu, info, i, j,jo
    real(dp)  tiny
    data tiny  /1.0d-10/
    ! 
    no=S%dim1 
    allocate(eigen(no))
    call m_allocate(SD01,no,no,m_storage)
    call m_allocate(SD02,no,no,m_storage)
    ! Takes overlap matrix S in dense form and returns its eigenvalues
    ! in eigen(*) and eigenvectors in S.
    call geteigen(S,eigen,m_operation)
    !
    do j=1,nu
      call LocalToGlobalOrb(j,Node,Nodes,jo)
      eig01=dsqrt(dabs(eigen(jo)))
      eig02=1.0d0/(eig01+tiny)
      do i=1,no
        varaux = S%zval(i,j)
        call m_set_element(SD01,i,jo,eig01*varaux,complx_0,m_operation)
        call m_set_element(SD02,i,jo,eig02*varaux,complx_0,m_operation)
      enddo
    enddo 
    !
    deallocate(eigen)
    
    call mm_multiply(SD01,'n',S,'c',sqrtS,cmplx(1.0_dp,0.0_dp,dp),&
                     cmplx(0.0_dp,0.0_dp,dp),m_operation)
    call mm_multiply(SD02,'n',S,'c',invsqS,cmplx(1.0_dp,0.0_dp,dp),&
                     cmplx(0.0_dp,0.0_dp,dp),m_operation)
    call m_deallocate(SD01)
    call m_deallocate(SD02)
    !
  END SUBROUTINE calculatesqrtS
END MODULE m_sankey_change_basis
