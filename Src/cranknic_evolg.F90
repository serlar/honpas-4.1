! ---
! Copyright (C) 1996-2016 The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

MODULE cranknic_evolg
 
  
  USE precision
  USE units,                 ONLY: eV, Ryd_time
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cn_evolg, Uphi
  complex(kind=dp), parameter :: cZERO = cmplx(0._dp, 0._dp, dp)
  complex(kind=dp), parameter :: cONE = cmplx(1._dp, 0._dp, dp)

CONTAINS

SUBROUTINE cn_evolg ( delt )
! ********************************************************************
! Subroutine to evolve TDKS wavefunctions using Crank-Nicolson method
! and to calculate the new Density Matrix from evolved states (including
! spin polarization). Gamma-point version.
! Written by D. Sanchez-Portal, November 2002-March 2003
! Rewritten by Rafi Ullah and Adiran Garaizar June-October 2015
! Making it parallel using Matrix Switch.
! ******************** OUTPUT **************************************
! real*8 Dnew(maxnd,nspin)    : Output New Density Matrix
! real*8 eo(maxo,nspin,1)       : Output instantaneous eigenvalues
!                              (only calculated if explicitly required
!                               by user and in the last MD step)
! New wavefunctions are calculated and stored for the next step.
! *********************************************************************
!
!  Modules
!
      use sys
      use parallel    
      use parallelsubs,          only: LocalToGlobalOrb, GetNodeOrbs
      use alloc
      use wavefunctions
      use MatrixSwitch
      use siesta_options,        only: eigen_time, ntded, extrapol_H_tdks, ntded_sub
      use sparse_matrices,       only: H, S, numh, listh, listhptr, Dscf 
      use m_eo,                  only: eo
      use m_steps,               only: fincoor, final
      use m_spin,                only: nspin
      use atomlist,              only: no_s, no_l, no_u
      use m_energies,            only: etot
#ifdef MPI
      use mpi_siesta
#endif
      !
      implicit none
      !
      real(dp)             :: delt
      !
      type(matrix)         :: Hauxms,Sauxms, wfaux1, wfaux2
      complex(dp)          :: cvar1, cvar2

#ifdef MPI
      character(len=5), parameter :: m_storage = 'pzdbc'
      character(len=3), parameter :: m_operation = 'lap'
#else
      character(len=5), parameter :: m_storage = 'szden'
      character(len=5), parameter :: m_operation = 'lap'
#endif

      ! 
      integer              :: i, j, io, jo, ie, ispin, ind, nocc 
      !
      real(dp)             :: eigv
      logical, save        :: firsttime = .true.

#ifdef DEBUG
      call write_debug( '    PRE cn_evolg' )
#endif
      !
      call timer( 'cn_evolg', 1 )
      !
      call m_allocate( Sauxms,no_u,no_u,m_storage)
      !
      do ispin = 1,nspin
        nocc = wavef_ms(1,ispin)%dim2
        call m_allocate(wfaux1,no_u,nocc,m_storage)
        call m_allocate(wfaux2,nocc,nocc,m_storage)
        call m_allocate( Hauxms,no_u,no_u,m_storage)
        ! H, S from sparse to dense matrix swtich distribution
        do i = 1,no_l
          call LocalToGlobalOrb (i,Node, Nodes, io)
          do j = 1,numh(i)
            ind = listhptr(i) + j
            jo = listh(ind)
            cvar1 = cmplx(H(ind,ispin),0.0,dp)
             call m_set_element(Hauxms, jo, io, cvar1,complx_0, m_operation)
             if (ispin .eq. 1 ) then
               cvar2 = cmplx(S(ind),0.0, dp)
               call m_set_element(Sauxms, jo, io, cvar2, complx_0,m_operation)
             end if
          enddo
        enddo
        ! Evolve wavefunctions.............................................
        call evol1new(Hauxms, Sauxms, no_u, ispin,              &
                      delt,extrapol_H_tdks,ntded_sub)
        ! The following is now parallel and working as it should.
        if (eigen_time) then 
          call mm_multiply(Hauxms,'n',wavef_ms(1,ispin),'n',wfaux1,cONE,cZERO,m_operation)
          call mm_multiply(wavef_ms(1,ispin),'c',wfaux1,'n',wfaux2,cONE,cZERO,m_operation)
          DO io=1,wavef_ms(1,ispin)%dim2
            eo(io,ispin,1)= real(wfaux2%zval(io,io)) + aimag(wfaux2%zval(io,io))
          END DO
        endif
        !
        call m_deallocate(Hauxms)
        call m_deallocate(wfaux1)
        call m_deallocate(wfaux2)
      enddo ! ispin

      IF (firsttime) THEN
        firsttime = .false.
      ELSE
      IF (IOnode)  write(6,"(/a,f14.4)")  'siesta: E_KS(eV) =        ', Etot/eV 
      END IF
      ! Calculating denisty matrix.
      IF (IOnode) WRITE(*,'(/a)') 'cn_evolg: Computing DM after evolving TDKS wavefunctions'
      call compute_tddm(Dscf)
      call m_deallocate(Sauxms)
      !
      call timer( 'cn_evolg', 2 )
#ifdef DEBUG
      call write_debug( '    POS cn_evolg' )
#endif

END SUBROUTINE cn_evolg
!---------------------------------------------------------------------------------!
    subroutine Uphi(H, S, phi, no, deltat)
!*************************************************************************
!Subroutine that calculates the new wavefunction, given the old 
!wavefunction by using for the time evolution. Gamma-point 
!version. Written by A. Tsolakidis, May 2000 
!Modified by D. Sanchez-Portal, July 2002
!Modified by D. Sanchez-Portal,  2008. 
!This version is limited to first order expansion
!and avoids the inversion of the overlap.
!Modified by Rafi Ullah, October, 2015.
!Parallelized using Matrix Swtich.
!*************************** INPUT ***************************************
!integer no                  : Number of basis orbitals
!integer nocc                : Number of occupied wavefunctions per spin
!real*8 H(no,no)             : Hamiltonian matrix
!real*8 S(no,no)             : Overlap matrix
!complex*16  Phi(no,nocc)    : Old wavefunctions
!*************************** OUTPUT *************************************
!complex*16  Phi(no,nocc)    : New wavefunctions
!*************************************************************************
! Modules
!
      use MatrixSwitch 
      use matswinversion,   only: getinverse
      use m_inversemm,       only: inversemm
      use siesta_options,    only: td_inverse_linear
#ifdef MPI
      use mpi_siesta
#endif
      use parallel
!*************************************************************************   
  
 implicit none 

 integer               :: no
 real(kind=dp)         :: deltat
 type(matrix)          :: H,S,phi
 
 ! Internal variables 
 type(matrix)          :: LHS, RHS
 complex(kind=dp)      :: alpha

 
#ifdef MPI
 character(len=5), parameter :: m_storage = 'pzdbc'
 character(len=3), parameter :: m_operation = 'lap'
#else
 character(len=5), parameter :: m_storage = 'szden'
 character(len=5), parameter :: m_operation = 'lap'
#endif
 
 ! First order expansion for the evolution operator
 alpha = -0.5_dp * cmplx(0.0_dp,1.0_dp,dp) * deltat

 ! Allocate work arrays
 call m_allocate(LHS,no,no,m_storage)
 call m_allocate(RHS,phi%dim1,phi%dim2,m_storage)

 ! Setup S - alpha * H
 call m_add(S, 'n', LHS, cone, cZERO, m_operation)
 call m_add(H, 'n', LHS, alpha, cONE, m_operation)
 
 ! Calculate
 !   (S - alpha * H) psi
 call mm_multiply(LHS, 'n', phi, 'n', RHS, cONE, cZERO, m_operation)

 ! Setup S + alpha * H
 call m_add(S, 'n', LHS, cONE, cZERO, m_operation)
 call m_add(H, 'n', LHS, -alpha, cONE, m_operation)

 !---------------------------------------------------------------------------------------!
 ! There are two ways to compute inverse. One is to first obtain
 ! LU-factorization using pzgetrf and then the inverse using pzgetri.
 ! The other ways is to solve AX = B system using pzgesv. The advantage of this
 ! is way is it does inverse and a matrix multiplication at the same time. Which
 ! is desirable and in the first method is done separately. However, it is not
 ! very clear which way offers better performance when the system size is not
 ! huge. This should be thoroughly tested. For the time being the second method is
 ! hard-wired by the inversemm_linear flag, while keeping first one for testing
 ! etc.
 !----------------------------------------------------------------------------------------!
 if ( td_inverse_linear ) then
   ! Calculating (S + alpha * H)^-1 * (S - alpha * H) * phi 
   call inversemm(LHS, RHS)
   ! Copying phi_evolved back to phi
   call m_add(RHS, 'n', phi, cONE, cZERO, m_operation)
 else
   ! Calculating inverse of (S + alpha * H)
   call getinverse(LHS)
   ! (S + alpha * H)^-1 * (S - alpha * H) * phi 
   call mm_multiply(LHS, 'n', RHS, 'n', phi, cONE, cZERO, m_operation)
 endif

 call m_deallocate(LHS)
 call m_deallocate(RHS)

 END SUBROUTINE Uphi
 !------------------------------------------------------------------------------------!
 SUBROUTINE evol1new(Hauxms, Sauxms, no, ispin,         & 
             delt, extrapol, nstp)
!*************************************************************************
!Subroutine that calculates the new wavefunction, given the old 
!wavefunction by using the formula for the time evolution. Gamma-point 
!version. Written by A. Tsolakidis, May 2000 
!Modified by D. Sanchez-Portal, July 2002
!This version is limited to first order expansion
!and avoids the inversion of the overlap
!Re-written by Rafi Ullah and Adiran Garaizar June-October 2015.
!Making it parallel using Matrix Swtich.
!*************************************************************************
! Modules
!
      use fdf
      use wavefunctions,         only: wavef_ms
      use m_spin,                only: nspin
      use parallel
      use MatrixSwitch
#ifdef MPI
      use mpi_siesta
#endif
!*************************************************************************   
  implicit none 
  !
  integer                 :: no,  ispin,  nstp
  real(dp)                :: deltat, delt, varaux
  !
  type(matrix),intent(in)         :: Hauxms, Sauxms 
  type(matrix),allocatable,save   :: Hsve(:)
#ifdef MPI
  character(len=5), parameter :: m_storage = 'pzdbc'
  character(len=3), parameter :: m_operation = 'lap'
#else
  character(len=5), parameter :: m_storage = 'szden'
  character(len=5), parameter :: m_operation = 'lap'
#endif
  logical                         :: extrapol
  ! Internal variables ...
  integer                :: i, l
  complex(dp)            ::  alpha
  logical, save          :: fsttim(2) = (/.true. , .true./)
  logical, save          :: frsttime = .true.
  save                   ::  deltat

  if (frsttime) then
    !nstp is the number of "substeps" in the electronic evolution
    !the evolution operator is applied in each substep although
    !an extrapolated Hamiltonian is used "rather" than 
    !a SCF Hamiltonian
    deltat=(delt * Ryd_time) / dble(nstp)
    if (Node.eq.0) then
       write(6,'(/a,f16.6)') 'cn_evolg: TDED time step (fs)      = ',delt
       if (extrapol) then
          write(6,'(a,f16.6)') 'cn_evolg: TDED time sub-step (fs)  = ',delt/nstp
       end if
    end if
    IF (extrapol) THEN
      allocate(Hsve(nspin))
      do i=1, nspin
        call m_allocate(Hsve(i),no,no,m_storage)
      end do
    END IF
      frsttime=.false.
  endif     ! frsttime
  !
  IF (extrapol) THEN
    do l=1,nstp
      if(fsttim(ispin))then
        call Uphi(Hauxms, Sauxms, wavef_ms(1,ispin), no, deltat)
      else
        varaux=(l-0.5_dp)/dble(nstp)
        call m_add(Hauxms,'n',Hsve(ispin),cONE,cmplx(-1.0,0.0,dp),m_operation)
        call m_add(Hauxms,'n',Hsve(ispin),cONE,cmplx(varaux,0.0,dp),m_operation) 
        call Uphi(Hsve(ispin), Sauxms, wavef_ms(1,ispin),no, deltat)
      endif
    enddo
    !Storing Hamitonian for extrapolation and later correction    
    call m_add(Hauxms,'n',Hsve(ispin),cONE,cZERO,m_operation)
    fsttim(ispin)=.false.
  ELSE
    call Uphi(Hauxms, Sauxms, wavef_ms(1,ispin), no, delt*Ryd_time)
  END IF
  !
  END SUBROUTINE evol1new
!---------------------------------------------------------------------------------------!

END MODULE cranknic_evolg
