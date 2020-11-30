module m_initwf
  use precision
  use wavefunctions, only:iowavef, wavef_ms, complx_0
  use MatrixSwitch
  use sparse_matrices, only: numh, listhptr, listh, H, S, xijo
  use m_eo,            only: eo, qo
#ifdef MPI
  use mpi_siesta, only: MPI_Comm_World
#endif
  !
  implicit none
  !
  private
  !
  public :: initwf

  CONTAINS
  !
  subroutine initwf( istpp,totime)
! *********************************************************************
! Subroutine to calculate the eigenvalues and eigenvectors,
! for given Hamiltonian and Overlap matrices (including
! spin polarization), providing the initial wavefunctions 
! for a time dependent electonic simulations.
! Written by D. Sanchez-Portal, November 2002-March 2003 after
! subroutine diagon by J.Soler, P.Ordejon, and J. D. Gale (1997-1999)
! Modified by Rafi Ullah, CIC nanoGUNE, October 2015 
! to make it  work with parallel TDDFT-siesta using MatrixSwitch
! **************************** INPUT **********************************
! integer no                  : Number of basis orbitals
! integer nspin               : Spin polarization (1 or 2)
! integer maxspn              : Second dimension of eo and qo
! integer maxnh               : Maximum number of orbitals interacting  
! integer maxnd               : Maximum number of nonzero elements of 
!                               each row of density matrix
! integer maxo                : First dimension of eo and qo
! integer numh(nuo)           : Number of nonzero elements of each row 
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to each row (-1) of the
!                               hamiltonian matrix
! integer listh(maxlh)        : Nonzero hamiltonian-matrix element  
!                               column indexes for each matrix row
! integer numd(nuo)           : Number of nonzero elements of each row 
!                               of density matrix
! integer listdptr(nuo)       : Pointer to each row (-1) of the
!                               density matrix
! integer listd(maxnh)        : Nonzero density-matrix element column 
!                               indexes for each matrix row
! real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not used if only gamma point)
! integer indxuo(no)          : Index of equivalent orbital in unit cell
!                               Unit cell orbitals must be the first in
!                               orbital lists, i.e. indxuo.le.nuo, with
!                               nuo the number of orbitals in unit cell
! integer nk                  : Number of k points
! real*8  kpoint(3,nk)        : k point vectors
! real*8  wk(nk)              : k point weights (must sum one)
! integer nuotot              : total number of orbitals in unit cell 
!                               over all processors
! *************************** OUTPUT **********************************
! real*8 eo(maxo,maxspn,nk)   : Eigenvalues
! real*8 qo(maxo,maxspn,nk)   : Occupations of eigenstates
! *************************** UNITS ***********************************
! xij and kpoint must be in reciprocal coordinates of each other.
! temp and H must be in the same energy units.
! eo, Enew and ef returned in the units of H.
! *************************** Parallel ********************************
!
!  Modules
!
      use parallel,        only : Node, Nodes, BlockSize
      use parallelsubs,    only : GlobalToLocalOrb, GetNodeOrbs
      use m_diag_option,   only : ParallelOverK
      use fdf
      use densematrix,     only : Haux, Saux, psi
      use sparse_matrices, only : maxnh
      use kpoint_scf_m,    only : kpoint_scf, gamma_scf
      use atomlist,        only : no_s, no_l, no_u, qtot, indxuo
      use m_spin,          only : nspin
      use alloc
      use m_memory
      use m_fermid,      only : fermid
      use sys,           only : die
#ifdef MPI
      use mpi_siesta,   only : mpi_bcast, mpi_comm_world, &
                               mpi_logical, mpi_double_precision
#endif
      !
      implicit none
      !
      integer, intent (inout) :: istpp
      !
      real(dp), intent (inout) :: totime
      ! 
      external           :: io_assign, io_close, readsp
      !
      character          :: sname*30, fname*33, m_storage*5
      !
      logical            :: fixspin, ioifound, degen
      logical, save      :: spiral
      logical, save      :: frstme = .true.
      !
      !
      integer            :: io, iuo, iu, nhs, npsi, nuo, nocc(2), ispin,ik,     &
                            i, j, sumqo,ikmax,iomax,ioi,iof,ispinmax

      real(dp)           :: qspiral(3), ef, temp,nelect, entrp,qomax,qtol
#ifdef MPI
      integer            :: MPIerror
#endif
!     Dynamic arrays
      integer, dimension(:),     allocatable, save :: muo
      integer, dimension(:,:),   allocatable, save :: nocck
      logical, dimension(:,:,:), allocatable, save :: occup
#ifdef DEBUG
      call write_debug('    PRE initwf')
#endif
!     First call initialisation
      if (frstme) then
#ifdef MPI
        if(ParallelOverk) then
          call die ('initwf: TDDFT is not parallelized over k-points.')
        end if
#endif
!       Read spin-spiral wavevector (if defined)
        call readsp( qspiral, spiral )
        if (spiral.and.Node.eq.0) then
          if (gamma_scf) write(6,*) &
            'diagon: WARNING: spin-spiral requires k sampling'
          if (nspin.ne.4) write(6,*) &
            'diagon: WARNING: spin-spiral requires nspin=4'
        end if
        frstme = .false.
      end if
!     Get Node number and calculate local orbital range
#ifdef MPI
      call GetNodeOrbs(no_u,Node,Nodes,nuo)
#else
      nuo = no_u
#endif
!     Start time counter ................................................
      call timer( 'initwf', 1 )
!     Check internal dimensions ..........................................
      if (nspin.le.2 .and. gamma_scf) then
        nhs  = no_u * nuo
        npsi = no_u * no_l * nspin
      else if (nspin.le.2 .and. .not.gamma_scf) then
        nhs  = 2 * no_u * nuo
        npsi = 2 * no_u * nuo
      else if (nspin.eq.4) then 
          call die ('initwf: TDDFT not yet      &
                      &implemented for non-collinear spin' ) 
      else
        call die('diagon: ERROR: incorrect value of nspin')
      end if
!     Allocate local arrays
      call re_alloc(Haux,1,nhs,name='Haux',routine='initwf')
      call re_alloc(Saux,1,nhs,name='Saux',routine='initwf')
      call re_alloc(psi,1,npsi,name='psi',routine='initwf')
      allocate(muo(nuo),stat=mem_stat)
      call memory('A','I',nuo,'initwf',stat=mem_stat)
      allocate(nocck(kpoint_scf%N,nspin),stat=mem_stat)
      call memory('A','I',kpoint_scf%N*nspin,'initwf',stat=mem_stat)
      allocate(occup(no_u,nspin,kpoint_scf%N),stat=mem_stat)
      call memory('A','L',nuo*kpoint_scf%N*nspin,'initwf',stat=mem_stat)
!     Check indxuo .......................................................
      do iuo = 1,nuo
        muo(iuo) = 0
      end do
      do io = 1,no_s
        iuo = indxuo(io)
        if (iuo.le.0 .or. iuo.gt.no_u) then
          if (Node.eq.0) then
            write(6,*) 'initwf: ERROR: invalid index: io, indxuo =',io, indxuo(io)
          end if
            call die('initwf: ERROR: invalid indxuo')
        end if
        call GlobalToLocalOrb(indxuo(io),Node,Nodes,iuo)
        if (iuo.gt.0) then
          muo(iuo) = muo(iuo) + 1
        end if
      end do
      do iuo = 1,nuo
        if (muo(iuo) .ne. muo(1)) then
          if (Node.eq.0) then
            write(6,'(/,2a,3i6)') 'initwf: ERROR: inconsistent indxuo', &
             '. iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
          end if
            call die ( 'initwf: ERROR: inconsistent indxuo.')
        end if
      end do
! ............................................................................!
!     Determine the number of occupied states this is not necessarily trivial !
!     if we have a metal. Only the occupied KS are saved and subsequently     !
!     evolved by integrating TDKS equations.                                  !
! ............................................................................!
      temp=1.0d-6
      call fermid( nspin, nspin, kpoint_scf%N, kpoint_scf%w, no_u, no_u, eo, &
                   temp, qtot, qo, ef, entrp )
      nocc(1) = 0
      nocc(2) = 0
      nelect=0.0d0
      degen= .false.
      !
      !
      do ik=1,kpoint_scf%N
        do ispin=1,nspin
          nocck(ik,ispin)=0
          do io=1,no_u
            occup(io,ispin,ik)=.false.
            if(dabs(qo(io,ispin,ik)-2.0d0*kpoint_scf%w(ik)/nspin).le.    &
               1.0d-2*dabs(2.0d0**kpoint_scf%w(ik)/nspin))  then
              nocc(ispin)=nocc(ispin)+1
              nocck(ik,ispin)=nocck(ik,ispin)+1
!             Accounting the number of electrons corresponding the states being marked
!             as occupied.
              nelect=nelect+dabs(2.0d0*kpoint_scf%w(ik)/nspin)
              occup(io,ispin,ik)=.true.
            else
              if ( dabs( qo(io,ispin,ik)) .gt.1.0d-2*dabs(2.0d0*kpoint_scf%w(ik)/nspin)) then
                IF (Node .eq. 0) THEN
                  IF(.not. degen) write(6,fmt="(/,a,tr3,a,tr3,a,tr3,a)") "initwf:","ik", &
                         "occupancy","maximum occupancy"
                  write(6,"(tr2,I10,tr3,f8.6,tr4,f8.6)") ik, qo(io,ispin,ik), &
                         2.0d0*kpoint_scf%w(ik)/nspin
                END IF
                degen = .true.  
              end if
            end if
          end do
        end do
      end do
!-------------------------------------------------------------------------------!
! Systems with odd number of electrons in spin-unpolarized calctions            !
! may encounter situations with partial occupations. In such case the           !
! program will stop.                                                            !
!-------------------------------------------------------------------------------!
!     Stop if the system has degenracy.
   if(Node .eq. 0) then
     write(6,'(/,(a,F18.6))') "initwf: No. of electrons corresponding occupied states =  ",nelect, &
                              "initwf: (Total charge - charge in selected states)     =  ",qtot-nelect
   end if
   !
   if (degen) then
     IF(Node .eq. 0) THEN
       Write(6,'(/,a,/,a)') "initwf: ERROR: System has degeneracy.", &
            "Change spin polarization, k-point sampling or shift to avoid it"
     END IF
     call die ('initwf: TDDFT doesnot allow degeneracy')
   end if
!..............
#ifdef MPI
      call ms_scalapack_setup(mpi_comm_world,1,'c',BlockSize)
      m_storage='pzdbc'
#else
      m_storage='szden'
#endif
      allocate(wavef_ms(1:kpoint_scf%N,1:nspin)) ! allocate (nkpnt*npsin) matrices inside wavef_ms
      do i=1,kpoint_scf%N !for every value of nkpnt and nspin, allocate a matrix of size (no_u x nocck(i,j))
        do j=1,nspin
          call m_allocate(wavef_ms(i,j),no_u,nocck(i,j),m_storage)
        end do
      end do
!     Call apropriate routine .............................................
      if (nspin.le.2 .and. gamma_scf) then
        call diaggiwf( nspin, nuo, no_l, maxnh, no_u,                     &
                    Haux, Saux, psi, no_u, occup)
      else if (nspin.le.2 .and. .not.gamma_scf) then
          call diagkiwf( nspin, nuo, no_s, nspin, no_l, maxnh,                 &
                         no_u, indxuo, kpoint_scf%N, kpoint_scf%k, Haux, Saux, &
                         psi, no_u, occup)
      else 
         call die('initwf: ERROR: non-collinear spin options for TDDFT not yet implemented')
      end if
!     Write/save wavefunction in .TDWF file to use for TDDFT calculation.
      IF (Node .eq. 0) WRITE(6,'(a)') 'initwf: Saving wavefunctions & 
                  &in <systemlabel>.TDWF file.'
      call  iowavef('write',wavef_ms,no_u,kpoint_scf%N,nspin)
!     Free local arrays
      call memory('D','I',size(muo),'initwf',stat=mem_stat)
      deallocate(muo,stat=mem_stat)
      call memory('D','I',size(nocck),'initwf',stat=mem_stat)
      deallocate(nocck,stat=mem_stat)
      call memory('D','L',size(occup),'initwf',stat=mem_stat)
      deallocate(occup,stat=mem_stat)

      call de_alloc( Haux, 'Haux', 'initwf')
      call de_alloc( Saux, 'Saux', 'initwf')
      call de_alloc( psi,  'phi',  'initwf')

#ifdef DEBUG
      call write_debug('    POS initwf')
#endif
!     Stop time counter ...................................................
      call timer( 'initwf', 2 )
      !
  end subroutine initwf
  ! Gamma point: solve KS by diagonalisation and store the occupied wavefunctions in wavef
  subroutine diaggiwf(nspin,nuo,maxuo,maxnh, maxo,Haux,Saux,psi,           &
                      nuotot,occup)
    use parallel, only : BlockSize
#ifdef MPI
    use m_diag, only: diag_descinit
#endif
      !
      implicit none
      !
      integer, intent(in)         :: maxnh, maxuo, maxo, nuo, nspin, nuotot
      real(dp), intent(inout)  :: Haux(nuotot,nuo), Saux(nuotot,nuo), psi(nuotot,maxuo,nspin)
      logical, intent(inout)      :: occup(nuotot,nspin,1)
      ! Internal variables
      integer                     :: ie, io, ispin, j, jo, ind, ierror, ioc, indwf
      real(dp)                    :: element
#ifdef MPI
      integer                     :: desch(9)
#endif
      !
#ifdef MPI
      call diag_descinit(nuotot,nuotot,BlockSize,desch)
#endif
      !
      indwf=0
      do ispin=1,nspin
        do ie=1,10
          Saux(1:nuotot,1:nuo) = 0.0d0
          Haux(1:nuotot,1:nuo) = 0.0d0
          do io=1,nuo
            do j=1,numh(io)
              ind=listhptr(io)+j
              jo=listh(ind)
              Saux(jo,io)=Saux(jo,io)+S(ind)
              Haux(jo,io)=Haux(jo,io)+H(ind,ispin)
            end do
          end do
          call rdiag(Haux,Saux,nuotot,nuo,nuotot,eo,psi(1,1,ispin),nuotot,1,ierror,BlockSize)
          if (ierror .eq. 0) then
            exit
          else if ((ierror .ne. -1) .or. (ie .eq. 10)) then
          call die('Terminating due to failed diagonalisation')
          end if
        end do     ! ie
!.............................       
        ioc=0
        do ie=1,nuotot
          if (occup(ie,ispin,1)) then 
            ioc=ioc+1
            indwf=indwf+1
            do j=1,nuotot
#ifdef MPI
              call pdelget('a',' ',element,psi(:,:,ispin),j,ie,desch)
#else
              element=psi(j,ie,ispin)
#endif
              call m_set_element(wavef_ms(1,ispin),j,ioc,cmplx(element,0.0,dp),complx_0,'lap')
            end do          ! j=1,nuotot
          end if            ! occup
        end do              ! ie=1,nuotot
      end do                ! ispin
      !
  end subroutine diaggiwf
  ! k points: solve KS by diagonalisation and store the occupied wavefunctions in wavef
  subroutine diagkiwf(nspin,nuo,no,maxspn,maxuo,maxnh, maxo, indxuo,nk,        &
                      kpoint, Haux,Saux,psi,nuotot,occup)
#ifdef MPI
      use m_diag, only: diag_descinit
#endif
      use parallel, only : BlockSize
      !
      implicit none
      !
      integer, intent(in)      :: maxnh, maxuo, maxo, no, nspin, nuo,          &
                                  nuotot, nk, maxspn, indxuo(no)
      real(dp), intent(in)     :: kpoint(3,nk)
      real(dp), intent(inout)  :: Haux(2,nuotot,nuo), Saux(2,nuotot,nuo), psi(2,nuotot,nuo)
      logical, intent(inout)   :: occup(nuotot,nspin,nk)
      ! Internal variables
      integer                  :: ispin, ie, ierror, ik, ind, iuo, j, jo, juo, indwf, ioc
      real(dp)                 :: ckxij, kxij, skxij, varaux(2)
      complex(dp)              :: element
#ifdef MPI
      integer                  :: desch(9)
#endif
      !
#ifdef MPI
      call diag_descinit(nuotot,nuotot,BlockSize,desch)
#endif
      !
      indwf=0
      do ik=1,nk
        do ispin=1,nspin
          Saux(1:2,1:nuotot,1:nuo)=0.0d0
          Haux(1:2,1:nuotot,1:nuo)=0.0d0
          do iuo=1,nuo
            do j=1,numh(iuo)
              ind=listhptr(iuo)+j
              jo=listh(ind)
              juo=indxuo(jo)
              kxij=kpoint(1,ik)*xijo(1,ind)+                      &
              kpoint(2,ik)*xijo(2,ind)+                           &
              kpoint(3,ik)*xijo(3,ind)
              ckxij=cos(kxij)
              skxij=sin(kxij)
!              Note: sign of complex part changed to match change in order of iuo/juo
              Saux(1,juo,iuo)=Saux(1,juo,iuo)+S(ind)*ckxij
              Saux(2,juo,iuo)=Saux(2,juo,iuo)-S(ind)*skxij
              Haux(1,juo,iuo)=Haux(1,juo,iuo)+H(ind,ispin)*ckxij
              Haux(2,juo,iuo)=Haux(2,juo,iuo)-H(ind,ispin)*skxij
            end do
          end do
          !
          call cdiag(Haux,Saux,nuotot,nuo,nuotot,eo(1,ispin,ik),psi,nuotot,1,ierror,BlockSize)
          if (ierror .ne. 0) then
          call die('Terminating due to failed diagonalisation')
          end if
          !.....................
          ioc=0
          do ie=1,nuotot
            if (occup(ie,ispin,ik)) then 
              ioc=ioc+1
              indwf=indwf+1
              do j=1,nuotot
#ifdef MPI
                call pzelget('a',' ',element,psi(:,:,:),j,ie,desch)
#else
                element=cmplx(psi(1,j,ie),psi(2,j,ie),dp)
#endif
                call m_set_element(wavef_ms(ik,ispin),j,ioc,element,complx_0,'lap')
              end do ! j
            end if
          end do ! ie
        end do  ! do ispin
      end do    ! do ikmax
  end subroutine diagkiwf


end module m_initwf
