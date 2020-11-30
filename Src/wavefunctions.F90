MODULE  wavefunctions
  USE  precision
  USE  MatrixSwitch
  !
  implicit none 
  !    
  PRIVATE
  PUBLIC       :: iowavef, compute_tddm, compute_tdEdm
  type(matrix), allocatable, save, public :: wavef_ms(:,:)
  TYPE(matrix), ALLOCATABLE, SAVE, public :: Hsave(:,:)
  ! 
  complex(kind=dp), parameter, public :: complx_1 = cmplx(1.0,0.0,dp)
  complex(kind=dp), parameter, public :: complx_0 = cmplx(0.0,0.0,dp)

CONTAINS
  !
  subroutine iowavef(task,wavef_rw,nuotot,nk,nspin)
 
    ! To read and write the TDKS orbitals
    ! 
    ! Written by Daniel Sanchez-Portal
    ! Re-written by Rafi Ullah, CIC nanoGUNE, October 16, 2015 to work
    ! with parallel TDDFT flavor of siesta using MatrixSwtich.
    !
    ! Although it is prepared to work with parallel TDDFT-siesta the
    ! reading/writing process itself is not parallel.
    use precision
    use parallel
    use fdf
    use MatrixSwitch
#ifdef MPI
    use mpi_siesta
#endif
    !
    implicit none
    !
    integer, intent(in)             :: nuotot, nk, nspin
    character*(*), intent(in)       :: task
    type(matrix), intent(inout)     :: wavef_rw(nk,nspin)
    ! Internal variables and arrays
    character         :: fname*33, sname*30, m_storage*5
    logical           ::  exist1, frstme
    integer           :: unit1, ie,nuototread,nkread,nspinread,dim2read
    integer           :: nwf, ik, ispin, mxnwf, io, ix,i,j
    external          :: chkdim, io_assign, io_close, timer,memory
    complex(kind=dp)  :: varaux
    save              :: frstme, fname
    data                 frstme /.true./
    !
#ifdef MPI
    INTEGER :: MPIerror
#endif
    !
#ifdef MPI
    call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
    call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
    Node = 0
    Nodes = 1
#endif
    !
#ifdef MPI
      m_storage='pzdbc'
#else
      m_storage='szden'
#endif
    ! Find file name
    if (frstme) then
      if (Node.eq.0) then
        sname = fdf_string('SystemLabel','siesta')
      endif
#ifdef MPI
      call MPI_Bcast(sname,30,MPI_character,0,MPI_Comm_World,MPIerror)
#endif
      fname = trim(sname)//'.TDWF'
      frstme = .false.
    endif
    !
    if (task.eq.'read' .or. task.eq.'READ') then      ! task read
      if (Node.eq.0) then
        inquire (file=fname, exist=exist1)
        if(exist1) then
          write(6,'(/,a)')                                            &
          'iowavef: reading initial KS wavefunctions from file'
         else
          write(6,'(/,a)')                                            &
          'iowavef: file containing the KS orbitals not found'
          write(6,'(a)')                                              &
          'iowavef: simulation can not be started or restarted'
          stop 
        endif
        call io_assign(unit1)
        open( unit1, file=fname, form='unformatted',                 &
        status='old' )
      endif
      if (Node.eq.0) then
        read(unit1) nuototread, nkread, nspinread
        if(nkread.ne.nk) stop 'iowavef: Nunber of K-points inconsistent '
        if(nuototread.ne.nuotot) stop 'iowavef: Number of KS orbitals inconsistent ' 
        if(nspinread.ne.nspin) stop 'iowavef: Spin inconsistent '
      endif
      !   
      do ispin=1,nspin                                 
        do ik=1,nk
          if (Node.eq.0) read(unit1) dim2read
#ifdef MPI
          call MPI_Bcast(dim2read,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
          call m_allocate(wavef_rw(ik,ispin),nuotot,dim2read,m_storage)
        enddo
      enddo
      !
      do ispin=1,nspin
        do ik=1,nk
          do i=1,nuotot
            do j=1,wavef_rw(ik,ispin)%dim2
              if (Node==0) read(unit1) varaux
#ifdef MPI 
              call MPI_Bcast(varaux,1,MPI_double_complex,0,MPI_Comm_World,MPIerror)
#endif
              call m_set_element(wavef_rw(ik,ispin),i,j,varaux,complx_0,'lap')
            end do
          end do
        enddo
      enddo
      if (Node.eq.0)   call io_close(unit1)
      !......................................................................   
    elseif(task.eq.'write'.or.task.eq.'WRITE') then     ! task write
      if (Node.eq.0) then
        call io_assign(unit1)
        open( unit1, file=fname, form='unformatted',status='replace' )
        rewind(unit1)
        write(unit1) nuotot,nk,nspin
        !
        do ispin=1,nspin
          do ik=1,nk
            write(unit1) (wavef_rw(ik,ispin)%dim2)
          enddo
        enddo
        !
      end if
      !
      do ispin=1,nspin
        do ik=1,nk
          do i=1,nuotot
            do j=1,wavef_rw(ik,ispin)%dim2
              call m_get_element(wavef_rw(ik,ispin),i,j,varaux,'lap')
              if (Node==0) write(unit1) varaux
            enddo
          enddo
        enddo 
      enddo 
      !
      if (Node .eq. 0) then
        call io_close(unit1)
      endif
    endif            ! end task 
    !
  end subroutine iowavef
  !
  subroutine  compute_tddm (Dnew)
  
    ! To calculate the time-dependent density matrix from time-dependent
    ! wavefunctions. 
    ! Written by Rafi Ullah January 2017
    !***********************************
    use sparse_matrices,      only: numh, maxnh, listh, listhptr,     &
                                    xijo   
    use kpoint_scf_m,         only: kpoint_scf, gamma_scf
    use atomlist,             only: no_l, no_u, indxuo
    use m_spin,               only: nspin
    integer                      :: ispin, nuo, nuotot
    integer                      :: io,jo, juo, j, ind, ik
    real(dp)                     :: Dnew (maxnh, nspin), wk, kxij
    real(dp)                     :: ckxij, skxij
    complex(dp)                  :: varaux
    type(matrix)                 :: Daux    
    character                    :: m_storage*5, m_operation*3
#ifdef MPI
    m_storage   ='pzdbc'
    m_operation = 'lap'
#else
    m_storage='szden'
    m_operation = 'lap'
#endif
    Dnew(1:maxnh,1:nspin) = 0.0_dp
    call m_allocate(Daux,no_u,no_u,m_storage)
    Do ik = 1, kpoint_scf%N
      Do ispin =1, nspin
      wk = 2.00_dp*kpoint_scf%w(ik)/dble(nspin)
      ! Calculating density matrix using MatrixSwitch.
      call mm_multiply(wavef_ms(ik,ispin),'n',wavef_ms(ik,ispin),'c',Daux,           & 
                    cmplx(wk,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),m_operation)
      ! Passing DM from dense to sparse form.  
      DO io=1,no_l
        DO j=1,numh(io)
          ind=listhptr(io) + j
          juo = listh(ind)
          jo  = indxuo(juo)
          IF (.NOT. gamma_scf) THEN
            kxij = kpoint_scf%k(1,ik) * xijo(1,ind) + &
            kpoint_scf%k(2,ik) * xijo(2,ind) +        &
            kpoint_scf%k(3,ik) * xijo(3,ind) 
            ckxij =  cos(kxij)
            skxij = -sin(kxij)
          ELSE
            ckxij = 1.0_dp
            skxij = 0.0_dp
          END IF
          varaux = real(Daux%zval(jo,io))*ckxij +     &
                   aimag(Daux%zval(jo,io))*skxij
          Dnew(ind,ispin)  = Dnew(ind,ispin) + varaux
        END DO
      END DO
    END DO
  END DO
    call m_deallocate (Daux)
  end subroutine compute_tddm
!----------------------------------------------------------------------------------!


!----------------------------------------------------------------------------------!
      subroutine compute_tdEdm(Enew)
!**********************************************************************************!
! Re-written by Rafi Ullah in June 2017.                                           !
! In this subroutine we calculate the energy density matrix in parallel using      !
! MatrixSwitch.                                                                    !
!**********************************************************************************!
      use parallel
      use precision
      use parallelsubs,     only: LocalToGlobalOrb
      use sparse_matrices,  only: numh, maxnh, listh, listhptr, xijo
      use sparse_matrices,  only: H, S
      use kpoint_scf_m,     only: kpoint_scf, gamma_scf
      use atomlist,         only: no_l, no_u, indxuo
      use m_spin,           only: nspin
      use MatrixSwitch
      use matswinversion,   only: getinverse
      !
      implicit none
      !
      real(dp)               :: Enew(maxnh, nspin)
      real(dp)               :: wk, kxij, ckxij, skxij
      integer                :: ik, ispin, i, j, io, jo, ind, juo, nocc
      logical, save          :: frstime = .true.
      character              :: m_storage*5, m_operation*3
      complex(dp)            :: cvar1, cvar2
      !
      type(matrix)           :: Sauxms, Hauxms,psi,Eaux
      type(matrix)           :: S_1
      !
#ifdef MPI
      m_storage='pzdbc'
      m_operation='lap'
#else
      m_storage='szden'
      m_operation='lap'
#endif
      Enew(1:maxnh,1:nspin) = 0.0_dp
      call m_allocate(S_1,no_u,no_u,m_storage)
      call m_allocate(Eaux,no_u,no_u,m_storage)
      Do ik=1,kpoint_scf%N
        Do ispin=1,nspin
          wk = 2.00_dp*kpoint_scf%w(ik)/dble(nspin)
          nocc = wavef_ms(ik,ispin)%dim2
          call m_allocate(psi,nocc,no_u,m_storage)
          call m_allocate (Hauxms,no_u, no_u, m_storage)
          call m_allocate (Sauxms,no_u, no_u, m_storage)
          Do i = 1, no_l
            call LocalToGlobalOrb(i, Node, Nodes, io)
            Do j = 1, numh(i)
              ind = listhptr(i) + j
              juo = listh(ind)
              jo  = indxuo(juo)
              if( .not. gamma_scf ) then
                kxij = kpoint_scf%k(1,ik)*xijo(1,ind) + kpoint_scf%k(2,ik)*xijo(2,ind) +         &
                       kpoint_scf%k(3,ik)*xijo(3,ind)
                ckxij =  cos(kxij)
                skxij = -sin(kxij)
              else
                ckxij = 1.0_dp
                skxij = 0.0_dp
              end if
              cvar1 = cmplx(H(ind,ispin)*ckxij,H(ind,ispin)*skxij,dp)
              cvar2 = cmplx(S(ind)*ckxij,S(ind)*skxij,dp)
              call m_set_element(Hauxms,jo,io,cvar1,complx_1,m_operation)
              call m_set_element(Sauxms,jo,io,cvar2,complx_1,m_operation)
            end do
          end do
          !
          if (ispin .eq. 1) then
            ! copy Sauxms to S_1
            call m_add (Sauxms,'n',S_1,complx_1,complx_0,m_operation)
            ! Invert the overlap matrix.
            call getinverse(S_1,m_operation)
          end if
          call mm_multiply(Hauxms,'n',S_1,'n',Eaux,complx_1,complx_0,m_operation)
          call mm_multiply(wavef_ms(ik,ispin),'c',Eaux,'n',psi,complx_1,            &
                           complx_0,m_operation)
          call mm_multiply(wavef_ms(ik,ispin),'n',psi,'n',Eaux,cmplx(wk,0.0,dp),    &
                           complx_0,m_operation)
          ! Passing Energy-DM from dense to sparse form.  
          DO io=1,no_l
            DO j=1,numh(io)
              ind=listhptr(io) + j
              juo = listh(ind)
              jo  = indxuo(juo)
              IF (.NOT. gamma_scf) THEN
                kxij = kpoint_scf%k(1,ik) * xijo(1,ind) + &
                kpoint_scf%k(2,ik) * xijo(2,ind) +        &
                kpoint_scf%k(3,ik) * xijo(3,ind) 
                ckxij =  cos(kxij)
                skxij = -sin(kxij)
              ELSE
                ckxij = 1.0_dp
                skxij = 0.0_dp
              END IF
              cvar1 =  real(Eaux%zval(jo,io))*ckxij +     &
                       aimag(Eaux%zval(jo,io))*skxij
              Enew(ind,ispin)  = Enew(ind,ispin) + cvar1
            END DO
          END DO
          call m_deallocate(Sauxms)
          call m_deallocate(Hauxms)
          call m_deallocate(psi)
        end do
      end do
      call m_deallocate(S_1)
      call m_deallocate(Eaux)
      end subroutine compute_tdEdm 
  !
end module wavefunctions       
