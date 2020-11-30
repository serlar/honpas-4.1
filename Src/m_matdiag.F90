 module matdiagon
 use precision
 use MatrixSwitch

 implicit none
 private
 public 						:: geteigen
  
 contains
 
 subroutine geteigen(C,eig,label)

! *********************************************************************************************
! Subroutine to calculate the eigenvalues and eigenvectors,for a given  MatrixSwitch derived 
! type matrix. It was created to replace the Scalapack calls in D. Sanchez-Portal's changebasis
! Written by Adiran Garaizar, June 2015.
! **************************** INPUT **********************************************************
! character(3) label               : MatrixSwitch label identifying the storage format   
! **************************** INPUT and OUTPUT ***********************************************
! type(matrix) C                   : MatrixSwitch derived type matrix, eigenvectors as output 
! *************************** OUTPUT **********************************************************
! real(dp) eig                     : Vector conaining the eigenvalues
! *************************** PARALLEL ********************************************************
! When running in parallel this routine uses Scalapack to perform a
! parallel matrix diagonalisation.
! *********************************************************************************************

 implicit none

 character(3), intent(in), optional 			:: label
 type(matrix), intent(inout) 				:: C 
 real(kind=dp),dimension(:),intent(out)		:: eig 

! We make a disctintion between symmetric and hermitician matrices

 if (C%is_real) then
     call m_deigen(C,eig,label)
 else
     call m_zeigen(C,eig,label)
 end if
 
 end subroutine geteigen 
 
 subroutine m_deigen(C,eig,label)

 implicit none

 character(3), intent(in), optional			:: label 
 type(matrix), intent(inout) 				:: C 
 type(matrix)  						:: eigvec
 real(kind=dp),dimension(C%dim1),intent(out)		:: eig 
 integer						:: no,lwork,info
 real(kind=dp),dimension(:),allocatable                 :: work

 no=C%dim1
 
#ifdef MPI
! First call to get work matrices

 	call m_allocate(eigvec,no,no,'pddbc')
        allocate(work(1))
        call pdsyev ('V', 'U', no, C%dval, 1, 1, C%iaux1, eig,eigvec%dval,&
 	 1, 1, eigvec%iaux1, work, -1, info )
        lwork=work(1)
        deallocate(work)

! Second call to get eigenvalues and eigenvectors

        allocate(work(lwork))
 	call pdsyev ('V', 'U', no, C%dval, 1, 1, C%iaux1, eig,eigvec%dval,&
 	 1, 1, eigvec%iaux1, work, lwork, info )
 	deallocate(work)
	call m_add(eigvec,'n',C,1.0_dp,0.0_dp,'lap')
        call m_deallocate(eigvec)
 	if(info.ne.0) stop 'Error calculatesqrtS'
#else
! First call to get work matrices

        allocate(work(1))
        call dsyev  ('V', 'U', no, C%dval,no, eig,work, -1, info)
        lwork=work(1)
        deallocate(work)

! Second call to get eigenvalues and eigenvectors

        allocate(work(lwork))
 	call dsyev  ('V', 'U', no, C%dval,no, eig,work, lwork, info)
        deallocate(work)
       	if(info.ne.0) stop 'Error calculatesqrtS'
#endif

 end subroutine m_deigen
 
 subroutine m_zeigen (C,eig,label)

 implicit none

 character(3), intent(in), optional			:: label 
 type(matrix), intent(inout) 				:: C 
 type(matrix) 						:: eigvec 
 real(kind=dp),dimension(:),intent(out)		:: eig 
 integer						:: no,lwork,lrwork,info
 complex(kind=dp),dimension(:),allocatable              :: work,zrwork
 real(kind=dp),dimension(:),allocatable                 :: rwork

 no=C%dim1
  
#ifdef MPI
! First call to get work matrices
  	call m_allocate(eigvec,no,no,'pzdbc')
        allocate(work(1))
        allocate(zrwork(1))
    call pzheev('V', 'U',no, C%zval, 1, 1, C%iaux1,eig,eigvec%zval,&
 	 1, 1, eigvec%iaux1, work, -1, zrwork, -1, info )
        lwork=work(1)
        lrwork=zrwork(1)
        deallocate(zrwork)
        deallocate(work)

! Second call to get eigenvalues and eigenvectors

        allocate(work(lwork))
        allocate(zrwork(lrwork))
 	call pzheev('V', 'U', no, C%zval, 1, 1, C%iaux1, eig,eigvec%zval,&
 	 1, 1, eigvec%iaux1, work, lwork, zrwork, lrwork, info )
        deallocate(zrwork)
        deallocate(work)
	call m_add(eigvec,'n',C,1.0_dp,0.0_dp,'lap')
        call m_deallocate(eigvec)
       	if(info.ne.0) stop 'Error calculatesqrtS'
#else
! First call to get work matrices

        allocate(rwork((max(1,3*no-2))))
        allocate(work(1))
 	call zheev('V','U',no,C%zval,no,eig,work,-1,rwork,info)
        lwork=work(1)
        deallocate(work)
! Second call to get eigenvalues and eigenvectors

        allocate(work(lwork))
 	call zheev('V','U',no,C%zval,no,eig,work,lwork,rwork,info)
	deallocate(work)
        deallocate(rwork)
       	if(info.ne.0) stop 'Error calculatesqrtS'
#endif 

 end subroutine m_zeigen
 end module matdiagon
