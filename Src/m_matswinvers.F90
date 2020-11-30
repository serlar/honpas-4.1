 module matswinversion
 use precision
 use MatrixSwitch
 implicit none
 private
 public 						:: getinverse

 contains
  
 subroutine getinverse(C,label)
 implicit none
 character(3),intent(in), optional                      :: label
 type(matrix), intent(inout) 				:: C 

! We make a disctintion between symmetric and hermitician matrices

 if (C%is_real) then
     call m_dinver(C)
 else
     call m_zinver(C)
 end if

 end subroutine getinverse

subroutine m_dinver(C,label)

  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot, lwork, i, j, info
  integer, allocatable :: ipiv(:)
  real(dp), allocatable :: work(:)
#ifdef MPI
  integer :: liwork
  integer, allocatable :: iwork(:)
#else
  !**** EXTERNAL ********************************!
  integer, external :: ilaenv
#endif

  !**********************************************!

  if (.not. C%is_square) call die('m_dinver: matrix C is not square')

#ifdef MPI
      allocate(ipiv(C%iaux1(3)+C%iaux1(5)))
      call pdgetrf(C%dim1,C%dim2,C%dval,1,1,C%iaux1,ipiv,info)
      if (info/=0) call die('m_dinver: error in pdgetrf')
      allocate(work(1))
      allocate(iwork(1))
      call pdgetri(C%dim1,C%dval,1,1,C%iaux1,ipiv,work,-1,iwork,-1,info)
      if (info/=0) call die('m_dinver: error in pdgetri')
      liwork=iwork(1)
      deallocate(iwork)
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork))
      allocate(iwork(liwork))
      call pdgetri(C%dim1,C%dval,1,1,C%iaux1,ipiv,work,lwork,iwork,liwork,info)
      if (info/=0) call die('m_dinver: error in pdgetri')
      deallocate(iwork)
      deallocate(work)
      deallocate(ipiv)
#else
      allocate(ipiv(C%dim1))
      lwork=C%dim1*ilaenv(1,'dsytrf','u',C%dim1,-1,-1,-1)
      allocate(work(lwork))
      call dsytrf('u',C%dim1,C%dval,C%dim1,ipiv,work,lwork,info)
      if (info/=0) call die('m_dinver: error in dsytrf')
      deallocate(work)
      allocate(work(C%dim1))
      call dsytri('u',C%dim1,C%dval,C%dim1,ipiv,work,info)
      if (info/=0) call die('m_dinver: error in dsytri')
      deallocate(work)
      deallocate(ipiv)
      do i=1,C%dim1-1
        do j=i+1,C%dim1
          C%dval(j,i)=C%dval(i,j)
        end do
      end do

#endif

end subroutine m_dinver

subroutine m_zinver(C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot, lwork, i, j, info
  complex(dp), allocatable :: work(:)
  integer, allocatable :: ipiv(:)
#ifdef MPI
  integer :: liwork
  integer, allocatable :: iwork(:)
#else
  !**** EXTERNAL ********************************!
  integer, external :: ilaenv
#endif

  !**********************************************!

  if (.not. C%is_square) call die('m_zinver: matrix C is not square')

#ifdef MPI



      allocate(ipiv(C%iaux1(3)+C%iaux1(5)))
      call pzgetrf(C%dim1,C%dim2,C%zval,1,1,C%iaux1,ipiv,info)
      if (info/=0) call die('m_zinver: error in pzgetrf')
      allocate(work(1))
      allocate(iwork(1))

      call pzgetri(C%dim1,C%zval,1,1,C%iaux1,ipiv,work,-1,iwork,-1,info)
      if (info/=0) call die('m_zinver: error in pzgetri')
      liwork=iwork(1)
      deallocate(iwork)
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork))
      allocate(iwork(liwork))
 




      call pzgetri(C%dim1,C%zval,1,1,C%iaux1,ipiv,work,lwork,iwork,liwork,info)
 

      if (info/=0) call die('m_zinver: error in pzgetri')
      deallocate(iwork)
      deallocate(work)
      deallocate(ipiv)


#else

      allocate(ipiv(C%dim1))
     call zgetrf(C%dim1,C%dim2,C%zval,C%dim1,ipiv,info)
      if (info/=0) call die('m_zinver: error in zhetrf')
      allocate(work(1))
      call zgetri(C%dim1,C%zval,C%dim1,ipiv,work,-1,info)
      lwork=work(1) 
      deallocate(work)
      allocate(work(lwork))
      call zgetri(C%dim1,C%zval,C%dim1,ipiv,work,lwork,info)
      if (info/=0) call die('m_zinver: error in zhetri')
      deallocate(work)
      deallocate(ipiv)
#endif
 end subroutine m_zinver
 
 subroutine die(message)
  implicit none

  !**** INPUT ***********************************!

  character(*), intent(in), optional :: message

  !**********************************************!

  print*, 'getinverse:', message
  stop

end subroutine die
 
 end module matswinversion
