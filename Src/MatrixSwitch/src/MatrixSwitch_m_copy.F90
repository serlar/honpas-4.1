#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==============================================================================!
!> @brief Implementations of \a m_copy.
!==============================================================================!
module MatrixSwitch_m_copy
  use MatrixSwitch_ops

  implicit none

  !**** INTERFACES ********************************!

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Copy external matrix (sparse coordinate list from dense block
  !!        cyclic, parallel distribution).
  !!
  !! Initializes a TYPE(MATRIX) variable with \c p?coo format by copying the
  !! information and element values from pre-existing matrix data with \c p?dbc
  !! format. This is a dense to sparse conversion.
  !!
  !! @param[inout] m_name    The matrix to be allocated.
  !! @param[in]    A         The values of the local matrix elements for the
  !!                         matrix to be copied, stored as a two-dimensional
  !!                         array.
  !! @param[in]    desc      BLACS array descriptor for the matrix to be
  !!                         copied.
  !! @param[in]    threshold Tolerance for zeroing elements. Elements with an
  !!                         absolute value below this threshold will be
  !!                         omitted.
  !============================================================================!
  interface m_copy_external_pdbcpcoo
     module procedure m_copy_external_pddbcpdcoo
     module procedure m_copy_external_pzdbcpzcoo
  end interface m_copy_external_pdbcpcoo
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Copy external matrix (compressed sparse column from dense block
  !!        cyclic, parallel distribution).
  !!
  !! Initializes a TYPE(MATRIX) variable with \c p?csc format by copying the
  !! information and element values from pre-existing matrix data with \c p?dbc
  !! format. This is a dense to sparse conversion.
  !!
  !! @param[inout] m_name    The matrix to be allocated.
  !! @param[in]    A         The values of the local matrix elements for the
  !!                         matrix to be copied, stored as a two-dimensional
  !!                         array.
  !! @param[in]    desc      BLACS array descriptor for the matrix to be
  !!                         copied.
  !! @param[in]    threshold Tolerance for zeroing elements. Elements with an
  !!                         absolute value below this threshold will be
  !!                         omitted.
  !============================================================================!
  interface m_copy_external_pdbcpcsc
     module procedure m_copy_external_pddbcpdcsc
     module procedure m_copy_external_pzdbcpzcsc
  end interface m_copy_external_pdbcpcsc
#endif

  !************************************************!

contains

  !============================================================================!
  !> @brief Copy matrix (simple dense, serial distribution, reference
  !!        implementation).
  !============================================================================!
  subroutine m_copy_sdensdenref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_sdensdenref

  !============================================================================!
  !> @brief Copy matrix (simple dense, serial distribution, reference
  !!        implementation with thresholding).
  !============================================================================!
  subroutine m_copy_sdensdenref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          do i=1,m_name%dim1
             do j=1,m_name%dim2
                if (abs(A%dval(i,j))>abs_threshold) then
                   m_name%dval(i,j)=A%dval(i,j)-soft_threshold*A%dval(i,j)/abs(A%dval(i,j))
                 else
                   m_name%dval(i,j)=0.0_dp
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          do i=1,m_name%dim1
             do j=1,m_name%dim2
                if (abs(A%zval(i,j))>abs_threshold) then
                   m_name%zval(i,j)=A%zval(i,j)-soft_threshold*A%zval(i,j)/abs(A%zval(i,j))
                else
                   m_name%zval(i,j)=cmplx_0
                end if
             end do
          end do
       end if

  end subroutine m_copy_sdensdenref_thre

  !============================================================================!
  !> @brief Copy matrix (dense block cyclic, parallel distribution, reference
  !!        implementation).
  !============================================================================!
  subroutine m_copy_pdbcpdbcref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       allocate(m_name%iaux1(9))
       m_name%iaux1_is_allocated=.true.
       allocate(m_name%iaux2(2))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux1=A%iaux1
       m_name%iaux2=A%iaux2
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),m_name%iaux2(2)))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),m_name%iaux2(2)))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_pdbcpdbcref

  !============================================================================!
  !> @brief Copy matrix (dense block cyclic, parallel distribution, reference
  !!        implementation with thresholding).
  !============================================================================!
  subroutine m_copy_pdbcpdbcref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j

    !**********************************************!

       allocate(m_name%iaux1(9))
       m_name%iaux1_is_allocated=.true.
       allocate(m_name%iaux2(2))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux1=A%iaux1
       m_name%iaux2=A%iaux2
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),m_name%iaux2(2)))
          m_name%dval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=1,m_name%iaux2(2)
                if (abs(A%dval(i,j))>abs_threshold) then
                   m_name%dval(i,j)=A%dval(i,j)-soft_threshold*A%dval(i,j)/abs(A%dval(i,j))
                else
                   m_name%dval(i,j)=0.0_dp
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),m_name%iaux2(2)))
          m_name%zval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=1,m_name%iaux2(2)
                if (abs(A%zval(i,j))>abs_threshold) then
                   m_name%zval(i,j)=A%zval(i,j)-soft_threshold*A%zval(i,j)/abs(A%zval(i,j))
                else
                   m_name%zval(i,j)=cmplx_0
                end if
             end do
          end do
       end if

  end subroutine m_copy_pdbcpdbcref_thre

  !============================================================================!
  !> @brief Copy matrix (sparse coordinate list, serial distribution, reference
  !!        implementation).
  !============================================================================!
  subroutine m_copy_scooscooref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       m_name%iaux3=A%iaux3
       m_name%iaux4=A%iaux4
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_scooscooref

  !============================================================================!
  !> @brief Copy matrix (m_name = sparse coordinate list, A = simple dense,
  !!        serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_sdenscooref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=m_name%dim1*m_name%dim2
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          k=0
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                k=k+1
                m_name%iaux3(k)=j
                m_name%iaux4(k)=i
                m_name%dval(k,1)=A%dval(j,i)
             end do
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          k=0
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                k=k+1
                m_name%iaux3(k)=j
                m_name%iaux4(k)=i
                m_name%zval(k,1)=A%zval(j,i)
             end do
          end do
       end if

  end subroutine m_copy_sdenscooref

  !============================================================================!
  !> @brief Copy matrix (m_name = sparse coordinate list, A = simple dense,
  !!        serial distribution, reference implementation with thresholding).
  !============================================================================!
  subroutine m_copy_sdenscooref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=0
       if (m_name%is_real) then
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%dval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%iaux2(1)))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%dval(m_name%iaux2(1),1))
             m_name%dval_is_allocated=.true.
             k=0
             do i=1,m_name%dim2
                do j=1,m_name%dim1
                   if (abs(A%dval(j,i))>abs_threshold) then
                      k=k+1
                      m_name%iaux3(k)=j
                      m_name%iaux4(k)=i
                      m_name%dval(k,1)=A%dval(j,i)-soft_threshold*A%dval(j,i)/abs(A%dval(j,i))
                   end if
                end do
             end do
          end if
       else
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%zval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%iaux2(1)))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%zval(m_name%iaux2(1),1))
             m_name%zval_is_allocated=.true.
             k=0
             do i=1,m_name%dim2
                do j=1,m_name%dim1
                   if (abs(A%zval(j,i))>abs_threshold) then
                      k=k+1
                      m_name%iaux3(k)=j
                      m_name%iaux4(k)=i
                      m_name%zval(k,1)=A%zval(j,i)-soft_threshold*A%zval(j,i)/abs(A%zval(j,i))
                   end if
                end do
             end do
          end if
       end if

  end subroutine m_copy_sdenscooref_thre

  !============================================================================!
  !> @brief Copy matrix (compressed sparse column, serial distribution,
  !!        reference implementation).
  !============================================================================!
  subroutine m_copy_scscscscref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%dim2+1))
       m_name%iaux4_is_allocated=.true.
       m_name%iaux3=A%iaux3
       m_name%iaux4=A%iaux4
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_scscscscref

  !============================================================================!
  !> @brief Copy matrix (compressed sparse row, serial distribution, reference
  !!        implementation).
  !============================================================================!
  subroutine m_copy_scsrscsrref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%dim1+1))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       m_name%iaux3=A%iaux3
       m_name%iaux4=A%iaux4
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_scsrscsrref

  !============================================================================!
  !> @brief Copy matrix (m_name = compressed sparse column, A = simple dense,
  !!        serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_sdenscscref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=m_name%dim1*m_name%dim2
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%dim2+1))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          k=0
          m_name%iaux4(1)=0
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                k=k+1
                m_name%iaux3(k)=j
                m_name%dval(k,1)=A%dval(j,i)
             end do
             m_name%iaux4(i+1)=k
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          k=0
          m_name%iaux4(1)=0
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                k=k+1
                m_name%iaux3(k)=j
                m_name%zval(k,1)=A%zval(j,i)
             end do
             m_name%iaux4(i)=k+1
          end do
       end if

  end subroutine m_copy_sdenscscref

  !============================================================================!
  !> @brief Copy matrix (m_name = compressed sparse column, A = simple dense,
  !!        serial distribution, reference implementation with thresholding).
  !============================================================================!
  subroutine m_copy_sdenscscref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=0
       if (m_name%is_real) then
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%dval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%dim2+1))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%dval(m_name%iaux2(1),1))
             m_name%dval_is_allocated=.true.
             k=0
             m_name%iaux4(1)=0
             do i=1,m_name%dim2
                do j=1,m_name%dim1
                   if (abs(A%dval(j,i))>abs_threshold) then
                      k=k+1
                      m_name%iaux3(k)=j
                      m_name%dval(k,1)=A%dval(j,i)-soft_threshold*A%dval(j,i)/abs(A%dval(j,i))
                   end if
                end do
                m_name%iaux4(i+1)=k
             end do
          end if
       else
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%zval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%dim2+1))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%zval(m_name%iaux2(1),1))
             m_name%zval_is_allocated=.true.
             k=0
             m_name%iaux4(1)=0
             do i=1,m_name%dim2
                do j=1,m_name%dim1
                   if (abs(A%zval(j,i))>abs_threshold) then
                      k=k+1
                      m_name%iaux3(k)=j
                      m_name%zval(k,1)=A%zval(j,i)-soft_threshold*A%zval(j,i)/abs(A%zval(j,i))
                   end if
                end do
                m_name%iaux4(i)=k+1
             end do
          end if
       end if

  end subroutine m_copy_sdenscscref_thre

  !============================================================================!
  !> @brief Copy matrix (m_name = compressed sparse row, A = simple dense,
  !!        serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_sdenscsrref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=m_name%dim1*m_name%dim2
       allocate(m_name%iaux3(m_name%dim1+1))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          k=0
          m_name%iaux3(1)=0
          do i=1,m_name%dim1
             do j=1,m_name%dim2
                k=k+1
                m_name%iaux4(k)=j
                m_name%dval(k,1)=A%dval(i,j)
             end do
             m_name%iaux3(i+1)=k
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          k=0
          m_name%iaux3(1)=0
          do i=1,m_name%dim1
             do j=1,m_name%dim2
                k=k+1
                m_name%iaux4(k)=j
                m_name%zval(k,1)=A%zval(i,j)
             end do
             m_name%iaux3(i)=k+1
          end do
       end if

  end subroutine m_copy_sdenscsrref

  !============================================================================!
  !> @brief Copy matrix (m_name = compressed sparse row, A = simple dense,
  !!        serial distribution, reference implementation with thresholding).
  !============================================================================!
  subroutine m_copy_sdenscsrref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=0
       if (m_name%is_real) then
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%dval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%dim1+1))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%iaux2(1)))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%dval(m_name%iaux2(1),1))
             m_name%dval_is_allocated=.true.
             k=0
             m_name%iaux3(1)=0
             do i=1,m_name%dim1
                do j=1,m_name%dim2
                   if (abs(A%dval(i,j))>abs_threshold) then
                      k=k+1
                      m_name%iaux4(k)=j
                      m_name%dval(k,1)=A%dval(i,j)-soft_threshold*A%dval(i,j)/abs(A%dval(i,j))
                   end if
                end do
                m_name%iaux3(i+1)=k
             end do
          end if
       else
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%zval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%dim1+1))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%iaux2(1)))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%zval(m_name%iaux2(1),1))
             m_name%zval_is_allocated=.true.
             k=0
             m_name%iaux3(1)=0
             do i=1,m_name%dim1
                do j=1,m_name%dim2
                   if (abs(A%zval(i,j))>abs_threshold) then
                      k=k+1
                      m_name%iaux4(k)=j
                      m_name%zval(k,1)=A%zval(i,j)-soft_threshold*A%zval(i,j)/abs(A%zval(i,j))
                   end if
                end do
                m_name%iaux3(i)=k+1
             end do
          end if
       end if

  end subroutine m_copy_sdenscsrref_thre

  !============================================================================!
  !> @brief Copy matrix (m_name = simple dense, A = sparse coordinate list,
  !!        serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_scoosdenref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%iaux2(1)
             m_name%dval(A%iaux3(i),A%iaux4(i))=A%dval(i,1)
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%iaux2(1)
             m_name%zval(A%iaux3(i),A%iaux4(i))=A%zval(i,1)
          end do
       end if

  end subroutine m_copy_scoosdenref

  !============================================================================!
  !> @brief Copy matrix (m_name = simple dense, A = sparse coordinate list,
  !!        serial distribution, reference implementation with thresholding).
  !============================================================================!
  subroutine m_copy_scoosdenref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%iaux2(1)
             if (abs(A%dval(i,1))>abs_threshold) then
                m_name%dval(A%iaux3(i),A%iaux4(i))=A%dval(i,1)-soft_threshold*A%dval(i,1)/abs(A%dval(i,1))
             else
                m_name%dval(A%iaux3(i),A%iaux4(i))=0.0_dp
             end if
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%iaux2(1)
             if (abs(A%zval(i,1))>abs_threshold) then
                m_name%zval(A%iaux3(i),A%iaux4(i))=A%zval(i,1)-soft_threshold*A%zval(i,1)/abs(A%zval(i,1))
             else
                m_name%zval(A%iaux3(i),A%iaux4(i))=cmplx_0
             end if
          end do
       end if

  end subroutine m_copy_scoosdenref_thre

  !============================================================================!
  !> @brief Copy matrix (m_name = simple dense, A = compressed sparse column,
  !!        serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_scscsdenref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%dim2
             do j=1,A%iaux4(i+1)-A%iaux4(i)
                k=A%iaux4(i)+j
                m_name%dval(A%iaux3(k),i)=A%dval(k,1)
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%dim2
             do j=1,A%iaux4(i+1)-A%iaux4(i)
                k=A%iaux4(i)+j
                m_name%zval(A%iaux3(k),i)=A%zval(k,1)
             end do
          end do
       end if

  end subroutine m_copy_scscsdenref

  !============================================================================!
  !> @brief Copy matrix (m_name = simple dense, A = compressed sparse column,
  !!        serial distribution, reference implementation with thresholding).
  !============================================================================!
  subroutine m_copy_scscsdenref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%dim2
             do j=1,A%iaux4(i+1)-A%iaux4(i)
                k=A%iaux4(i)+j
                if (abs(A%dval(k,1))>abs_threshold) then
                   m_name%dval(A%iaux3(k),i)=A%dval(k,1)-soft_threshold*A%dval(k,1)/abs(A%dval(k,1))
                else
                   m_name%dval(A%iaux3(k),i)=0.0_dp
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%dim2
             do j=1,A%iaux4(i+1)-A%iaux4(i)
                k=A%iaux4(i)+j
                if (abs(A%zval(k,1))>abs_threshold) then
                   m_name%zval(A%iaux3(k),i)=A%zval(k,1)-soft_threshold*A%zval(k,1)/abs(A%zval(k,1))
                else
                   m_name%zval(A%iaux3(k),i)=0.0_dp
                end if
             end do
          end do
       end if

  end subroutine m_copy_scscsdenref_thre

  !============================================================================!
  !> @brief Copy matrix (m_name = simple dense, A = compressed sparse row,
  !!        serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_scsrsdenref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%dim1
             do j=1,A%iaux3(i+1)-A%iaux3(i)
                k=A%iaux3(i)+j
                m_name%dval(i,A%iaux4(k))=A%dval(k,1)
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%dim1
             do j=1,A%iaux3(i+1)-A%iaux3(i)
                k=A%iaux3(i)+j
                m_name%zval(i,A%iaux4(k))=A%zval(k,1)
             end do
          end do
       end if

  end subroutine m_copy_scsrsdenref

  !============================================================================!
  !> @brief Copy matrix (m_name = simple dense, A = compressed sparse row,
  !!        serial distribution, reference implementation with thresholding).
  !============================================================================!
  subroutine m_copy_scsrsdenref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%dim1
             do j=1,A%iaux3(i+1)-A%iaux3(i)
                k=A%iaux3(i)+j
                if (abs(A%dval(k,1))>abs_threshold) then
                   m_name%dval(i,A%iaux4(k))=A%dval(k,1)-soft_threshold*A%dval(k,1)/abs(A%dval(k,1))
                else
                   m_name%dval(i,A%iaux4(k))=0.0_dp
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%dim1
             do j=1,A%iaux3(i+1)-A%iaux3(i)
                k=A%iaux3(i)+j
                if (abs(A%zval(k,1))>abs_threshold) then
                   m_name%zval(i,A%iaux4(k))=A%zval(k,1)-soft_threshold*A%zval(k,1)/abs(A%zval(k,1))
                else
                   m_name%zval(i,A%iaux4(k))=0.0_dp
                end if
             end do
          end do
       end if

  end subroutine m_copy_scsrsdenref_thre

  !============================================================================!
  !> @brief Copy matrix (m_name = sparse coordinate list, A = compressed sparse
  !!        column, serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_scscscooref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if
       do i=1,A%dim2
          do j=1,A%iaux4(i+1)-A%iaux4(i)
             k=A%iaux4(i)+j
             m_name%iaux3(k)=A%iaux3(k)
             m_name%iaux4(k)=i
          end do
       end do

  end subroutine m_copy_scscscooref

  !============================================================================!
  !> @brief Copy matrix (m_name = sparse coordinate list, A = compressed sparse
  !!        row, serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_scsrscooref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if
       do i=1,A%dim1
          do j=1,A%iaux3(i+1)-A%iaux3(i)
             k=A%iaux3(i)+j
             m_name%iaux3(k)=i
             m_name%iaux4(k)=A%iaux4(k)
          end do
       end do

  end subroutine m_copy_scsrscooref

  !============================================================================!
  !> @brief Copy matrix (m_name = compressed sparse column, A = sparse
  !!        coordinate list, serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_scooscscref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k
    integer, allocatable :: sort_temp(:)

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%dim2+1))
       m_name%iaux4_is_allocated=.true.
       m_name%iaux4=0
       do i=1,m_name%iaux2(1)
          m_name%iaux4(A%iaux4(i)+1)=m_name%iaux4(A%iaux4(i)+1)+1
       end do
       do i=1,m_name%dim2
          m_name%iaux4(i+1)=m_name%iaux4(i+1)+m_name%iaux4(i)
       end do
       m_name%iaux3=0
       allocate(sort_temp(m_name%dim2))
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=m_name%iaux4(A%iaux4(i))+1,m_name%iaux4(A%iaux4(i)+1)
                if (m_name%iaux3(j)==0) then
                   sort_temp(A%iaux4(i))=j
                   m_name%iaux3(j)=A%iaux3(i)
                   m_name%dval(j,1)=A%dval(i,1)
                   exit
                else if (m_name%iaux3(j)>A%iaux3(i)) then
                   k=sort_temp(A%iaux4(i))
                   m_name%iaux3(j+1:k+1)=m_name%iaux3(j:k)
                   m_name%dval(j+1:k+1,1)=m_name%dval(j:k,1)
                   sort_temp(A%iaux4(i))=k+1
                   m_name%iaux3(j)=A%iaux3(i)
                   m_name%dval(j,1)=A%dval(i,1)
                   exit
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=m_name%iaux4(A%iaux4(i))+1,m_name%iaux4(A%iaux4(i)+1)
                if (m_name%iaux3(j)==0) then
                   sort_temp(A%iaux4(i))=j
                   m_name%iaux3(j)=A%iaux3(i)
                   m_name%zval(j,1)=A%zval(i,1)
                   exit
                else if (m_name%iaux3(j)>A%iaux3(i)) then
                   k=sort_temp(A%iaux4(i))
                   m_name%iaux3(j+1:k+1)=m_name%iaux3(j:k)
                   m_name%zval(j+1:k+1,1)=m_name%zval(j:k,1)
                   sort_temp(A%iaux4(i))=k+1
                   m_name%iaux3(j)=A%iaux3(i)
                   m_name%zval(j,1)=A%zval(i,1)
                   exit
                end if
             end do
          end do
       end if
       deallocate(sort_temp)

  end subroutine m_copy_scooscscref

  !============================================================================!
  !> @brief Copy matrix (m_name = compressed sparse row, A = sparse
  !!        coordinate list, serial distribution, reference implementation).
  !============================================================================!
  subroutine m_copy_scooscsrref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k
    integer, allocatable :: sort_temp(:)

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       allocate(m_name%iaux3(m_name%dim1+1))
       m_name%iaux3_is_allocated=.true.
       m_name%iaux3=0
       do i=1,m_name%iaux2(1)
          m_name%iaux3(A%iaux3(i)+1)=m_name%iaux3(A%iaux3(i)+1)+1
       end do
       do i=1,m_name%dim1
          m_name%iaux3(i+1)=m_name%iaux3(i+1)+m_name%iaux3(i)
       end do
       m_name%iaux4=0
       allocate(sort_temp(m_name%dim1))
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=m_name%iaux3(A%iaux3(i))+1,m_name%iaux3(A%iaux3(i)+1)
                if (m_name%iaux4(j)==0) then
                   sort_temp(A%iaux3(i))=j
                   m_name%iaux4(j)=A%iaux4(i)
                   m_name%dval(j,1)=A%dval(i,1)
                   exit
                else if (m_name%iaux4(j)>A%iaux4(i)) then
                   k=sort_temp(A%iaux3(i))
                   m_name%iaux4(j+1:k+1)=m_name%iaux4(j:k)
                   m_name%dval(j+1:k+1,1)=m_name%dval(j:k,1)
                   sort_temp(A%iaux3(i))=k+1
                   m_name%iaux4(j)=A%iaux4(i)
                   m_name%dval(j,1)=A%dval(i,1)
                   exit
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=m_name%iaux3(A%iaux3(i))+1,m_name%iaux3(A%iaux3(i)+1)
                if (m_name%iaux4(j)==0) then
                   sort_temp(A%iaux3(i))=j
                   m_name%iaux4(j)=A%iaux4(i)
                   m_name%zval(j,1)=A%zval(i,1)
                   exit
                else if (m_name%iaux4(j)>A%iaux4(i)) then
                   k=sort_temp(A%iaux3(i))
                   m_name%iaux4(j+1:k+1)=m_name%iaux4(j:k)
                   m_name%zval(j+1:k+1,1)=m_name%zval(j:k,1)
                   sort_temp(A%iaux3(i))=k+1
                   m_name%iaux4(j)=A%iaux4(i)
                   m_name%zval(j,1)=A%zval(i,1)
                   exit
                end if
             end do
          end do
       end if
       deallocate(sort_temp)

  end subroutine m_copy_scooscsrref

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Copy external matrix (sparse coordinate list from dense block
  !!        cyclic, parallel distribution, real version).
  !============================================================================!
  subroutine m_copy_external_pddbcpdcoo(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9)

    real(dp), intent(in) :: A(:,:)
    real(dp), intent(in), optional :: threshold

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    allocate(m_name%iaux1(size(desc)))
    m_name%iaux1_is_allocated=.true.
    m_name%iaux1=desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='coo'
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.

    call psp_den2sp_m(A,desc,m_name%spm,m_name%str_type,threshold)

    m_name%is_initialized=.true.

  end subroutine m_copy_external_pddbcpdcoo
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Copy external matrix (sparse coordinate list from dense block
  !!        cyclic, parallel distribution, complex version).
  !============================================================================!
  subroutine m_copy_external_pzdbcpzcoo(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9)

    complex(dp), intent(in) :: A(:,:)
    real(dp), intent(in), optional :: threshold

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    allocate(m_name%iaux1(size(desc)))
    m_name%iaux1_is_allocated=.true.
    m_name%iaux1=desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='coo'
    m_name%is_serial=.false.
    m_name%is_real=.false.
    m_name%is_sparse=.true.

    call psp_den2sp_m(A,desc,m_name%spm,m_name%str_type,threshold)

    m_name%is_initialized=.true.

  end subroutine m_copy_external_pzdbcpzcoo
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Copy external matrix (compressed sparse column from dense block
  !!        cyclic, parallel distribution, real version).
  !============================================================================!
  subroutine m_copy_external_pddbcpdcsc(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9)

    real(dp), intent(in) :: A(:,:)
    real(dp), intent(in), optional :: threshold

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    allocate(m_name%iaux1(size(desc)))
    m_name%iaux1_is_allocated=.true.
    m_name%iaux1=desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='csc'
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.

    call psp_den2sp_m(A,desc,m_name%spm,m_name%str_type,threshold)

    m_name%is_initialized=.true.

  end subroutine m_copy_external_pddbcpdcsc
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Copy external matrix (compressed sparse column from dense block
  !!        cyclic, parallel distribution, complex version).
  !============================================================================!
  subroutine m_copy_external_pzdbcpzcsc(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9)

    complex(dp), intent(in) :: A(:,:)
    real(dp), intent(in), optional :: threshold

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    allocate(m_name%iaux1(size(desc)))
    m_name%iaux1_is_allocated=.true.
    m_name%iaux1=desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='csc'
    m_name%is_serial=.false.
    m_name%is_real=.false.
    m_name%is_sparse=.true.

    call psp_den2sp_m(A,desc,m_name%spm,m_name%str_type,threshold)

    m_name%is_initialized=.true.

  end subroutine m_copy_external_pzdbcpzcsc
#endif

end module MatrixSwitch_m_copy
