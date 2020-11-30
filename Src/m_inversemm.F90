MODULE m_inversemm
USE precision
USE MatrixSwitch
!
IMPLICIT NONE
PRIVATE
PUBLIC      :: inversemm
!
CONTAINS
! 
SUBROUTINE inversemm(C,D)
!
! This subroutine performs matrix inversion and multiplication
! at the same time by invoking scalapack subroutine pzgesv that solves the 
! A*X = B problem. This trick saves time and memory. 
! Previously inversion and multiplication were being performed 
! Separately.
! Written by Rafi Ullah, 2017
!
IMPLICIT NONE
!**************** INPUT ***********************************!
!
TYPE(matrix), INTENT(INOUT) :: C  
! C: Matrix to be inverted before multiplication with D.
!
!**************** INOUT ***********************************!
!
TYPE(matrix), INTENT(INOUT) :: D 
! D: Matrix to be multiplied with inverse of C
! D=C^-1*D
!
!**************** INTERNAL ********************************!
!
INTEGER, ALLOCATABLE :: ipiv(:)
INTEGER              :: info
 !
 IF (.NOT. C%is_square) CALL die('ERROR: inversemm: matrix C is not square')
 IF (C%dim2 .NE. D%dim1) CALL die('ERROR: inversemm: order mismatch')
 !
#ifdef MPI      
 !
 ALLOCATE(ipiv(C%iaux1(3)+C%iaux1(5)))
 CALL pzgesv(C%dim1,D%dim2,C%zval,1,1,C%iaux1,ipiv,D%zval,1,1,D%iaux1,info)
 IF (info .NE. 0) CALL die('ERROR: error in pzgesv')
#else
 ALLOCATE(ipiv(C%dim1))
 CALL zgesv(C%dim1,D%dim2,C%zval,C%dim1,ipiv,D%zval,D%dim1,info)
 IF (info .NE. 0) CALL die('ERROR: error in zgsesv')
#endif

 DEALLOCATE(ipiv)
!
END SUBROUTINE inversemm 
!
END MODULE m_inversemm
