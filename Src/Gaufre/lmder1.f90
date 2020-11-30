subroutine lmder1 ( fcn, m, n, x, fvec, fjac, ldfjac, tol, info )

!*****************************************************************************80
!
!! LMDER1 minimizes M functions in N variables by Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDER1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    This is done by using the more general least-squares solver LMDER.
!    The user must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
!      integer ( kind = 4 ) ldfjac
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjac(ldfjac,n)
!      real ( kind = 8 ) fvec(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.  
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!      P' * ( JAC' * JAC ) * P = R' * R,
!    where P is a permutation matrix and JAC is the final calculated
!    jacobian.  Column J of P is column IPVT(J) of the identity matrix.
!    The lower trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, is the leading dimension of FJAC,
!    which must be no less than M.
!
!    Input, real ( kind = 8 ) TOL.  Termination occurs when the algorithm
!    estimates either that the relative error in the sum of squares is at
!    most TOL or that the relative error between X and the solution is at
!    most TOL.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares is
!       possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) factor
  external fcn
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) ftol
  real ( kind = 8 ) fvec(m)
  real ( kind = 8 ) gtol
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) njev
  integer ( kind = 4 ) nprint
  real ( kind = 8 ) qtf(n)
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtol

  info = 0

  if ( n <= 0 ) then
    return
  else if ( m < n ) then
    return
  else if ( ldfjac < m ) then
    return
  else if ( tol < 0.0D+00 ) then
    return
  end if

  factor = 100.0D+00
  maxfev = 100 * ( n + 1 )
  ftol = tol
  xtol = tol
  gtol = 0.0D+00
  mode = 1
  nprint = 0

  call lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
    diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

  if ( info == 8 ) then
    info = 4
  end if

  return
end subroutine lmder1
