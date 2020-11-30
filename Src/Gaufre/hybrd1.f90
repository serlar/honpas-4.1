subroutine hybrd1 ( fcn, n, x, fvec, tol, info )

!*****************************************************************************80
!
!! HYBRD1 seeks a zero of N nonlinear equations in N variables.
!
!  Discussion:
!
!    HYBRD1 finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  This is done by using the
!    more general nonlinear equation solver HYBRD.  The user must provide a
!    subroutine which calculates the functions.  The jacobian is then
!    calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    19 August 2016
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
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( n, x, fvec, iflag )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fvec(n)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(N), the functions evaluated at the output X.
!
!    Input, real ( kind = 8 ) TOL.  Termination occurs when the algorithm
!    estimates that the relative error between X and the solution is at
!    most TOL.  TOL should be nonnegative.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the
!    description of FCN.
!    Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    2, number of calls to FCN has reached or exceeded 200*(N+1).
!    3, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, the iteration is not making good progress.
!
  implicit none

  integer ( kind = 4 ) lwa
  integer ( kind = 4 ) n

  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) epsfcn
  real ( kind = 8 ) factor
  external fcn
  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) nprint
  real ( kind = 8 ) qtf(n)
  real ( kind = 8 ) r((n*(n+1))/2)
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtol

  if ( n <= 0 ) then
    info = 0
    return
  end if

  if ( tol < 0.0D+00 ) then
    info = 0
    return
  end if

  xtol = tol
  maxfev = 200 * ( n + 1 )
  ml = n - 1
  mu = n - 1
  epsfcn = 0.0D+00
  diag(1:n) = 1.0D+00
  mode = 2
  factor = 100.0D+00
  nprint = 0
  info = 0
  nfev = 0
  fjac(1:n,1:n) = 0.0D+00
  ldfjac = n
  r(1:(n*(n+1))/2) = 0.0D+00
  lr = ( n * ( n + 1 ) ) / 2
  qtf(1:n) = 0.0D+00

  call hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
    factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf )

  if ( info == 5 ) then
    info = 4
  end if

  return
end subroutine hybrd1
