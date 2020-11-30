! *** Module: nao2gto_fit ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Fitting of Natural Atomic Orbitals to Gaussians
!!
!! This modules provides tools to determine the coefficients of
!! zero-centered Gaussians which best represent Natural Atomic Orbitals.
!!
!! \note This module is fully independent from SIESTA, in order to allow
!!       its testing in several contexts.
!! 
!! \author Yann Pouillon
!!
!! \par History
!!      - 11.2017 Created [Yann Pouillon]
!!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
! *****************************************************************************
module nao2gto_fit

  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

  implicit none

  private

  integer, parameter :: dp = kind(0.0d0)

  public :: &
&   nao2gto_fit_dump, &
&   nao2gto_fit_eval, &
&   nao2gto_fit_find

  ! Data structure to store the fitting data and preserve its consistency
  type :: nao2gto_fit_t
    logical :: converged = .false.
    integer :: npts = 0
    integer :: ngfs = 0
    integer :: nprm = 2
    integer :: lwrk = 0
    real(dp) :: tol = 0.0_dp
    integer, allocatable :: iwrk(:)
    real(dp), allocatable :: fvec(:)
    real(dp), allocatable :: gwrk(:)
    real(dp), allocatable :: trial(:)
  end type

  ! Data to fit
  real(dp), allocatable, save :: xdata(:)
  real(dp), allocatable, save :: ydata(:)

contains

! *****************************************************************************
!> \brief Pretty-printer for the whole data structure
!!
!! \par History
!!      - 11.2017 Created [Yann Pouillon]
!!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
!!
!! \param[in] gdata: data structure containing fitting Gaussian parameters
! *****************************************************************************
  subroutine nao2gto_fit_dump(gdata)

    implicit none

    ! Arguments
    type(nao2gto_fit_t), intent(in) :: gdata

    ! Local variables
    logical :: data_ok
    integer :: igfs

    ! -------------------------------------------------------------------------

    ! Banner
    write(unit=output_unit, fmt=*) &
&     "DUMPING GAUSSIAN FIT DATA STRUCTURE"
    write(unit=output_unit, fmt=*)

    ! Show scalar components
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,L12)') &
&     "CONVERGED", gdata%converged
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,I12)') &
&     "NGFS", gdata%ngfs
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,I12)') &
&     "NPRM", gdata%nprm
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,I12)') &
&     "NPTS", gdata%npts
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,I12)') &
&     "LWRK", gdata%lwrk
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,E12.5)') &
&     "TOL", gdata%tol

    ! Preview allocated array components
    if ( allocated(gdata%iwrk) ) then
      write(unit=output_unit, fmt='(1X,A16,1X,"= [",4(1X,I12),1X,"... ]")') &
        "IWRK(NPRM*NGFS)", gdata%iwrk(1:4)
    end if
    if ( allocated(gdata%fvec) ) then
      write(unit=output_unit, fmt='(1X,A16,1X,"= [",4(1X,E12.5),1X,"... ]")') &
        "FVEC(NPTS)", gdata%fvec(1:4)
    end if
    if ( allocated(gdata%gwrk) ) then
      write(unit=output_unit, fmt='(1X,A16,1X,"= [",4(1X,E12.5),1X,"... ]")') &
        "GWRK(LWRK)", gdata%gwrk(1:4)
    end if
    if ( allocated(gdata%trial) ) then
      write(unit=output_unit, fmt=*)
      write(unit=output_unit, fmt='(1X,A)') &
&       "GAUSSIAN PROJECTION PARAMETERS", &
&       "    g_{orb}(r) = SUM_{i=1}^{n} \alpha_i e^{-\beta_i r^2}"
      write(unit=output_unit, fmt=*)
      write(unit=output_unit, fmt='(5X,A5,2(1X,A20))') &
&       "INDEX", "ALPHA", "BETA"
      do igfs=1,gdata%ngfs
        write(unit=output_unit, fmt='(5X,I5,2(1X,E20.8))') &
&         igfs, gdata%trial(2*igfs-1), gdata%trial(2*igfs)
      end do
    end if

    ! Check for improper input parameters (N .LE. 0, M .LT. N,
    ! or TOL .LT. 0.D0, or LWA .LT. M*N+5*N+M)
    data_ok = .true.
    write(unit=output_unit, fmt=*)
    write(unit=output_unit, fmt=*) "BEGIN CONSISTENCY CHECKS"
    if ( gdata%nprm .ne. 2 ) then
      write(unit=output_unit, fmt=*) &
&       "  * ERROR: invalid number of Gaussian parameters (NPRM=", &
&       gdata%nprm, ", should be 2)"
      data_ok = .false.
    end if
    if ( gdata%ngfs .le. 0 ) then
      write(unit=output_unit, fmt=*) &
&       "  * ERROR: negative number of Gaussian functions (NGFS=", &
&       gdata%ngfs, ")"
      data_ok = .false.
    end if
    if ( gdata%npts .le. gdata%nprm*gdata%ngfs ) then
      write(unit=output_unit, fmt=*) &
&       "  * ERROR: insufficient number of points (NPTS=", &
&       gdata%npts, " < ", gdata%nprm*gdata%ngfs, ")"
      data_ok = .false.
    end if
    if ( gdata%tol .lt. 0.0_dp ) then
      write(unit=output_unit, fmt=*) &
&       "  * ERROR: negative tolerance (TOL=", &
&       gdata%tol, ")"
      data_ok = .false.
    end if
    if ( gdata%lwrk .lt. &
&        gdata%nprm*gdata%ngfs*gdata%npts + 5*gdata%nprm*gdata%ngfs + &
&        gdata%npts ) then
      write(unit=output_unit, fmt=*) &
&       "  * ERROR: insufficient work array size (LWRK=", &
&       gdata%lwrk, ")"
      data_ok = .false.
    end if
    if ( data_ok ) then
      write(unit=output_unit, fmt=*) &
&       "  * All parameters OK"
    end if
    write(unit=output_unit, fmt=*) "END CONSISTENCY CHECKS"
    write(unit=output_unit, fmt=*)

  end subroutine nao2gto_fit_dump

! *****************************************************************************
!> \brief Evaluates a Gaussian expansion on a radial grid
!!
!! \par History
!!      - 11.2017 Created [Yann Pouillon]
!!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
!!
!! \param[in] ngfs: number of Gaussian functions
!! \param[in] npts: number of grid points
!! \param[in] coeffs: coefficients of the Gaussian expansion
!! \param[in] rad: points of the radial grid
!! \param[out] g_orb: values of the Gaussian expansion on the radial grid
! *****************************************************************************
  subroutine nao2gto_fit_eval(ngfs, npts, coeffs, rad, g_orb)

    implicit none

    ! Arguments
    integer, intent(in) :: ngfs, npts
    real(dp), intent(in) :: coeffs(2,ngfs)
    real(dp), intent(in) :: rad(npts)
    real(dp), intent(out) :: g_orb(npts)

    ! Local variables
    integer :: igfs, ipts

    ! -------------------------------------------------------------------------

    g_orb(:) = 0.0_dp
    do ipts=1,npts
      do igfs=1,ngfs
        g_orb(ipts) = g_orb(ipts) &
&         + coeffs(1,igfs) &
&         * exp(-coeffs(2,igfs)*(rad(ipts)**2))
      end do
    end do

  end subroutine nao2gto_fit_eval

! *****************************************************************************
!> \brief Finds the best possible Gaussian expansion of the input (x,f(x)) data
!!
!! This routine attempts to fit a real-space function with a sum of
!! zero-centered Gaussians using the LMDIF1 routine of MINPACK.
!!
!! The fitting is realized through a Levenberg-Marquardt procedure for
!! nonlinear optimization.
!!
!! \par History
!!      - 11.2017 Created [Yann Pouillon]
!!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
!!
!! \param[in] ngfs: number of Gaussian functions
!! \param[in] npts: number of grid points
!! \param[in] rad: points of the radial grid
!! \param[in] orb: values of the orbital to expand on the radial grid
!! \param[out] coeffs: coefficients of the Gaussian expansion
!! \param[out] residual: residual error after the fitting procedure
!! \return a boolean value telling whether the fitting procedure converged
! *****************************************************************************
  function nao2gto_fit_find(ngfs, npts, rad, orb, coeffs, residual) &
&   result(conv_ok)

    use minpack, only: enorm

    implicit none

    ! Arguments
    integer, intent(in) :: ngfs, npts
    real(dp), intent(in) :: rad(npts), orb(npts)
    real(dp), intent(out) :: coeffs(2,ngfs), residual

    ! Local variables
    logical :: conv_ok
    character(len=1024) :: dbgmsg
    integer :: info, igfs
    type(nao2gto_fit_t) :: gdata

    ! -------------------------------------------------------------------------

    ! Init
    coeffs(:,:) = 0.0_dp
    call gauss_data_init(gdata, ngfs, npts, rad, orb)

    ! Find the Gaussian expansion using the MINPACK modified
    ! Levenberg-Marquardt algorithm
    call lmdif1(fcn, gdata%npts, gdata%nprm*gdata%ngfs, gdata%trial, &
&     gdata%fvec, gdata%tol, info, gdata%iwrk, gdata%gwrk, gdata%lwrk)
    residual = enorm(gdata%npts, gdata%fvec)

    ! Check convergence
    !
    ! [Note]
    !   If you want to return before the end in case of error,
    !   just make sure you call gauss_data_free right before.
    !
    select case(info)
      case(0)
        conv_ok = .false.
        write(unit=error_unit, fmt=*) &
&         "nao2gto_fit_find: Error: improper or unrealistic input parameters"
      case(1:3)
        conv_ok = .true.
      case(4)
        conv_ok = .true.
        write(unit=error_unit, fmt=*) &
&         "nao2gto_fit_find: Warning: actual convergence should be checked"
      case(5:7)
        conv_ok = .false.
        write(unit=error_unit, fmt=*) &
&         "nao2gto_fit_find: Warning: could not converge Gaussian fit"
      case default
        conv_ok = .false.
        write(unit=error_unit, fmt=*) &
&         "nao2gto_fit_find: Warning: unknown convergence return status"
    end select
    gdata%converged = conv_ok

    ! Translate data structure into Gaussian parameters
    if ( conv_ok ) then
      do igfs=1,ngfs
        coeffs(1,igfs) = gdata%trial(2*igfs-1)
        coeffs(2,igfs) = gdata%trial(2*igfs)
      end do
    end if

    ! Dump data structure
    call nao2gto_fit_dump(gdata)
    write(unit=output_unit, fmt='(1X,A," = ",I12)') &
&     "Return value from LMDIF1", info
    write(unit=output_unit, fmt='(1X,A," = ",E12.5)') &
&     "L2-norm of the residuals", residual

    ! Clean-up the mess
    call gauss_data_free(gdata)

  end function nao2gto_fit_find

! *****************************************************************************
! *** Private routines                                                      ***
! *****************************************************************************

! *****************************************************************************
!> \brief Constructor for the Gaussian expansion data structure
!!
!! This routine allocates all the necessary memory and initializes the
!! corresponding parameters before running a Gaussian fitting procedure.
!! It must be called before using the specified data structure.
!!
!! \note We want to be able to use subsets of grid points and orbital values
!!       not necessarily including boundaries and without interfering with
!!       the original variables, this is why we explictly copy them for now.
!!
!! \par History
!!      - 11.2017 Created [Yann Pouillon]
!!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
!!
!! \param[out] gdata: data structure containing fitting Gaussian parameters
!! \param[in] ngfs: number of desired Gaussian functions
!! \param[in] npts: number of desired radial grid points
!! \param[in] x2fit: grid points of the data to fit
!! \param[in] y2fit: values of the data to fit at x2fit
! *****************************************************************************
  subroutine gauss_data_init(gdata, ngfs, npts, x2fit, y2fit)

    implicit none

    ! Arguments
    integer, intent(in)              :: ngfs, npts
    real(dp), intent(in)             :: x2fit(npts), y2fit(npts)
    type(nao2gto_fit_t), intent(out) :: gdata

    ! The current implementation only works with 2 coefficients per Gaussian
    integer, parameter :: nprm = 2
    integer :: igfs, ipts

    ! -------------------------------------------------------------------------

    ! Check that internal variables are free for a new fitting procedure
    if ( allocated(gdata%iwrk) ) then
      write(unit=error_unit, fmt=*) &
&       "gauss_data_init: Error: indexing array already allocated"
      return
    end if
    if ( allocated(gdata%fvec) ) then
      write(unit=error_unit, fmt=*) &
&       "gauss_data_init: Error: residual vector already allocated"
      return
    end if
    if ( allocated(gdata%gwrk) ) then
      write(unit=error_unit, fmt=*) &
&       "gauss_data_init: Error: working array already allocated"
      return
    end if
    if ( allocated(gdata%trial) ) then
      write(unit=error_unit, fmt=*) &
&       "gauss_data_init: Error: trial vector already allocated"
      return
    end if
    if ( allocated(xdata) ) then
      write(unit=error_unit, fmt=*) &
&       "gauss_data_init: Error: X data vector already allocated"
      return
    end if
    if ( allocated(ydata) ) then
      write(unit=error_unit, fmt=*) &
&       "gauss_data_init: Error: Y data vector already allocated"
      return
    end if

    ! Check input parameters
    if ( npts < 2 * nprm * ngfs ) then
      write(unit=error_unit, fmt=*) &
&       "gauss_data_init: Warning: low number of points wrt Gaussians"
      write(unit=error_unit, fmt=*) &
&       "                           the fitting may have trouble converging"
    end if

    ! Set tolerance to the square root of the machine precision (unless high
    ! precision solutions are required, this is the recommended setting)
    !gdata%tol = sqrt(epsilon(1.0_dp))

    ! FIXME: need to increase tolerance quite a lot to get a result
    !gdata%tol = 5.0e-3_dp
    gdata%tol = 5.0e-4_dp

    ! Import problem dimensions and data
    !
    ! [Note]
    !   We give the minimum value allowed by MINPACK to LWRK, that is:
    !   M*N+5*N+M, where M is NPTS and N is 2*NGFS (2 coefficients per
    !   Gaussian function).
    !
    gdata%ngfs = ngfs
    gdata%nprm = nprm
    gdata%npts = npts
    gdata%lwrk = nprm*ngfs*npts + 5*nprm*ngfs + npts
    allocate(gdata%iwrk(nprm*ngfs))
    gdata%iwrk(:) = 0
    allocate(gdata%fvec(gdata%npts))
    gdata%fvec(:) = 0.0_dp
    allocate(gdata%gwrk(gdata%lwrk))
    gdata%gwrk(:) = 0.0_dp

    ! A vector of ones is usually a good starting guess
    allocate(gdata%trial(nprm*ngfs))
    call random_number(gdata%trial(:))
    do igfs=1,nprm*ngfs,2
      gdata%trial(igfs) = 2.0_dp * gdata%trial(igfs) - 1.0_dp
    end do

    ! Import data to fit
    ! Note: It must be available to the fcn function called by LMDIF1
    allocate(xdata(npts))
    allocate(ydata(npts))
    do ipts=1,npts
      xdata(ipts) = x2fit(ipts)
      ydata(ipts) = y2fit(ipts)
    end do

  end subroutine gauss_data_init

! *****************************************************************************
!> \brief Destructor for the Gaussian expansion data structure
!!
!! This routine frees all memory allocated for the Gaussian fitting procedure
!! and should be called once the data structure is not needed anymore.
!!
!! \note This routine is idempotent and can be run as many times as desired.
!!
!! \par History
!!      - 11.2017 Created [Yann Pouillon]
!!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
!!
!! \param[in] gdata: data structure containing fitting Gaussian parameters
! *****************************************************************************
  subroutine gauss_data_free(gdata)

    implicit none

    type(nao2gto_fit_t), intent(inout) :: gdata

    gdata%converged = .false.
    gdata%ngfs = 0
    gdata%nprm = 2
    gdata%npts = 0
    gdata%lwrk = 0
    gdata%tol = 0.0_dp
    if ( allocated(gdata%iwrk) ) deallocate(gdata%iwrk)
    if ( allocated(gdata%fvec) ) deallocate(gdata%fvec)
    if ( allocated(gdata%gwrk) ) deallocate(gdata%gwrk)
    if ( allocated(gdata%trial) ) deallocate(gdata%trial)
    if ( allocated(xdata) ) deallocate(xdata)
    if ( allocated(ydata) ) deallocate(ydata)

  end subroutine gauss_data_free

! *****************************************************************************
!> \brief Routine called internally by LMDIF1
!!
!! Details at: https://www.math.utah.edu/software/minpack/minpack/lmdif1.html
!!
!! \par History
!!      - 11.2017 Created [Yann Pouillon]
!!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
!!
!! \note The interface of this routine is imposed by MINPACK.
!!
!! \param[in] M: number of data points (npts in the rest of this module)
!! \param[in] N: number of variables (2*ngfs in the rest of this module)
!! \param[in] X: trial set of variables
!! \param[in,out] FVEC: computed residuals
!! \param[in,out] IFLAG: controls what to compute (1=residuals,
!!                       2=Jacobian residuals)
! *****************************************************************************
  subroutine fcn(m, n, x, fvec, iflag)

    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

    implicit none

    integer, intent(in) :: m, n
    integer, intent(inout) :: iflag
    double precision, intent(in) :: x(n)
    double precision, intent(inout) :: fvec(m)

    integer :: ipts, ivar, ngfs

    ! Make sure to return zero in case of error
    fvec(:) = 0.0d0

    ! Reference data must be initialized to compute function values
    if ( .not. (allocated(xdata) .and. allocated(ydata)) ) then
      write(unit=error_unit, fmt=*) &
&           "nao2gto_fit: Error: data structure is not initialized"
      iflag = -1
      return
    end if

    ! The number of coefficients must be divisible by 2 (a, b)
    if ( modulo(n, 2) /= 0 ) then
      write(unit=error_unit, fmt=*) &
&       "nao2gto_fit: Error: the number of coefficients is not 2*p"
      iflag = -1
      return
    end if

    ! The value of IFLAG defines what to do (should be 1)
    select case(iflag)

      case(0)
        write(unit=output_unit, fmt='(A)') &
&         "GAUSSIAN PROJECTION PARAMETERS: ", &
&         "g_{orb}(r) = SUM_{i=1}^{n} \alpha_i e^{-\beta_i r^2}"
        write(unit=output_unit, fmt=*)
        write(unit=output_unit, fmt='(A5,2(1X,A20))') &
&         "INDEX", "ALPHA", "BETA"

        ngfs = n/2
        do ivar=1,ngfs
          write(unit=output_unit, fmt='(I5,2(1X,E20.8))') &
&           ivar, x(2*ivar-1), x(2*ivar)
        end do
        write(unit=output_unit, fmt=*)

      case(1:2)
        do ipts=1,m
          do ivar=1,n-1,2
            fvec(ipts) = fvec(ipts) &
&               + x(ivar) * exp(-x(ivar+1)*(xdata(ipts)**2))
          end do
          fvec(ipts) = fvec(ipts) - ydata(ipts)
        end do

      case default
        write(unit=error_unit, fmt=*) &
&           "nao2gto_fit: Error: unsupported value IFLAG = ", iflag
        return

    end select

  end subroutine fcn

end module nao2gto_fit
