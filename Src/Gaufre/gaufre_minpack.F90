! *** Module: gaufre_minpack ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Fitting of Natural Atomic Orbitals to Gaussians using MINPACK
!!
!! This module provides an interface to MINPACK that can be used to determine
!! the coefficients of zero-centered Gaussians which best represent
!! Natural Atomic Orbitals.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 11.2017 Created [Yann Pouillon]
!!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
!!      - 04.2018 Enhanced fitting capabilities [Yann Pouillon]
! *****************************************************************************
module gaufre_minpack

  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

  use gaufre_common

  implicit none

  private

  ! Data structure to store the fitting data and preserve its consistency
  type, public :: gaufre_minpack_t
    logical :: converged = .false.
    integer :: npts = 0
    integer :: ngfs = 0
    integer :: lwrk = 0
    real(gfdp) :: residue
    real(gfdp) :: tol = 0.0_gfdp
    integer, allocatable :: iwrk(:)
    real(gfdp), allocatable :: fvec(:)
    real(gfdp), allocatable :: gwrk(:)
    real(gfdp), allocatable :: trial(:)
  end type gaufre_minpack_t

  ! Data to fit
  real(gfdp), pointer :: xdata(:) => null()
  real(gfdp), pointer :: ydata(:) => null()

  ! Public resources
  public :: &
&   gaufre_minpack_dump, &
&   gaufre_minpack_free, &
&   gaufre_minpack_init, &
&   gaufre_minpack_wrap

contains

  ! ***************************************************************************
  !> \brief Pretty-printer for the MINPACK data structure
  !!
  !! \par History
  !!      - 11.2017 Created [Yann Pouillon]
  !!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
  !!
  !! \param[in] fit_mpck: data structure containing fitting Gaussian parameters
  ! ***************************************************************************
  subroutine gaufre_minpack_dump(fit_mpck)

    implicit none

    ! Arguments
    type(gaufre_minpack_t), intent(in) :: fit_mpck

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
&     "CONVERGED", fit_mpck%converged
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,I12)') &
&     "NGFS", fit_mpck%ngfs
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,I12)') &
&     "NPTS", fit_mpck%npts
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,I12)') &
&     "LWRK", fit_mpck%lwrk
    write(unit=output_unit, fmt='(1X,A16,1X,"=",1X,E12.5)') &
&     "TOL", fit_mpck%tol

    ! Preview allocated array components
    if ( allocated(fit_mpck%iwrk) ) then
      write(unit=output_unit, fmt='(1X,A16,1X,"= [",4(1X,I12),1X,"... ]")') &
        "IWRK(NPRM*NGFS)", fit_mpck%iwrk(1:min(4,2*fit_mpck%ngfs))
    end if
    if ( allocated(fit_mpck%fvec) ) then
      write(unit=output_unit, fmt='(1X,A16,1X,"= [",4(1X,E12.5),1X,"... ]")') &
        "FVEC(NPTS)", fit_mpck%fvec(1:4)
    end if
    if ( allocated(fit_mpck%gwrk) ) then
      write(unit=output_unit, fmt='(1X,A16,1X,"= [",4(1X,E12.5),1X,"... ]")') &
        "GWRK(LWRK)", fit_mpck%gwrk(1:4)
    end if
    if ( allocated(fit_mpck%trial) ) then
      write(unit=output_unit, fmt=*)
      write(unit=output_unit, fmt='(1X,A)') &
&       "GAUSSIAN PROJECTION PARAMETERS", &
&       "    g_{orb}(r) = SUM_{i=1}^{n} \alpha_i e^{-\beta_i r^2}"
      write(unit=output_unit, fmt=*)
      write(unit=output_unit, fmt='(5X,A5,2(1X,A20))') &
&       "INDEX", "ALPHA", "BETA"
      do igfs=1,fit_mpck%ngfs
        write(unit=output_unit, fmt='(5X,I5,2(1X,E20.8))') &
&         igfs, fit_mpck%trial(2*igfs-1), fit_mpck%trial(2*igfs)
      end do
    end if

    ! Check for improper input parameters (N .LE. 0, M .LT. N,
    ! or TOL .LT. 0.D0, or LWA .LT. M*N+5*N+M)
    data_ok = .true.
    write(unit=output_unit, fmt=*)
    write(unit=output_unit, fmt=*) "BEGIN CONSISTENCY CHECKS"
    if ( fit_mpck%ngfs .le. 0 ) then
      write(unit=output_unit, fmt=*) &
&       "  * ERROR: negative number of Gaussian functions (NGFS=", &
&       fit_mpck%ngfs, ")"
      data_ok = .false.
    end if
    if ( fit_mpck%npts .le. 2*fit_mpck%ngfs ) then
      write(unit=output_unit, fmt=*) &
&       "  * ERROR: insufficient number of points (NPTS=", &
&       fit_mpck%npts, " < ", 2*fit_mpck%ngfs, ")"
      data_ok = .false.
    end if
    if ( fit_mpck%tol .lt. 0.0_gfdp ) then
      write(unit=output_unit, fmt=*) &
&       "  * ERROR: negative tolerance (TOL=", &
&       fit_mpck%tol, ")"
      data_ok = .false.
    end if
    if ( fit_mpck%lwrk .lt. &
&        2*fit_mpck%ngfs*fit_mpck%npts + 5*2*fit_mpck%ngfs + &
&        fit_mpck%npts ) then
      write(unit=output_unit, fmt=*) &
&       "  * ERROR: insufficient work array size (LWRK=", &
&       fit_mpck%lwrk, ")"
      data_ok = .false.
    end if
    if ( data_ok ) then
      write(unit=output_unit, fmt=*) &
&       "  * All parameters OK"
    end if
    write(unit=output_unit, fmt=*) "END CONSISTENCY CHECKS"
    write(unit=output_unit, fmt=*)

  end subroutine gaufre_minpack_dump

  ! ***************************************************************************
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
  !! \param[in] fit_mpck: data structure containing fitting Gaussian parameters
  ! ***************************************************************************
  subroutine gaufre_minpack_free(fit_mpck)

    implicit none

    type(gaufre_minpack_t), intent(inout) :: fit_mpck

    fit_mpck%converged = .false.
    fit_mpck%ngfs = 0
    fit_mpck%npts = 0
    fit_mpck%lwrk = 0
    fit_mpck%tol = 0.0_gfdp
    if ( allocated(fit_mpck%iwrk) ) deallocate(fit_mpck%iwrk)
    if ( allocated(fit_mpck%fvec) ) deallocate(fit_mpck%fvec)
    if ( allocated(fit_mpck%gwrk) ) deallocate(fit_mpck%gwrk)
    if ( allocated(fit_mpck%trial) ) deallocate(fit_mpck%trial)
    xdata => null()
    ydata => null()

  end subroutine gaufre_minpack_free

  ! ***************************************************************************
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
  !! \param[out] fit_mpck: MINPACK-specific data structure
  !! \param[in] fit_data: method-agnostic initialised data structure
  ! ***************************************************************************
  subroutine gaufre_minpack_init(fit_mpck, fit_data, fit_tol)

    use gaufre_data, only: gaufre_data_t

    implicit none

    ! Arguments
    type(gaufre_minpack_t), intent(out) :: fit_mpck
    type(gaufre_data_t), intent(in) :: fit_data
    real(gfdp), intent(in) :: fit_tol

    ! Local variables
    integer :: ngfs, npts

    ! -------------------------------------------------------------------------

    ! Set convenience variables
    ngfs = fit_data%ngfs
    npts = fit_data%npts

    ! Check that internal variables are free for a new fitting procedure
    if ( allocated(fit_mpck%iwrk) ) then
      write(unit=error_unit, fmt=*) &
&       "gaufre_minpack_init: Error: indexing array already allocated"
      return
    end if
    if ( allocated(fit_mpck%fvec) ) then
      write(unit=error_unit, fmt=*) &
&       "gaufre_minpack_init: Error: residual vector already allocated"
      return
    end if
    if ( allocated(fit_mpck%gwrk) ) then
      write(unit=error_unit, fmt=*) &
&       "gaufre_minpack_init: Error: working array already allocated"
      return
    end if
    if ( allocated(fit_mpck%trial) ) then
      write(unit=error_unit, fmt=*) &
&       "gaufre_minpack_init: Error: trial vector already allocated"
      return
    end if
    if ( associated(xdata) ) then
      write(unit=error_unit, fmt=*) &
&       "gaufre_minpack_init: Error: X data vector already set"
      return
    end if
    if ( associated(ydata) ) then
      write(unit=error_unit, fmt=*) &
&       "gaufre_minpack_init: Error: Y data vector already set"
      return
    end if

    ! Check input parameters
    if ( npts < 2*2*ngfs ) then
      write(unit=error_unit, fmt=*) &
&       "gaufre_minpack_init: Warning: low number of points wrt Gaussians"
      write(unit=error_unit, fmt=*) &
&       "                           the fitting may have trouble converging"
    end if

    ! Import problem dimensions and data
    !
    ! [Note]
    !   We give the minimum value allowed by MINPACK to LWRK, that is:
    !   M*N+5*N+M, where M is NPTS and N is 2*NGFS (2 coefficients per
    !   Gaussian function).
    !
    fit_mpck%ngfs = ngfs
    fit_mpck%npts = npts
    fit_mpck%tol = fit_tol
    fit_mpck%lwrk = 2*ngfs*npts + 5*2*ngfs + npts
    allocate(fit_mpck%iwrk(2*ngfs))
    fit_mpck%iwrk(:) = 0
    allocate(fit_mpck%fvec(npts))
    fit_mpck%fvec(:) = 0.0_gfdp
    allocate(fit_mpck%gwrk(fit_mpck%lwrk))
    fit_mpck%gwrk(:) = 0.0_gfdp

    ! A vector of ones is usually a good starting guess
    allocate(fit_mpck%trial(2*ngfs))
    fit_mpck%trial(:) = 1.0_gfdp

    ! Import data to fit
    ! Note: It must be available to the fcn function called by LMDIF1
    allocate(xdata(npts))
    allocate(ydata(npts))
    xdata(:) = fit_data%x2fit(:)
    ydata(:) = fit_data%y2fit(:)

  end subroutine gaufre_minpack_init

  ! ***************************************************************************
  !> \brief Finds the best possible Gaussian expansion of the input
  !!        (x,f(x)) data
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
  !! \param[inout] fit_data: method-agnostic Gaussian-fitting data structure
  ! ***************************************************************************
  subroutine gaufre_minpack_wrap(fit_data, fit_tol)

    use gaufre_data, only: gaufre_data_t, gaufre_data_set_coeffs

    implicit none

    ! Arguments
    type(gaufre_data_t), intent(inout) :: fit_data
    real(gfdp), intent(in) :: fit_tol

    ! Local variables
    logical :: conv_ok
    integer :: info, igfs, itry
    real(gfdp) :: xguess
    type(gaufre_minpack_t) :: fit_mpck

    ! External MINPACK routines
    external :: lmdif1
    real(gfdp), external :: enorm

    ! -------------------------------------------------------------------------

    ! Init
    fit_data%coeffs(:,:) = 0.0_gfdp
    call gaufre_minpack_init(fit_mpck, fit_data, fit_tol)

    ! We may have to try more than once with different intial guesses
    ! before finding a suitable expansion (10 times at most)
    do itry=1,10

      ! Find the Gaussian expansion using the modified
      ! Levenberg-Marquardt algorithm from MINPACK
      call lmdif1(fcn, fit_mpck%npts, 2*fit_mpck%ngfs, fit_mpck%trial, &
&       fit_mpck%fvec, fit_mpck%tol, info, fit_mpck%iwrk, fit_mpck%gwrk, &
&       fit_mpck%lwrk)
      fit_data%residue = enorm(fit_mpck%npts, fit_mpck%fvec)

      ! Check convergence
      !
      ! [Note]
      !   If you want to return before the end in case of error,
      !   just make sure you call gaufre_minpack_free right before.
      !
      select case(info)
        case(0)
          conv_ok = .false.
        case(1:4)
          conv_ok = .true.
        case(5:7)
          conv_ok = .false.
        case default
          conv_ok = .false.
      end select

      ! LMDIF1 sometimes returns a successful convergence but the results
      ! are wrong, thus we check the L2-norm of the residue as well
      ! FIXME: hard-coded absolute residue tolerance
      if ( fit_data%residue > 1.0e-3_gfdp ) then
        conv_ok = .false.
      end if

      ! Propagate convergence status
      fit_mpck%converged = conv_ok
      fit_data%converged = conv_ok

      ! Dump data structure when debugging
#if ( DEBUG_LEVEL > 0 )
      call gaufre_minpack_dump(fit_mpck)
      write(unit=output_unit, fmt='(1X,A," = ",I12)') &
&       "Return value from LMDIF1", info
      write(unit=output_unit, fmt='(1X,A," = ",E12.5)') &
&       "L2-norm of the residuals", fit_data%residue
      write(unit=output_unit, fmt=*)
#endif

      ! Translate data structure into Gaussian parameters
      if ( conv_ok ) then
        call gaufre_data_set_coeffs(fit_data, fit_mpck%trial)
        exit
      end if

      ! Reset the trial vector with random values if convergence was not
      ! reached
      call random_number(fit_mpck%trial(:))
      do igfs=1,fit_data%ngfs
        call random_number(xguess)
        fit_mpck%trial(2*igfs-1) = abs(fit_mpck%trial(2*igfs-1) * xguess*itry)
        fit_mpck%trial(2*igfs) = (fit_mpck%trial(2*igfs) - 2.0_gfdp*xguess)*itry
      end do

    end do   ! itry

    ! Clean-up the mess
    call gaufre_minpack_free(fit_mpck)

    ! Display warnings if needed
    select case(info)
      case(0)
        write(unit=error_unit, fmt=*) &
&         "gaufre_minpack_find: Error: improper or unrealistic input parameters"
      case(1:3)
      case(4)
        write(unit=error_unit, fmt=*) &
&         "gaufre_minpack_find: Warning: actual convergence should be checked"
      case(5:7)
        write(unit=error_unit, fmt=*) &
&         "gaufre_minpack_find: Warning: could not converge Gaussian fit"
      case default
        write(unit=error_unit, fmt=*) &
&         "gaufre_minpack_find: Warning: unknown convergence return status", &
&         info
    end select

  end subroutine gaufre_minpack_wrap

  ! ***************************************************************************
  ! *** Private routines                                                    ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Routine called internally by LMDIF1
  !!
  !! This routine calculates the residues required by LMDIF1 to complete
  !! the multi-parameter non-linear optimisation. These residues are defined
  !! in \cite Shang2011 as follows:
  !!
  !!   \f$ \sum_i \left( \frac{\chi(r_i)}{r_i^l} - \frac{\phi_{nl}(r_i)}{r_i^l}
  !!       \right)^2 \f$
  !!
  !! For more details about LMDIF1, please see:
  !! https://www.math.utah.edu/software/minpack/minpack/lmdif1.html
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
  ! ***************************************************************************
  subroutine fcn(m, n, x, fvec, iflag)

    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

    implicit none

    integer, intent(in) :: m, n
    integer, intent(in) :: iflag
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: fvec(m)

    integer :: ipts, ivar, ngfs

    ! Make sure to return zero in case of error
    fvec(:) = 0.0d0

    ! Reference data must be initialized to compute function values
    if ( .not. (associated(xdata) .and. associated(ydata)) ) then
      write(unit=error_unit, fmt=*) &
&           "gaufre_minpack: Error: data structure is not initialized"
      !iflag = -1
      return
    end if

    ! The number of coefficients must be divisible by 2 (a, b)
    if ( modulo(n, 2) /= 0 ) then
      write(unit=error_unit, fmt=*) &
&       "gaufre_minpack: Error: the number of coefficients is not 2*p"
      !iflag = -1
      return
    end if

    ! The value of IFLAG defines what to do (should be 1)
    select case(iflag)

      ! Print current parameters
      case(0)
        write(unit=output_unit, fmt='(A)') &
&         "GAUSSIAN PROJECTION PARAMETERS: ", &
&         "g_{orb}(r) = SUM_{i=1}^{n} \beta_i e^{-\alpha_i r^2}"
        write(unit=output_unit, fmt=*)
        write(unit=output_unit, fmt='(A5,2(1X,A20))') &
&         "INDEX", "ALPHA", "BETA"

        ngfs = n/2
        do ivar=1,ngfs
          write(unit=output_unit, fmt='(I5,2(1X,E20.8))') &
&           ivar, x(2*ivar-1), x(2*ivar)
        end do
        write(unit=output_unit, fmt=*)

      ! Calculate residue
      case(1:2)
        fvec(:) = 0.0d0
        do ipts=1,m
          if ( xdata(ipts) > 1.0d-13 ) then
            do ivar=1,n-1,2
              fvec(ipts) = fvec(ipts) &
&                 + x(ivar+1) &
&                 * exp(-abs(x(ivar))*(xdata(ipts)**2))
            end do
            fvec(ipts) = (fvec(ipts) - ydata(ipts))
            fvec(ipts) = fvec(ipts) * fvec(ipts)
          end if
        end do

      case default
        write(unit=error_unit, fmt=*) &
&           "gaufre_minpack: Error: unsupported value IFLAG = ", iflag
        return

    end select

  end subroutine fcn

end module gaufre_minpack
