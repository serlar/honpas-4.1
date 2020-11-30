! *** Module: gaufre_data ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Data types for the fitting of Natural Atomic Orbitals to Gaussians
!!
!! This module defines algorithm-independent data types and utilities
!! to fit Natural Atomic Orbitals with Gaussians.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 04.2018 Created [Yann Pouillon]
! *****************************************************************************
module gaufre_data

  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

  use gaufre_common, only: gfdp

  implicit none

  private

  ! ***************************************************************************
  !> \brief Data structure to store the fitting parameters
  !!
  !! This data type stores the information needed to perform and assess
  !! the fitting of a function by a sum of zero-centered Gaussians. The
  !! user provides the number of desired Gaussian functions (NGFS) and
  !! the number of data points (NPTS), as well as the positions of the data
  !! points (X2FIT) and their values (Y2FIT). The results are then stored
  !! in the COEFFS(2,NGFS) array along with their convergence status
  !! (CONVERGED) and the residual norm of the fitting (RESIDUE). The lower
  !! the residue, the better.
  !!
  !! COEFFS(1,:) stores the exponents of the Gaussians, while COEFFS(2,:)
  !! contains the corresponding scaling factors. The fitted data hence reads:
  !! \f[ f(r) = \sum_{i=1}^{NGFS} COEFFS(2,i) \times
  !!            e^{-COEFFS(1,i) \cdot r^2} \f]
  !!
  ! ***************************************************************************
  type, public :: gaufre_data_t
    logical :: converged = .false.   !< Whether the fitting procedure converged
    integer :: ngfs = 0   !< Number of Gaussian functions
    integer :: npts = 0   !< Number of data points to fit
    real(gfdp) :: residue = 0.0_gfdp   !< Residual norm of the fit
    real(gfdp), pointer :: coeffs(:,:) => null()   !< Gaussian coefficients
    real(gfdp), pointer :: x2fit(:) => null()   !< Mesh points
    real(gfdp), pointer :: y2fit(:) => null()   !< Data values at mesh points
  end type gaufre_data_t

  ! Public resources
  public :: &
&   gaufre_data_dump, &
&   gaufre_data_eval, &
&   gaufre_data_free, &
&   gaufre_data_init, &
&   gaufre_data_set_coeffs, &
&   gaufre_data_set_trial

contains

  ! ***************************************************************************
  ! *** Public routines                                                     ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Pretty-printer for the gaufre_data_t data structure
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[in] fit_data: data structure containing fitting Gaussian parameters
  ! ***************************************************************************
  subroutine gaufre_data_dump(fit_data)

    implicit none

    ! Arguments
    type(gaufre_data_t), intent(in) :: fit_data

    ! Local variables
    logical :: data_ok
    integer :: igfs

    ! -------------------------------------------------------------------------

    ! Banner
    write(unit=output_unit, fmt='(A)') "GAUFRE DATA STRUCTURE BEGIN"
    write(unit=output_unit, fmt='(A)')

    ! Show scalar components
    write(unit=output_unit, fmt='(2X,"* ",A16,1X,"=",1X,L12)') &
&     "CONVERGED", fit_data%converged
    write(unit=output_unit, fmt='(2X,"* ",A16,1X,"=",1X,I12)') &
&     "NGFS", fit_data%ngfs
    write(unit=output_unit, fmt='(2X,"* ",A16,1X,"=",1X,I12)') &
&     "NPTS", fit_data%npts
    write(unit=output_unit, fmt='(2X,"* ",A16,1X,"=",1X,I12)') &
&     "RESIDUE", fit_data%residue

    ! Preview allocated array components
    if ( associated(fit_data%x2fit) ) then
      write(unit=output_unit, fmt='(1X,A16,1X,"= [",4(1X,E12.5),1X,"... ]")') &
        "X2FIT(NPTS)", fit_data%x2fit(1:4)
    else
      write(unit=output_unit, fmt='(1X,A16,1X,"= NULL")') "X2FIT"
    end if
    if ( associated(fit_data%y2fit) ) then
      write(unit=output_unit, fmt='(1X,A16,1X,"= [",4(1X,E12.5),1X,"... ]")') &
        "Y2FIT(NPTS)", fit_data%y2fit(1:4)
    else
      write(unit=output_unit, fmt='(1X,A16,1X,"= NULL")') "Y2FIT"
    end if

    ! Show Gaussian fitting coefficients
    if ( associated(fit_data%coeffs) ) then
      write(unit=output_unit, fmt='(A)')
      write(unit=output_unit, fmt='(2X,A)') &
&       "GAUSSIAN FITTING COEFFICIENTS (COEFFS)", &
&       "    fit(r) = SUM_{i=1}^{n} \beta_i e^{-\alpha_i r^2}"
      write(unit=output_unit, fmt='(A)')
      write(unit=output_unit, fmt='(4X,A5,2(1X,A20))') &
&       "INDEX", "ALPHA", "BETA"
      do igfs=1,fit_data%ngfs
        write(unit=output_unit, fmt='(4X,I5,2(1X,E20.8))') &
&         igfs, fit_data%coeffs(1,igfs), fit_data%coeffs(2,igfs)
      end do
    else
      write(unit=output_unit, fmt='(4X,A16,1X,"= NULL")') "COEFFS"
    end if

    ! Check for improper input parameters (NGFS .LE. 0, NPTS .LT. 2*NGFS)
    data_ok = .true.
    write(unit=output_unit, fmt='(A)')
    write(unit=output_unit, fmt='(2X,A)') "CONSISTENCY CHECKS BEGIN"
    if ( fit_data%ngfs <= 0 ) then
      write(unit=output_unit, fmt='(4X,"* ",A)') &
&       "ERROR: negative number of Gaussian functions (NGFS=", &
&       fit_data%ngfs, ")"
      data_ok = .false.
    end if
    if ( fit_data%npts <= 2*fit_data%ngfs ) then
      write(unit=output_unit, fmt='(4X,"* ",A)') &
&       "ERROR: insufficient number of points (NPTS=", &
&       fit_data%npts, " < ", 2*fit_data%ngfs, ")"
      data_ok = .false.
    end if
    if ( data_ok ) then
      write(unit=output_unit, fmt='(4X,"* ",A)') &
&       "All parameters OK"
    end if
    write(unit=output_unit, fmt='(2X,A)') "CONSISTENCY CHECKS END"

    ! Footer
    write(unit=output_unit, fmt='(A)')
    write(unit=output_unit, fmt='(A)') "GAUFRE DATA STRUCTURE END"

  end subroutine gaufre_data_dump

  ! ***************************************************************************
  !> \brief Evaluates a Gaussian expansion on a radial mesh
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[in] fit_data: data structure containing the Gaussian expansion
  !! \param[in] rad_mesh: points of the radial mesh
  !! \param[out] fit_values: values of the Gaussian expansion on the radial
  !!                         mesh
  ! ***************************************************************************
  subroutine gaufre_data_eval(fit_data, rad_mesh, fit_values)

    implicit none

    ! Arguments
    type(gaufre_data_t), intent(in) :: fit_data
    real(gfdp), intent(in) :: rad_mesh(fit_data%npts)
    real(gfdp), intent(out) :: fit_values(fit_data%npts)

    ! Local variables
    integer :: igfs, ipts

    ! -------------------------------------------------------------------------

    fit_values(:) = 0.0_gfdp
    do ipts=1,fit_data%npts
      do igfs=1,fit_data%ngfs
        fit_values(ipts) = fit_values(ipts) &
&         + fit_data%coeffs(2,igfs) &
&         * exp(-fit_data%coeffs(1,igfs)*(rad_mesh(ipts)**2))
      end do
    end do

  end subroutine gaufre_data_eval

  ! ***************************************************************************
  !> \brief Destructor for the Gaussian expansion data structure
  !!
  !! This routine frees all memory allocated for the Gaussian fitting procedure
  !! and should be called once the data structure is not needed anymore.
  !!
  !! \note This routine is idempotent and can be run as many times as desired.
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[in] fit_data: data structure containing fitting Gaussian parameters
  ! ***************************************************************************
  subroutine gaufre_data_free(fit_data)

    implicit none

    ! Arguments
    type(gaufre_data_t), intent(inout) :: fit_data

    ! -------------------------------------------------------------------------

    ! Free allocated arrays
    if ( associated(fit_data%coeffs) ) then
      deallocate(fit_data%coeffs)
      fit_data%coeffs => null()
    end if
    if ( associated(fit_data%x2fit) ) then
      deallocate(fit_data%x2fit)
      fit_data%x2fit => null()
    end if
    if ( associated(fit_data%y2fit) ) then
      deallocate(fit_data%y2fit)
      fit_data%y2fit => null()
    end if

    ! Reset scalar variables
    fit_data%converged = .false.
    fit_data%ngfs = 0
    fit_data%npts = 0
    fit_data%residue = 0.0_gfdp

  end subroutine gaufre_data_free

  ! ***************************************************************************
  !> \brief Constructor for the Gaussian expansion data structure
  !!
  !! This routine allocates all the necessary memory and initializes the
  !! corresponding parameters before running a Gaussian fitting procedure.
  !! It must be called before using the specified data structure.
  !!
  !! \note We want to be able to use subsets of mesh points and orbital values
  !!       not necessarily including boundaries and without interfering with
  !!       the original variables, this is why we explictly copy them for now.
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[out] fit_data: data structure containing Gaussian fitting
  !!                       parameters
  !! \param[in] ngfs: number of desired Gaussian functions
  !! \param[in] npts: number of desired radial mesh points
  !! \param[in] x2fit: mesh points of the data to fit
  !! \param[in] y2fit: values of the data to fit at x2fit
  ! ***************************************************************************
  subroutine gaufre_data_init(fit_data, ngfs, npts, x2fit, y2fit)

    implicit none

    ! Arguments
    integer, intent(in) :: ngfs, npts
    real(gfdp), intent(in) :: x2fit(npts), y2fit(npts)
    type(gaufre_data_t), intent(out) :: fit_data

    ! Local variables
    logical :: data_ok

    ! -------------------------------------------------------------------------

    ! Check data consistency
    ! Note: catching all errors in one pass before returning
    data_ok = .true.
    if ( ngfs <= 0 ) then
      write(unit=error_unit, fmt=*) &
&       "Error: insufficient number of Gaussian functions (NGFS=", ngfs, ")"
      data_ok = .false.
    end if
    if ( npts <= 2*ngfs ) then
      write(unit=error_unit, fmt=*) &
&       "Error: insufficient number of points (NPTS=", npts, " < ", 2*ngfs, ")"
      data_ok = .false.
    end if
    if ( size(x2fit, 1) /= npts ) then
      write(unit=error_unit, fmt=*) &
&       "Error: invalid X2FIT size (NPTS=", npts, " != ", size(x2fit, 1), ")"
      data_ok = .false.
    end if
    if ( size(y2fit, 1) /= npts ) then
      write(unit=error_unit, fmt=*) &
&       "Error: invalid Y2FIT size (NPTS=", npts, " != ", size(y2fit, 1), ")"
      data_ok = .false.
    end if
    if ( associated(fit_data%coeffs) ) then
      write(unit=error_unit, fmt=*) &
&       "Error: Gaussian coefficients already allocated (size=", &
&       shape(fit_data%coeffs), ")"
      data_ok = .false.
    end if
    if ( associated(fit_data%x2fit) ) then
      write(unit=error_unit, fmt=*) &
&       "Error: internal X2FIT already allocated (size=", &
&       size(fit_data%x2fit, 1), ")"
      data_ok = .false.
    end if
    if ( associated(fit_data%y2fit) ) then
      write(unit=error_unit, fmt=*) &
&       "Error: internal Y2FIT already allocated (size=", &
&       size(fit_data%y2fit, 1), ")"
      data_ok = .false.
    end if
    if ( .not. data_ok ) return

    ! Set internal fields
    fit_data%ngfs = ngfs
    fit_data%npts = npts

    ! Allocate and initialise internal arrays
    ! Note: we copy the data to make sure that it will remain untouched
    allocate(fit_data%coeffs(2,ngfs))
    fit_data%coeffs(:,:) = 0.0_gfdp
    allocate(fit_data%x2fit(npts))
    fit_data%x2fit(:) = x2fit(:)
    allocate(fit_data%y2fit(npts))
    fit_data%y2fit(:) = y2fit(:)

  end subroutine gaufre_data_init

  ! ***************************************************************************
  ! *** Utility routines                                                    ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Converts a trial vector to a Gaussian expansion
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[inout] fit_data: data structure containing Gaussian fitting
  !!                         parameters
  !! \param[in] trial: vector containing the flattened coefficients
  ! ***************************************************************************
  subroutine gaufre_data_set_coeffs(fit_data, trial)

    implicit none

    ! Arguments
    type(gaufre_data_t), intent(in) :: fit_data
    real(gfdp), intent(in) :: trial(2*fit_data%ngfs)

    ! Local variables
    integer :: igfs

    ! -------------------------------------------------------------------------

    do igfs = 1,fit_data%ngfs
      fit_data%coeffs(1,igfs) = abs(trial(2*igfs-1))
      fit_data%coeffs(2,igfs) = trial(2*igfs)
    end do

  end subroutine gaufre_data_set_coeffs

  ! ***************************************************************************
  !> \brief Converts a Gaussian expansion to a trial vector
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[in] fit_data: data structure containing Gaussian fitting
  !!                      parameters
  !! \param[out] trial: vector containing the flattened coefficients
  ! ***************************************************************************
  subroutine gaufre_data_set_trial(fit_data, trial)

    implicit none

    ! Arguments
    type(gaufre_data_t), intent(in) :: fit_data
    real(gfdp), intent(out) :: trial(2*fit_data%ngfs)

    ! Local variables
    integer :: igfs

    ! -------------------------------------------------------------------------

    do igfs = 1,fit_data%ngfs
      trial(2*igfs-1) = fit_data%coeffs(1,igfs)
      trial(2*igfs) = fit_data%coeffs(2,igfs)
    end do

  end subroutine gaufre_data_set_trial

end module gaufre_data
