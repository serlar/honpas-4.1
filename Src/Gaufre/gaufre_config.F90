! *** Module: gaufre_config ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Global configuration for the fitting of Natural Atomic Orbitals
!!        to Gaussians
!!
!! This module defines an algorithm-independent data structure to configure
!! the fitting of Natural Atomic Orbitals by Gaussians.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 04.2018 Created [Yann Pouillon]
! *****************************************************************************
module gaufre_config

  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

  use gaufre_common, only: gfdp, GAUFRE_EPSILON_DEFAULT

  implicit none

  private

  !> Data structure to store the Gaussian fitting configuration
  !!
  !! This data type defines the criteria and restrictions to apply
  !! when fitting (X,Y) data with a sum of Gaussians.
  !!
  !! Some algorithms (e.g. CMPFIT) support the definition of boundaries
  !! for the fitting coefficients. These can be defined by setting some
  !! elements of the \ref gaufre_conf_t::limits field to .true. and
  !! setting the corresponding elements of the \ref gaufre_conf_t::bounds
  !! field to the desired limits. The indices of these two arrays are
  !! defined as follows:
  !!   - (1,1): lower bound of the exponential exponents;
  !!   - (2,1): upper bound of the exponential exponents;
  !!   - (1,2): lower bound of the scaling factors;
  !!   - (2,2): upper bound of the scaling factors.
  !!
  !! \note Depending on the algorithm used to build the Gaussian expansion,
  !! some fields will or will not be taken into account. The relevant
  !! parameters for each algorithm are described in the corresponding
  !! documentation.
  !!
  type, public :: gaufre_conf_t
    real(gfdp) :: tolerance = 0.0_gfdp    !< Tolerance on the residue of
                                          !! the Gaussian expansion (main
                                          !! criterion for the fitting)
    integer  :: maxiter = 1000            !< Maximum number of algorithm
                                          !! iterations
    integer  :: nprint = 1                !< How often to print information
    logical  :: limits(2,2) = .false.     !< Whether to limit the values of
                                          !! the Gaussian coefficients (see
                                          !! the bounds array)
    real(gfdp) :: bounds(2,2) = 0.0_gfdp  !< Bounding intervals for the
                                          !! Gaussian coefficients
    real(gfdp) :: epsilon = &
&     GAUFRE_EPSILON_DEFAULT              !< Internal precision tolerance,
                                          !! for algorithms that support it
  end type

  ! Public resources
  public :: &
&   gaufre_conf_dump, &
&   gaufre_conf_free, &
&   gaufre_conf_init

contains

  ! ***************************************************************************
  !> \brief Pretty-printer for the gaufre_conf_t data structure
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[in] fit_conf: data structure containing fitting Gaussian parameters
  ! ***************************************************************************
  subroutine gaufre_conf_dump(fit_conf)

    implicit none

    ! Arguments
    type(gaufre_conf_t), intent(in) :: fit_conf

    ! Local variables
    logical :: data_ok

    ! -------------------------------------------------------------------------

    ! Header
    write(unit=output_unit, fmt='(A)') "GAUFRE CONFIGURATION BEGIN"
    write(unit=output_unit, fmt='(A)')

    ! Show scalar components
    write(unit=output_unit, fmt='(2X,"* ",A16,1X,"=",1X,I10)') &
&     "MAXITER", fit_conf%maxiter
    write(unit=output_unit, fmt='(2X,"* ",A16,1X,"=",1X,I10)') &
&     "NPRINT", fit_conf%nprint
    write(unit=output_unit, fmt='(2X,"* ",A16,1X,"=",4(1X,A2,":",1X,L10))') &
&     "LIMITS", "XL", fit_conf%limits(1,1), "XU", fit_conf%limits(2,1), &
&     "SL", fit_conf%limits(1,2), "SU", fit_conf%limits(2,2)
    write(unit=output_unit, fmt='(2X,"* ",A16,1X,"=",4(1X,A2,":",1X,E10.2))') &
&     "LIMITS", "XL", fit_conf%bounds(1,1), "XU", fit_conf%bounds(2,1), &
&     "SL", fit_conf%bounds(1,2), "SU", fit_conf%bounds(2,2)
    write(unit=output_unit, fmt='(2X,"* ",A16,1X,"=",1X,E10.2)') &
&     "TOLERANCE", fit_conf%tolerance

    ! Show parameter consistency
    write(unit=output_unit, fmt='(A)')
    write(unit=output_unit, fmt='(2X,A)') "CONSISTENCY CHECKS BEGIN"
    write(unit=output_unit, fmt='(A)')

    ! Check for improper input parameters
    data_ok = .true.
    write(unit=output_unit, fmt=*)
    write(unit=output_unit, fmt=*) "CONSISTENCY CHECKS BEGIN"
    if ( fit_conf%maxiter <= 0 ) then
      write(unit=output_unit, fmt='(4X,"* ",A)') &
&       "ERROR: invalid parameter (MAXITER=", fit_conf%maxiter, ")"
      data_ok = .false.
    end if
    if ( data_ok ) then
      write(unit=output_unit, fmt='(4X,"* ",A)') "All parameters OK"
    end if
    write(unit=output_unit, fmt='(A)')
    write(unit=output_unit, fmt='(2X,A)') "CONSISTENCY CHECKS END"

    ! Footer
    write(unit=output_unit, fmt='(A)')
    write(unit=output_unit, fmt='(A)') "GAUFRE CONFIGURATION END"

  end subroutine gaufre_conf_dump

  ! ***************************************************************************
  !> \brief Destructor for the Gaussian expansion configuration
  !!
  !! This routine frees all memory allocated for the Gaussian fitting procedure
  !! configuration and should be called once the data structure is not needed
  !! anymore.
  !!
  !! \note This routine is idempotent and can be run as many times as desired.
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[in] fit_conf: data structure containing a Gaussian fitting
  !!                      configuration
  ! ***************************************************************************
  subroutine gaufre_conf_free(fit_conf)

    implicit none

    type(gaufre_conf_t), intent(inout) :: fit_conf

    fit_conf%maxiter = 0
    fit_conf%nprint = 0
    fit_conf%limits(:,:) = .false.
    fit_conf%bounds(:,:) = 0.0_gfdp
    fit_conf%epsilon = GAUFRE_EPSILON_DEFAULT
    fit_conf%tolerance = 0.0_gfdp

  end subroutine gaufre_conf_free

  ! ***************************************************************************
  !> \brief Constructor for the Gaussian expansion configuration
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
  !! \param[out] fit_conf: data structure containing fitting Gaussian parameters
  !! \param[in] tolerance: tolerance on the fitting residue
  !! \param[in] maxiter: maximum number of algorithm iterations (optional)
  !! \param[in] nprint: how often to print information, for implementations
  !!                    supporting it (optional)
  !! \param[in] limits: whether to put boundaries to the expansion
  !!                    coefficients (optional, see \ref gaufre_conf_t::limits)
  !! \param[in] bounds: bounding values for the expansion coefficients
  !!                    (optional, see \ref gaufre_conf_t::bounds)
  ! ***************************************************************************
  subroutine gaufre_conf_init(fit_conf, tolerance, maxiter, nprint, &
&   limits, bounds)

    implicit none

    ! Arguments
    type(gaufre_conf_t), intent(out) :: fit_conf
    real(gfdp), intent(in) :: tolerance
    integer, optional, intent(in) :: maxiter
    integer, optional, intent(in) :: nprint
    logical, optional, intent(in) :: limits(2,2)
    real(gfdp), optional, intent(in) :: bounds(2,2)

    ! -------------------------------------------------------------------------

    ! Import problem dimensions and data to fit
    fit_conf%tolerance = tolerance
    if ( present(maxiter) ) fit_conf%maxiter = maxiter
    if ( present(nprint) ) fit_conf%nprint = nprint
    if ( present(limits) ) fit_conf%limits(:,:) = limits(:,:)
    if ( present(bounds) ) fit_conf%bounds(:,:) = bounds(:,:)

  end subroutine gaufre_conf_init

end module gaufre_config
