! *** Module: gaufre_driver ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Fitting of Natural Atomic Orbitals to Gaussians
!!
!! This module provides a unified high-level interface that can be used
!! to determine the coefficients of zero-centered Gaussians which best
!! represent Natural Atomic Orbitals. Different algorithms are available.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 11.2017 Created [Yann Pouillon]
!!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
!!      - 04.2018 Enhanced fitting capabilities [Yann Pouillon]
! *****************************************************************************
module gaufre_driver

  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

  use gaufre_common
  use gaufre_config
  use gaufre_data

  implicit none

  private

  !> Data structure to drive the fitting of Natural Atomic Orbitals
  !!
  !! This structure stores all the required information to fit radial data
  !! and track the status of the procedure.
  type, public :: gaufre_driver_t
    integer :: fit_method
    type(gaufre_conf_t) :: fit_conf
    type(gaufre_data_t) :: fit_data
  end type gaufre_driver_t

  ! Public resources
  public :: &
&   gaufre_driver_find, &
&   gaufre_driver_free, &
&   gaufre_driver_init

contains

  ! ***************************************************************************
  ! *** Public routines                                                     ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Finds a Gaussian expansion from numerical data
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[inout] driver: initialised fitting data structure
  !! \param[in] method: method (implementation) to use to fit the data
  ! ***************************************************************************
  subroutine gaufre_driver_find(driver)

    use gaufre_minpack

    implicit none

    ! Arguments
    type(gaufre_driver_t), intent(inout) :: driver

    ! Local variables
    type(gaufre_minpack_t) :: fit_mpck

    ! -------------------------------------------------------------------------

    ! Select the method
    select case(driver%fit_method)

      case(GAUFRE_LM_CMPFIT)
        write(error_unit, fmt=*) "Fitting method not implemented: CMPFIT"
        return

      case(GAUFRE_LM_LEVMAR)
        write(error_unit, fmt=*) "Fitting method not implemented: LEVMAR"
        return

      case(GAUFRE_LM_MINPACK)
        call gaufre_minpack_wrap(driver%fit_data, &
&         driver%fit_conf%tolerance)

      case(GAUFRE_LM_NL2SOL)
        write(error_unit, fmt=*) "Fitting method not implemented: LEVMAR"
        return

      case default
        write(error_unit, fmt=*) "Invalid fitting method: ", driver%fit_method
        return

    end select

  end subroutine gaufre_driver_find

  ! ***************************************************************************
  !> \brief Wipes clean a fitting data structure
  !!
  !! \note This routine is idempotent and can be called as many times as
  !!       desired.
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[inout] driver: fitting data structure to clean
  ! ***************************************************************************
  subroutine gaufre_driver_free(driver)

    implicit none

    ! Arguments
    type(gaufre_driver_t), intent(inout) :: driver

    ! -------------------------------------------------------------------------

    ! Reset scalar variables
    driver%fit_method = GAUFRE_LM_UNKNOWN

    ! Delegate the resetting of internal structures
    call gaufre_conf_free(driver%fit_conf)
    call gaufre_data_free(driver%fit_data)

  end subroutine gaufre_driver_free

  ! ***************************************************************************
  !> \brief Initialises a fitting data structure
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[inout] driver: fitting data structure to initialise
  ! ***************************************************************************
  subroutine gaufre_driver_init(driver, method, ngfs, npts, x2fit, y2fit)

    implicit none

    ! Arguments
    type(gaufre_driver_t), intent(inout) :: driver
    integer, intent(in) :: method
    integer, intent(in) :: ngfs
    integer, intent(in) :: npts
    real(gfdp), intent(in) :: x2fit(npts)
    real(gfdp), intent(in) :: y2fit(npts)

    ! Local variables
    logical :: data_ok

    ! -------------------------------------------------------------------------

    ! Set internal fields
    driver%fit_method = method

    ! FIXME: hard-coded configuration
    call gaufre_conf_init(driver%fit_conf, tolerance=1.0e-4_gfdp)

    ! Delegate data initialisation
    call gaufre_data_init(driver%fit_data, ngfs, npts, x2fit, y2fit)

  end subroutine gaufre_driver_init

end module gaufre_driver
