! *** Module: gaufre_common ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Global constants for the fitting of Natural Atomic Orbitals
!!        to Gaussians
!!
!! This module contains global constants and basic common data for
!! the fitting of Natural Atomic Orbitals by Gaussians.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 04.2018 Created [Yann Pouillon]
! *****************************************************************************
module gaufre_common

  implicit none

  private

                    ! ------------------------------------ !

  !
  ! Constants
  !

  !> Real precision selector
  integer, parameter, public :: gfdp = kind(1.0d0)

  !> Default nternal precision for algorithms that support it
  real(gfdp), parameter, public :: GAUFRE_EPSILON_DEFAULT = 2.2204460e-16

                    ! ------------------------------------ !

  !
  ! Data formats
  !

  !> Unknown data format, to detect uninitialised variables
  integer, parameter, public :: GAUFRE_DATA_FMT_UNKNOWN = -1

  !> Column-wise representation of (X,Y) data
  integer, parameter, public :: GAUFRE_DATA_FMT_COLUMNS = 1

  !> Flexible Data Format representation of data objects
  integer, parameter, public :: GAUFRE_DATA_FMT_FDF     = 2

  !> YAML representation of data objects
  integer, parameter, public :: GAUFRE_DATA_FMT_YAML    = 3

                    ! ------------------------------------ !

  !
  ! I/O operations
  !

  !> Unknown I/O mode, to detect uninitialised variables
  integer, parameter, public :: GAUFRE_IO_UNKNOWN = -1

  !> I/O mode to load data from a file
  integer, parameter, public :: GAUFRE_IO_LOAD = 1

  !> I/O mode to save data to a file
  integer, parameter, public :: GAUFRE_IO_SAVE = 2

  !> I/O file format for SIESTA
  integer, parameter, public :: GAUFRE_IO_FMT_SIESTA = 1

                    ! ------------------------------------ !

  !
  ! Gaussian fitting implementations
  !

  !> Unknown implementation, to detect uninitialised variables
  integer, parameter, public :: GAUFRE_LM_UNKNOWN = -1

  !> CMPFIT implementation of the Levenberg-Marquardt algorithm,
  !! available at https://www.physics.wisc.edu/~craigm/idl/cmpfit.html
  integer, parameter, public :: GAUFRE_LM_CMPFIT  = 1

  !> Levmar implementation of the Levenberg-Marquardt algorithm,
  !! available at http://users.ics.forth.gr/~lourakis/levmar/
  integer, parameter, public :: GAUFRE_LM_LEVMAR  = 2

  !> MINPACK implementation of the Levenberg-Marquardt algorithm,
  !! available at:
  !!   - http://www.netlib.org/minpack/
  !!   - http://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html
  !!   - https://github.com/certik/minpack
  !!
  integer, parameter, public :: GAUFRE_LM_MINPACK = 3

  !> NL2SOL implementation of the Levenberg-Marquardt algorithm,
  !! available at http://people.sc.fsu.edu/~jburkardt/f_src/nl2sol/nl2sol.html
  integer, parameter, public :: GAUFRE_LM_NL2SOL  = 4

  !> C/C++ MINPACK implementation of the Levenberg-Marquardt algorithm,
  !! available at http://devernay.free.fr/hacks/cminpack/
  integer, parameter, public :: GAUFRE_LM_CMINPACK  = 5

                    ! ------------------------------------ !

  ! Public resources
  public :: &
&   gaufre_has_method, &
&   gaufre_method_name

contains

  ! ***************************************************************************
  !> \brief Returns whether a method is implemented in GAUFRE
  !!
  !! \par History
  !!      - 05.2018 Created [Yann Pouillon]
  !!
  !! \param[in] method_index: index of the method (see \ref gaufre_common)
  ! ***************************************************************************
  function gaufre_has_method(method_index)

    implicit none

    ! Arguments
    integer, intent(in) :: method_index

    ! Local variables
    logical :: gaufre_has_method

    ! -------------------------------------------------------------------------

    ! We assume that a method is not implemented if not mentioned otherwise
    gaufre_has_method = .false.

    ! Explicitly mark actually implemented methods
    select case(method_index)
      case(GAUFRE_LM_MINPACK)
        gaufre_has_method = .true.
    end select

  end function gaufre_has_method

  ! ***************************************************************************
  !> \brief Returns a string representing the specified method index
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[in] method_index: index of the method (see \ref gaufre_common)
  ! ***************************************************************************
  function gaufre_method_name(method_index)

    implicit none

    ! Arguments
    integer, intent(in) :: method_index

    ! Local variables
    character(len=12) :: gaufre_method_name

    ! -------------------------------------------------------------------------

    select case(method_index)
      case(GAUFRE_LM_UNKNOWN)
        gaufre_method_name = "UNKNOWN     "
      case(GAUFRE_LM_CMPFIT)
        gaufre_method_name = "CMPFIT      "
      case(GAUFRE_LM_LEVMAR)
        gaufre_method_name = "LEVMAR      "
      case(GAUFRE_LM_MINPACK)
        gaufre_method_name = "MINPACK     "
      case(GAUFRE_LM_NL2SOL)
        gaufre_method_name = "NL2SOL      "
      case(GAUFRE_LM_CMINPACK)
        gaufre_method_name = "CMINPACK    "
      case default
        gaufre_method_name = "###ERR###   "
    end select

  end function gaufre_method_name

end module gaufre_common
