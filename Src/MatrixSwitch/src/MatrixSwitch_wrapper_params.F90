#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==============================================================================!
!> @brief Storage and auxiliary operations for MatrixSwitch C bindings.
!==============================================================================!
module MatrixSwitch_wrapper_params
  use MatrixSwitch, only: matrix
  use MatrixSwitch_ops, only : dp, die

  implicit none

  private

  !**** PARAMS ************************************!

  integer, parameter :: max_key_length=10 !< Maximum character length of keys.

  !**** VARIABLES *********************************!

  character(max_key_length), allocatable :: ms_keys(:) !< Array of keys.

  integer :: ms_num_matrices !< Number of matrices stored in the wrapper.

  type(matrix), allocatable :: ms_matrices(:) !< Array of MatrixSwitch matrices stored in the wrapper.

  !************************************************!

  public :: matrix
  public :: ms_keys
  public :: ms_num_matrices
  public :: ms_matrices
  public :: ms_lookup

contains

  !============================================================================!
  !> @brief MatrixSwitch wrapper map lookup.
  !!
  !! Finds the value bound to the given key, corresponding to the index of the
  !! matrix in the \p ms_matrices array.
  !!
  !! @param[in] key They key of the matrix.
  !! @return        The index of the matrix in \p ms_matrices.
  !============================================================================!
  integer function ms_lookup(key)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: key ! matrix key

    !**** INTERNAL ********************************!

    integer :: i

    !**********************************************!

    ms_lookup=0

    do i=1,ms_num_matrices
      if (key .eq. ms_keys(i)) then
        ms_lookup=i
        exit
      end if
    end do

    if (ms_lookup==0) call die('ms_lookup: key not found')

  end function ms_lookup

end module MatrixSwitch_wrapper_params
