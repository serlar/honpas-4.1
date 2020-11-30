#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==============================================================================!
!> @brief Storage and auxiliary operations for MatrixSwitch.
!==============================================================================!
module MatrixSwitch_ops
#ifdef HAVE_PSPBLAS
  use pspBLAS
#endif

  implicit none

  public

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300) !< Double precision.

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp) !< Complex 1.
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp) !< Complex i.
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp) !< Complex 0.

  !**** VARIABLES *********************************!

#ifdef HAVE_MPI
  character(1), save :: ms_lap_order !< Ordering of the BLACS process grid.

  integer, save :: ms_mpi_comm !< MPI communicator.
  integer, save :: ms_mpi_size !< Number of MPI processes.
  integer, save :: ms_mpi_rank !< MPI process rank.
  integer, save :: ms_lap_nprow !< Number of rows in the BLACS process grid.
  integer, save :: ms_lap_npcol !< Number of columns in the BLACS process grid.
  integer, save :: ms_lap_bs_def !< Default block size.
  integer, save :: ms_lap_bs_num !< Number of block size exceptions.
  integer, save :: ms_lap_icontxt !< BLACS context handle used by MatrixSwitch.
  integer, allocatable, save :: ms_lap_bs_list(:,:) !< Block size exception list.
#endif

  !**** TYPES *************************************!

  !============================================================================!
  !> @brief MatrixSwitch matrix type.
  !!
  !! This is the derived type that encapsulates all matrix storage
  !! possibilities and hides the details from the user. Typically, the elements
  !! below will never need to be accessed directly.
  !============================================================================!
  type matrix
     character(3) :: str_type !< Label identifying the storage format.

     logical :: is_initialized=.false. !< Has the matrix been initialized?
     logical :: is_serial !< Is the matrix serial or parallel distributed?
     logical :: is_real !< Is the matrix real or complex (both kind \p dp)?
     logical :: is_square !< Is the matrix square?
     logical :: is_sparse !< Is the matrix sparse?
     logical :: iaux1_is_allocated=.false. !< Is iaux1 directly allocated or it is a pointer?
     logical :: iaux2_is_allocated=.false. !< Is iaux2 directly allocated or it is a pointer?
     logical :: iaux3_is_allocated=.false. !< Is iaux3 directly allocated or it is a pointer?
     logical :: iaux4_is_allocated=.false. !< Is iaux4 directly allocated or it is a pointer?
     logical :: dval_is_allocated=.false. !< Is dval directly allocated or it is a pointer?
     logical :: zval_is_allocated=.false. !< Is zval directly allocated or it is a pointer?

     integer :: dim1 !< Row dimension size of the matrix.
     integer :: dim2 !< Column dimension size of the matrix.
     integer, pointer :: iaux1(:) => null() !< Auxiliary information for certain storage formats.
     integer, pointer :: iaux2(:) => null() !< Auxiliary information for certain storage formats.
     integer, pointer :: iaux3(:) => null() !< Auxiliary information for certain storage formats.
     integer, pointer :: iaux4(:) => null() !< Auxiliary information for certain storage formats.

     real(dp), pointer :: dval(:,:) => null() !< Matrix elements for a real matrix.

     complex(dp), pointer :: zval(:,:) => null() !< Matrix elements for a complex matrix.

#ifdef HAVE_PSPBLAS
     type(psp_matrix_spm) :: spm !< pspBLAS matrix type.
#endif
  end type matrix

  !**** INTERFACES ********************************!

  !============================================================================!
  !> @brief \p opM parameter converter.
  !!
  !! Converts the input parameters \p opA and \p opB in the subroutines
  !! \a mm_multiply and \a m_add from a character to a logical (real version)
  !! or integer (complex version) for internal use:
  !! \arg \c n / \c N mapped to \c .false. (real version) or \c 0 (complex
  !! version)
  !! \arg \c t / \c T mapped to \c .true. (real version) or \c 2 (complex
  !! version)
  !! \arg \c c / \c C mapped to \c .true. (real version) or \c 1 (complex
  !! version)
  !!
  !! @param[in]  seM Parameter 
  !! @param[out] luM 
  !============================================================================!
  interface process_opM
     module procedure process_lopM
     module procedure process_iopM
  end interface process_opM

  !************************************************!

contains

  !============================================================================!
  !> @brief \p opM parameter converter (real version).
  !============================================================================!
  subroutine process_lopM(opM,trM)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opM

    !**** INOUT ***********************************!

    logical, intent(inout) :: trM

    !**********************************************!

    if ((opM .eq. 'T') .or. &
         (opM .eq. 't') .or. &
         (opM .eq. 'C') .or. &
         (opM .eq. 'c')) then
       trM=.true.
    else if ((opM .eq. 'N') .or. &
         (opM .eq. 'n')) then
       trM=.false.
    else
       call die('process_lopM: invalid opM')
    end if

  end subroutine process_lopM

  !============================================================================!
  !> @brief \p opM parameter converter (complex version).
  !============================================================================!
  subroutine process_iopM(opM,tcM)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opM

    !**** INOUT ***********************************!

    integer, intent(inout) :: tcM

    !**********************************************!

    if ((opM .eq. 'T') .or. &
         (opM .eq. 't')) then
       tcM=2
    else if ((opM .eq. 'C') .or. &
         (opM .eq. 'c')) then
       tcM=1
    else if ((opM .eq. 'N') .or. &
         (opM .eq. 'n')) then
       tcM=0
    else
       call die('process_iopM: invalid opM')
    end if

  end subroutine process_iopM

  !============================================================================!
  !> @brief \p seC parameter converter.
  !!
  !! Converts the input parameter \p seC in the subroutine \a m_set from a
  !! character to an integer for internal use:
  !! \arg \c l / \c L mapped to \c 2
  !! \arg \c u / \c U mapped to \c 1
  !! \arg other: mapped to \c 0
  !!
  !! @param[in]  seM Parameter 
  !! @param[out] luM 
  !============================================================================!
  subroutine process_seM(seM,luM)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seM

    !**** INOUT ***********************************!

    integer, intent(inout) :: luM

    !**********************************************!

    if ((seM .eq. 'L') .or. &
         (seM .eq. 'l')) then
       luM=2
    else if ((seM .eq. 'U') .or. &
         (seM .eq. 'u')) then
       luM=1
    else
       luM=0
    end if

  end subroutine process_seM

  !============================================================================!
  !> @brief Code termination from error.
  !!
  !! Outputs a custom error message to file and then stops execution with a
  !! non-zero exit status.
  !!
  !! @param[in] message Error mesage to output before stopping.
  !============================================================================!
  subroutine die(message)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in), optional :: message

    !**** INTERNAL ********************************!

    integer :: err_unit

    !**********************************************!

    call io_assign(err_unit)
    open(unit=err_unit,file='MatrixSwitch.err',status='replace')
    write(err_unit,'(a)') 'FATAL ERROR in matrix_switch!'
    if (present(message)) write(err_unit,'(a)') message
#ifdef HAVE_MPI
    write(err_unit,'(a,1x,i5)') 'MPI rank:', ms_mpi_rank
#endif
    close(err_unit)
    stop 1

  end subroutine die

  subroutine io_assign(lun)
      integer, intent(out) :: lun
      logical used
      integer iostat

!     Looks for a free unit and assigns it to lun

      do lun= 10, 99
            inquire(unit=lun, opened=used, iostat=iostat)
            if (iostat .ne. 0) used = .true.
            if (.not. used) return
      enddo
  end subroutine io_assign

end module MatrixSwitch_ops
