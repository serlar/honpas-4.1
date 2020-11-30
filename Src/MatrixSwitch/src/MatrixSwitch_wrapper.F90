#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==============================================================================!
!> @brief MatrixSwitch C bindings.
!==============================================================================!
module MatrixSwitch_wrapper
  use MatrixSwitch_ops, only: dp
  use MatrixSwitch_wrapper_params
  use MatrixSwitch, only: &
    mm_multiply_orig => mm_multiply, &
    m_add_orig => m_add, &
    m_trace_orig => m_trace, &
    mm_trace_orig => mm_trace, &
    m_scale_orig => m_scale, &
    m_set_orig => m_set, &
    m_set_element_orig => m_set_element, &
    m_get_element_orig => m_get_element, &
    m_register_sden_orig => m_register_sden, &
    m_allocate_orig => m_allocate, &
    m_deallocate_orig => m_deallocate, &
    m_copy_orig => m_copy, &
    m_convert_orig => m_convert
#ifdef HAVE_MPI
  use MatrixSwitch, only: &
    m_register_pdbc_orig => m_register_pdbc, &
#ifdef HAVE_SCALAPACK
    ms_scalapack_setup, &
#endif
    ms_lap_icontxt
#endif
#ifdef HAVE_PSPBLAS
  use MatrixSwitch, only: &
    m_copy_external_pdbcpcoo_orig => m_copy_external_pdbcpcoo, &
    m_copy_external_pdbcpcsc_orig => m_copy_external_pdbcpcsc, &
    m_register_pcoo_orig => m_register_pcoo, &
    m_register_pcsc_orig => m_register_pcsc
#endif

  implicit none

  private

  !**** INTERFACES ********************************!

  !============================================================================!
  !> @brief Wrapper to matrix-matrix multiplication.
  !============================================================================!
  interface mm_multiply
     module procedure mm_dmultiply
     module procedure mm_zmultiply
  end interface mm_multiply

  !============================================================================!
  !> @brief Wrapper to matrix addition.
  !============================================================================!
  interface m_add
     module procedure m_dadd
     module procedure m_zadd
  end interface m_add

  !============================================================================!
  !> @brief Wrapper to matrix trace.
  !============================================================================!
  interface m_trace
     module procedure m_dtrace
     module procedure m_ztrace
  end interface m_trace

  !============================================================================!
  !> @brief Wrapper to matrix product trace.
  !============================================================================!
  interface mm_trace
     module procedure mm_dtrace
     module procedure mm_ztrace
  end interface mm_trace

  !============================================================================!
  !> @brief Wrapper to scale matrix.
  !============================================================================!
  interface m_scale
     module procedure m_dscale
     module procedure m_zscale
  end interface m_scale

  !============================================================================!
  !> @brief Wrapper to set matrix.
  !============================================================================!
  interface m_set
     module procedure m_dset
     module procedure m_zset
  end interface m_set

  !============================================================================!
  !> @brief Wrapper to set matrix element.
  !============================================================================!
  interface m_set_element
     module procedure m_dset_element
     module procedure m_zset_element
  end interface m_set_element

  !============================================================================!
  !> @brief Wrapper to get matrix element.
  !============================================================================!
  interface m_get_element
     module procedure m_dget_element
     module procedure m_zget_element
  end interface m_get_element

  !============================================================================!
  !> @brief Wrapper to register matrix (simple dense, serial distribution).
  !============================================================================!
  interface m_register_sden
     module procedure m_register_sdden
     module procedure m_register_szden
  end interface m_register_sden

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Wrapper to register matrix (dense block cyclic, parallel
  !!        distribution).
  !============================================================================!
  interface m_register_pdbc
     module procedure m_register_pddbc
     module procedure m_register_pzdbc
  end interface m_register_pdbc
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (sparse coordinate list from dense
  !!        block cyclic, parallel distribution).
  !============================================================================!
  interface m_copy_external_pdbcpcoo
     module procedure m_copy_external_pddbcpdcoo
     module procedure m_copy_external_pzdbcpzcoo
  end interface m_copy_external_pdbcpcoo

  !============================================================================!
  !> @brief Wrapper to copy external matrix (compressed sparse column from
  !!        dense block cyclic, parallel distribution).
  !============================================================================!
  interface m_copy_external_pdbcpcsc
     module procedure m_copy_external_pddbcpdcsc
     module procedure m_copy_external_pzdbcpzcsc
  end interface m_copy_external_pdbcpcsc

  !============================================================================!
  !> @brief Wrapper to register matrix (sparse coordinate list, parallel
  !!        distribution).
  !============================================================================!
  interface m_register_pcoo
     module procedure m_register_pdcoo
     module procedure m_register_pzcoo
  end interface m_register_pcoo

  !============================================================================!
  !> @brief Wrapper to register matrix (compressed sparse column, parallel
  !!        distribution).
  !============================================================================!
  interface m_register_pcsc
     module procedure m_register_pdcsc
     module procedure m_register_pzcsc
  end interface m_register_pcsc
#endif

  !************************************************!

  public :: ms_wrapper_open
  public :: ms_wrapper_close
  public :: ms_is_initialized
  public :: ms_is_serial
  public :: ms_is_real
  public :: ms_is_square
  public :: ms_is_sparse
  public :: ms_dim1
  public :: ms_dim2
  public :: mm_multiply
  public :: m_add
  public :: m_trace
  public :: mm_trace
  public :: m_scale
  public :: m_set
  public :: m_set_element
  public :: m_get_element
  public :: m_register_sden
  public :: m_allocate
  public :: m_deallocate
  public :: m_copy
  public :: m_convert
#ifdef HAVE_MPI
  public :: m_register_pdbc
#ifdef HAVE_SCALAPACK
  public :: ms_scalapack_setup
#endif
  public :: ms_lap_icontxt
#endif
#ifdef HAVE_PSPBLAS
  public :: m_copy_external_pdbcpcoo
  public :: m_copy_external_pdbcpcsc
  public :: m_register_pcoo
  public :: m_register_pcsc
#endif

contains

  !============================================================================!
  !> @brief MatrixSwitch wrapper setup.
  !!
  !! Sets up the MatrixSwitch wrapper with a fixed number of matrices, each of
  !! which is identified by a unique key.
  !!
  !! @param[in] num_matrices The number of matrices to store in the wrapper.
  !! @param[in] keys         The array of keys for accessing the stored
  !!                         matrices. Each key cannot exceed \p max_key_length
  !!                         characters.
  !============================================================================!
  subroutine ms_wrapper_open(num_matrices,keys)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: num_matrices ! number of matrices to use
    character(*), intent(in) :: keys(num_matrices) ! mapping keys


    !**********************************************!

    allocate(ms_matrices(num_matrices))
    allocate(ms_keys(num_matrices))

    ms_num_matrices=num_matrices
    ms_keys=keys

  end subroutine ms_wrapper_open

  !============================================================================!
  !> @brief Close the MatrixSwitch wrapper.
  !!
  !! All MatrixSwitch matrices stored in the wrapper are deallocated, and the
  !! map is deallocated.
  !============================================================================!
  subroutine ms_wrapper_close()
    implicit none

    !**** INTERNAL ********************************!

    integer :: i

    !**********************************************!

    if (allocated(ms_keys)) deallocate(ms_keys)
    if (allocated(ms_matrices)) then
      do i=1,ms_num_matrices
        call m_deallocate_orig(ms_matrices(i))
      end do
      deallocate(ms_matrices)
    end if

  end subroutine ms_wrapper_close

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_initialized).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_initialized
  !============================================================================!
  logical function ms_is_initialized(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_initialized=ms_matrices(ms_lookup(m_name))%is_initialized

  end function ms_is_initialized

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_serial).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_serial
  !============================================================================!
  logical function ms_is_serial(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_serial=ms_matrices(ms_lookup(m_name))%is_serial

  end function ms_is_serial

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_real).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_real
  !============================================================================!
  logical function ms_is_real(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_real=ms_matrices(ms_lookup(m_name))%is_real

  end function ms_is_real

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_square).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_square
  !============================================================================!
  logical function ms_is_square(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_square=ms_matrices(ms_lookup(m_name))%is_square

  end function ms_is_square

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_sparse).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_sparse
  !============================================================================!
  logical function ms_is_sparse(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_sparse=ms_matrices(ms_lookup(m_name))%is_sparse

  end function ms_is_sparse

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_dim1).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_dim1
  !============================================================================!
  integer function ms_dim1(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_dim1=ms_matrices(ms_lookup(m_name))%dim1

  end function ms_dim1

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_dim2).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_dim2
  !============================================================================!
  integer function ms_dim2(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_dim2=ms_matrices(ms_lookup(m_name))%dim2

  end function ms_dim2

  !============================================================================!
  !> @brief Wrapper to allocate matrix.
  !============================================================================!
  subroutine m_allocate(m_name,i,j,label)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label
    character(*), intent(in) :: m_name

    integer, intent(in) :: i
    integer, intent(in) :: j

    !**********************************************!

    call m_allocate_orig(ms_matrices(ms_lookup(m_name)),i,j,label)

  end subroutine m_allocate

  !============================================================================!
  !> @brief Wrapper to deallocate matrix.
  !============================================================================!
  subroutine m_deallocate(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_deallocate_orig(ms_matrices(ms_lookup(m_name)))

  end subroutine m_deallocate

  !============================================================================!
  !> @brief Wrapper to copy matrix.
  !============================================================================!
  subroutine m_copy(m_name,A,label,threshold,threshold_is_soft)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label

    logical, intent(in), optional :: threshold_is_soft

    real(dp), intent(in), optional :: threshold

    character(*), intent(in) :: A
    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_orig(ms_matrices(ms_lookup(m_name)),ms_matrices(ms_lookup(A)),label,threshold,threshold_is_soft)

  end subroutine m_copy

  !============================================================================!
  !> @brief Wrapper to convert matrix format.
  !============================================================================!
  subroutine m_convert(m_name,label,threshold,threshold_is_soft)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label

    logical, intent(in), optional :: threshold_is_soft

    real(dp), intent(in), optional :: threshold

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_convert_orig(ms_matrices(ms_lookup(m_name)),label,threshold,threshold_is_soft)

  end subroutine m_convert

  !============================================================================!
  !> @brief Wrapper to matrix-matrix multiplication (real version).
  !============================================================================!
  subroutine mm_dmultiply(A,opA,B,opB,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(1), intent(in) :: opB
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    character(*), intent(in) :: A
    character(*), intent(in) :: B
    character(*), intent(in) :: C

    !**********************************************!

    call mm_multiply_orig(ms_matrices(ms_lookup(A)),opA,ms_matrices(ms_lookup(B)),opB,ms_matrices(ms_lookup(C)),alpha,beta,label)

  end subroutine mm_dmultiply

  !============================================================================!
  !> @brief Wrapper to matrix-matrix multiplication (complex version).
  !============================================================================!
  subroutine mm_zmultiply(A,opA,B,opB,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(1), intent(in) :: opB
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    character(*), intent(in) :: A
    character(*), intent(in) :: B
    character(*), intent(in) :: C

    !**********************************************!

    call mm_multiply_orig(ms_matrices(ms_lookup(A)),opA,ms_matrices(ms_lookup(B)),opB,ms_matrices(ms_lookup(C)),alpha,beta,label)

  end subroutine mm_zmultiply

  !============================================================================!
  !> @brief Wrapper to matrix addition (real version).
  !============================================================================!
  subroutine m_dadd(A,opA,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    character(*), intent(in) :: A
    character(*), intent(in) :: C

    !**********************************************!

    call m_add_orig(ms_matrices(ms_lookup(A)),opA,ms_matrices(ms_lookup(C)),alpha,beta,label)

  end subroutine m_dadd

  !============================================================================!
  !> @brief Wrapper to matrix addition (complex version).
  !============================================================================!
  subroutine m_zadd(A,opA,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    character(*), intent(in) :: A
    character(*), intent(in) :: C

    !**********************************************!

    call m_add_orig(ms_matrices(ms_lookup(A)),opA,ms_matrices(ms_lookup(C)),alpha,beta,label)

  end subroutine m_zadd

  !============================================================================!
  !> @brief Wrapper to matrix trace (real version).
  !============================================================================!
  subroutine m_dtrace(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    character(*), intent(in) :: A

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**********************************************!

    call m_trace_orig(ms_matrices(ms_lookup(A)),alpha,label)

  end subroutine m_dtrace

  !============================================================================!
  !> @brief Wrapper to matrix trace (complex version).
  !============================================================================!
  subroutine m_ztrace(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    character(*), intent(in) :: A

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**********************************************!

    call m_trace_orig(ms_matrices(ms_lookup(A)),alpha,label)

  end subroutine m_ztrace

  !============================================================================!
  !> @brief Wrapper to matrix product trace (real version).
  !============================================================================!
  subroutine mm_dtrace(A,B,alpha,label)
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    character(*), intent(in) :: A
    character(*), intent(in) :: B

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**********************************************!

    call mm_trace_orig(ms_matrices(ms_lookup(A)),ms_matrices(ms_lookup(B)),alpha,label)

  end subroutine mm_dtrace

  !============================================================================!
  !> @brief Wrapper to matrix product trace (complex version).
  !============================================================================!
  subroutine mm_ztrace(A,B,alpha,label)
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    character(*), intent(in) :: A
    character(*), intent(in) :: B

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**********************************************!

    call mm_trace_orig(ms_matrices(ms_lookup(A)),ms_matrices(ms_lookup(B)),alpha,label)

  end subroutine mm_ztrace

  !============================================================================!
  !> @brief Wrapper to scale matrix (real version).
  !============================================================================!
  subroutine m_dscale(C,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    real(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_scale_orig(ms_matrices(ms_lookup(C)),beta,label)

  end subroutine m_dscale

  !============================================================================!
  !> @brief Wrapper to scale matrix (complex version).
  !============================================================================!
  subroutine m_zscale(C,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_scale_orig(ms_matrices(ms_lookup(C)),beta,label)

  end subroutine m_zscale

  !============================================================================!
  !> @brief Wrapper to set matrix (real version).
  !============================================================================!
  subroutine m_dset(C,seC,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_set_orig(ms_matrices(ms_lookup(C)),seC,alpha,beta,label)

  end subroutine m_dset

  !============================================================================!
  !> @brief Wrapper to set matrix (complex version).
  !============================================================================!
  subroutine m_zset(C,seC,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_set_orig(ms_matrices(ms_lookup(C)),seC,alpha,beta,label)

  end subroutine m_zset

  !============================================================================!
  !> @brief Wrapper to set matrix element (real version).
  !============================================================================!
  subroutine m_dset_element(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_set_element_orig(ms_matrices(ms_lookup(C)),i,j,alpha,beta,label)

  end subroutine m_dset_element

  !============================================================================!
  !> @brief Wrapper to set matrix element (complex version).
  !============================================================================!
  subroutine m_zset_element(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_set_element_orig(ms_matrices(ms_lookup(C)),i,j,alpha,beta,label)

  end subroutine m_zset_element

  !============================================================================!
  !> @brief Wrapper to get matrix element (real version).
  !============================================================================!
  subroutine m_dget_element(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    character(*), intent(in) :: C

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**********************************************!

    call m_get_element_orig(ms_matrices(ms_lookup(C)),i,j,alpha,label)

  end subroutine m_dget_element

  !============================================================================!
  !> @brief Wrapper to get matrix element (complex version).
  !============================================================================!
  subroutine m_zget_element(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    character(*), intent(in) :: C

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**********************************************!

    call m_get_element_orig(ms_matrices(ms_lookup(C)),i,j,alpha,label)

  end subroutine m_zget_element

  !============================================================================!
  !> @brief Wrapper to register matrix (simple dense, serial distribution, real
  !!        version).
  !============================================================================!
  subroutine m_register_sdden(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_sden_orig(ms_matrices(ms_lookup(m_name)),A)

  end subroutine m_register_sdden

  !============================================================================!
  !> @brief Wrapper to register matrix (simple dense, serial distribution,
  !!        complex version).
  !============================================================================!
  subroutine m_register_szden(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    complex(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_sden_orig(ms_matrices(ms_lookup(m_name)),A)

  end subroutine m_register_szden

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Wrapper to register matrix (dense block cyclic, parallel
  !!        distribution, real version).
  !============================================================================!
  subroutine m_register_pddbc(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pdbc_orig(ms_matrices(ms_lookup(m_name)),A,desc)

  end subroutine m_register_pddbc
#endif

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Wrapper to register matrix (dense block cyclic, parallel
  !!        distribution, complex version).
  !============================================================================!
  subroutine m_register_pzdbc(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9)

    complex(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pdbc_orig(ms_matrices(ms_lookup(m_name)),A,desc)

  end subroutine m_register_pzdbc
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (sparse coordinate list from dense
  !!        block cyclic, parallel distribution, real version).
  !============================================================================!
  subroutine m_copy_external_pddbcpdcoo(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)
    real(dp), intent(in), optional :: threshold

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pdbcpcoo_orig(ms_matrices(ms_lookup(m_name)),A,desc,threshold)

  end subroutine m_copy_external_pddbcpdcoo
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (sparse coordinate list from dense
  !!        block cyclic, parallel distribution, complex version).
  !============================================================================!
  subroutine m_copy_external_pzdbcpzcoo(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), optional :: threshold

    complex(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pdbcpcoo_orig(ms_matrices(ms_lookup(m_name)),A,desc,threshold)

  end subroutine m_copy_external_pzdbcpzcoo
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (compressed sparse column from
  !!        dense block cyclic, parallel distribution, real version).
  !============================================================================!
  subroutine m_copy_external_pddbcpdcsc(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)
    real(dp), intent(in), optional :: threshold

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pdbcpcsc_orig(ms_matrices(ms_lookup(m_name)),A,desc,threshold)

  end subroutine m_copy_external_pddbcpdcsc
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (compressed sparse column from
  !!        dense block cyclic, parallel distribution, complex version).
  !============================================================================!
  subroutine m_copy_external_pzdbcpzcsc(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), optional :: threshold

    complex(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pdbcpcsc_orig(ms_matrices(ms_lookup(m_name)),A,desc,threshold)

  end subroutine m_copy_external_pzdbcpzcsc
#endif

#ifdef HAVE_PSPBLAS

  !============================================================================!
  !> @brief Wrapper to register matrix (sparse coordinate list, parallel
  !!        distribution, real version).
  !============================================================================!
  subroutine m_register_pdcoo(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    real(dp), intent(in), target :: val(:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pcoo_orig(ms_matrices(ms_lookup(m_name)),idx1,idx2,val,desc)

  end subroutine m_register_pdcoo
#endif

#ifdef HAVE_PSPBLAS

  !============================================================================!
  !> @brief Wrapper to register matrix (sparse coordinate list, parallel
  !!        distribution, complex version).
  !============================================================================!
  subroutine m_register_pzcoo(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    complex(dp), intent(in), target :: val(:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pcoo_orig(ms_matrices(ms_lookup(m_name)),idx1,idx2,val,desc)

  end subroutine m_register_pzcoo
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to register matrix (compressed sparse column, parallel
  !!        distribution, real version).
  !============================================================================!
  subroutine m_register_pdcsc(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    real(dp), intent(in), target :: val(:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pcsc_orig(ms_matrices(ms_lookup(m_name)),idx1,idx2,val,desc)

  end subroutine m_register_pdcsc
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to register matrix (compressed sparse column, parallel
  !!        distribution, complex version).
  !============================================================================!
  subroutine m_register_pzcsc(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    complex(dp), intent(in), target :: val(:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pcsc_orig(ms_matrices(ms_lookup(m_name)),idx1,idx2,val,desc)

  end subroutine m_register_pzcsc
#endif

end module MatrixSwitch_wrapper
