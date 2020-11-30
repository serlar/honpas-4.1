! *** Module: libint1f ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Fortran interface to the Libint library (version 1.1.x)
!!
!! This module uses the libint1_glue library to access the resources provided
!! by Libint 1.1.x. In order to avoid memory management issues, quantities
!! are directly extracted from Libint whenever possible, in particular array
!! sizes.
!!
!! For completeness, here are the sizes of the core elements of the data
!! structures of Libint:
!!     - build_eri : LIBINT_MAX_AM - 1
!!     - Libint_t.vrr_classes : 2 * (LIBINT_MAX_AM - 1) + 1
!!     - prim_data.f : 4 * (LIBINT_MAX_AM - 1) + 1
!!
!! The LIBINT_MAX_AM parameter is a C preprocessing macro defined when
!! configuring the build of Libint. It corresponds to the value passed
!! to the configure script plus one. Its default value is 5. In order
!! to make it available from the Fortran side, the libint1_glue library
!! provides the libint1_get_libint_max_am function.
!!
!! The whole structure of Libint revolves around build_eri, a 4-dimensional
!! array of procedure pointers, which makes it relatively complex to safely
!! call the routines of Libint from Fortran programs. This is why the current
!! Libint interface fully encapsulates the internal structure of Libint_t as
!! well as build_eri.
!!
!! \note
!!     The F field of the \ref libint1f::libint1f_prim_t type is not
!!     deallocated automatically. You have to do it yourself for now.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 05.2018 Created [Yann Pouillon]
! *****************************************************************************
module libint1f

  use, intrinsic :: iso_c_binding

  implicit none

  private

  public :: &
&   libint1f_build_eris, &
&   libint1f_free, &
&   libint1f_init
!&   libint1f_dump_data, &

  ! Private internal constants
  integer, parameter :: dp = kind(0.0d0)   !< Default real precision

  !> \brief Fortran wrapper for the data structures of Libint
  !!
  !! This type hides away the internals of the data structures of Libint
  !! from Fortran programs. This is the recommended way to handle C data
  !! structures within the ESL (https://esl.cecam.org/).
  type :: libint1f_t
    type(c_ptr), private :: ptr = C_NULL_PTR
  end type libint1f_t

  !> \brief Fortran type to store a prim_data structure from Libint
  !!
  !! This is a convenience type to make it more comfortable to access
  !! a Libint prim_data structure from a Fortran program. We do this
  !! because this is the part of Libint_t that is most accessed when
  !! calling Libint from Fortran.
  type :: libint1f_prim_t
    real(dp), pointer :: f(:) => null()
    real(dp) :: oo2n = 0.0_dp
    real(dp) :: oo2p = 0.0_dp
    real(dp) :: oo2z = 0.0_dp
    real(dp) :: oo2zn = 0.0_dp
    real(dp) :: pon = 0.0_dp
    real(dp) :: poz = 0.0_dp
    real(dp) :: ss_r12_ss = 0.0_dp
    real(dp) :: twozeta_a = 0.0_dp
    real(dp) :: twozeta_b = 0.0_dp
    real(dp) :: twozeta_c = 0.0_dp
    real(dp) :: twozeta_d = 0.0_dp
    real(dp) :: u(3,6) = 0.0_dp
  end type libint1f_prim_t

  interface

    ! ************************************************************************
    ! Global getters
    ! ************************************************************************

    function libint1_get_libint_max_am() &
&                bind(c)
      import
      integer(kind=c_int) :: libint1_get_libint_max_am
    end function libint1_get_libint_max_am

    function libint1_get_libint_opt_am() &
&                bind(c)
      import
      integer(kind=c_int) :: libint1_get_libint_opt_am
    end function libint1_get_libint_opt_am

    ! ************************************************************************
    ! Getters for Libint_t
    ! ************************************************************************

    function libint1_get_libt_ab(ptr, itbl) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      integer(kind=c_int), value :: itbl
      real(kind=c_double) :: libint1_get_libt_ab
    end function libint1_get_libt_ab

    function libint1_get_libt_cd(ptr, itbl) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      integer(kind=c_int), value :: itbl
      real(kind=c_double) :: libint1_get_libt_cd
    end function libint1_get_libt_cd

    ! ************************************************************************
    ! Getters for prim_data
    ! ************************************************************************

    function libint1_get_prim_f(ptr, itbl) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      integer(kind=c_int), value :: itbl
      real(kind=c_double) :: libint1_get_prim_f
    end function libint1_get_prim_f

    function libint1_get_prim_oo2n(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_oo2n
    end function libint1_get_prim_oo2n

    function libint1_get_prim_oo2p(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_oo2p
    end function libint1_get_prim_oo2p

    function libint1_get_prim_oo2z(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_oo2z
    end function libint1_get_prim_oo2z

    function libint1_get_prim_oo2zn(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_oo2zn
    end function libint1_get_prim_oo2zn

    function libint1_get_prim_pon(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_pon
    end function libint1_get_prim_pon

    function libint1_get_prim_poz(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_poz
    end function libint1_get_prim_poz

    function libint1_get_prim_ss_r12_ss(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_ss_r12_ss
    end function libint1_get_prim_ss_r12_ss

    function libint1_get_prim_twozeta_a(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_twozeta_a
    end function libint1_get_prim_twozeta_a

    function libint1_get_prim_twozeta_b(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_twozeta_b
    end function libint1_get_prim_twozeta_b

    function libint1_get_prim_twozeta_c(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_twozeta_c
    end function libint1_get_prim_twozeta_c

    function libint1_get_prim_twozeta_d(ptr) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      real(kind=c_double) :: libint1_get_prim_twozeta_d
    end function libint1_get_prim_twozeta_d

    function libint1_get_prim_u(ptr, itbl, jtbl) &
&                bind(c)
      import
      type(c_ptr), value :: ptr
      integer(c_int), value :: itbl
      integer(c_int), value :: jtbl
      real(kind=c_double) :: libint1_get_prim_u
    end function libint1_get_prim_u

    ! ************************************************************************
    ! Setters for Libint_t
    ! ************************************************************************

    subroutine libint1_set_libt_ab(ptr, itbl, setval) &
&                bind(c)
      import
      type(c_ptr) :: ptr
      integer(kind=c_int), value :: itbl
      real(kind=c_double), value :: setval
    end subroutine libint1_set_libt_ab

    subroutine libint1_set_libt_cd(ptr, itbl, setval) &
&                bind(c)
      import
      type(c_ptr) :: ptr
      integer(kind=c_int), value :: itbl
      real(kind=c_double), value :: setval
    end subroutine libint1_set_libt_cd

    ! ************************************************************************
    ! Setters for prim_data
    ! ************************************************************************

    subroutine libint1_set_prim_f(ptr, itbl, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      integer(kind=c_int), value :: itbl
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_f

    subroutine libint1_set_prim_oo2n(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_oo2n

    subroutine libint1_set_prim_oo2p(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_oo2p

    subroutine libint1_set_prim_oo2z(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_oo2z

    subroutine libint1_set_prim_oo2zn(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_oo2zn

    subroutine libint1_set_prim_pon(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_pon

    subroutine libint1_set_prim_poz(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_poz

    subroutine libint1_set_prim_ss_r12_ss(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_ss_r12_ss

    subroutine libint1_set_prim_twozeta_a(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_twozeta_a

    subroutine libint1_set_prim_twozeta_b(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_twozeta_b

    subroutine libint1_set_prim_twozeta_c(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_twozeta_c

    subroutine libint1_set_prim_twozeta_d(ptr, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_twozeta_d

    subroutine libint1_set_prim_u(ptr, itbl, jtbl, setval) &
&                  bind(c)
      import
      type(c_ptr) :: ptr
      integer(c_int), value :: itbl
      integer(c_int), value :: jtbl
      real(kind=c_double), value :: setval
    end subroutine libint1_set_prim_u

    ! ************************************************************************
    ! ERI routines
    ! ************************************************************************

    subroutine libint1_get_eris(ptr, n_a, n_b, n_c, n_d, eri_buffer) &
&                bind(c)
      import
      type(c_ptr) :: ptr
      integer(kind=c_int), value :: n_a, n_b, n_c, n_d
      type(c_ptr) :: eri_buffer
    end subroutine libint1_get_eris

    ! ************************************************************************
    ! Utility routines
    ! ************************************************************************

    function libint1_init(ptr, max_am) &
&                bind(c)
      import
      type(c_ptr) :: ptr
      integer(kind=c_int), value :: max_am
      integer(kind=c_int) :: libint1_init
    end function libint1_init

    subroutine libint1_free(ptr) &
&                bind(c)
      import
      type(c_ptr) :: ptr
    end subroutine libint1_free

  end interface

contains

  ! **************************************************************************
  !> \brief Initialises a Libint_t data structure
  !!
  !! This routines initialises a Libint_t data structure and allocates the
  !! corresponding memory.
  !!
  !! \note The routine stops the program if Libint fails to initialise.
  !!
  !! \param[out] libint_data: the wrapped Libint_t structure to set up
  !! \param[in] max_am: maximum angular momentum to consider for the ERIs
  ! **************************************************************************
  subroutine libint1f_init(libint_data, max_am)

    implicit none

    ! Arguments
    type(libint1f_t), intent(out) :: libint_data
    integer, intent(in) :: max_am

    ! Local variables
    integer :: ierr

    ! ------------------------------------------------------------------------

    ierr = libint1_init(libint_data%ptr, max_am)

    if ( ierr /= 0 ) stop 1

  end subroutine libint1f_init

  ! **************************************************************************
  !> \brief Initialises a Libint_t data structure
  !!
  !! This routines initialises a Libint_t data structure and allocates the
  !! corresponding memory.
  !!
  !! \note The routine stops the program if Libint fails to initialise.
  !!
  !! \param[inout] libint_data: the wrapped Libint_t structure to set up
  ! **************************************************************************
  subroutine libint1f_free(libint_data)

    implicit none

    ! Arguments
    type(libint1f_t), intent(inout) :: libint_data

    ! ------------------------------------------------------------------------

    call libint1_free(libint_data%ptr)

  end subroutine libint1f_free

  ! **************************************************************************
  !> \brief Retrieves PrimQuartet data from a Libint_t data structure
  !!
  !! \param[in] libint_data: the wrapped Libint_t structure to query (C)
  !! \param[out] prim_data: the prim_data structure to fill-in (Fortran)
  ! **************************************************************************
  subroutine libint1f_get_prim(libint_data, prim_data)

    implicit none

    ! Arguments
    type(libint1f_t), intent(in) :: libint_data
    type(libint1f_prim_t), intent(out) :: prim_data

    ! Local variables
    integer :: itbl, jtbl, n_prim_f

    ! ------------------------------------------------------------------------

    ! Retrieve scalar fields (one-to-one correspondence)
    prim_data%oo2n = libint1_get_prim_oo2n(libint_data%ptr)
    prim_data%oo2p = libint1_get_prim_oo2p(libint_data%ptr)
    prim_data%oo2z = libint1_get_prim_oo2z(libint_data%ptr)
    prim_data%oo2zn = libint1_get_prim_oo2zn(libint_data%ptr)
    prim_data%pon = libint1_get_prim_pon(libint_data%ptr)
    prim_data%poz = libint1_get_prim_poz(libint_data%ptr)
    prim_data%ss_r12_ss = libint1_get_prim_ss_r12_ss(libint_data%ptr)
    prim_data%twozeta_a = libint1_get_prim_twozeta_a(libint_data%ptr)
    prim_data%twozeta_b = libint1_get_prim_twozeta_b(libint_data%ptr)
    prim_data%twozeta_c = libint1_get_prim_twozeta_c(libint_data%ptr)
    prim_data%twozeta_d = libint1_get_prim_twozeta_d(libint_data%ptr)

    ! Retrieve all elements from U, translating indices from C to Fortran
    do itbl=1,6
      do jtbl=1,3
        prim_data%u(jtbl,itbl) = libint1_get_prim_u(libint_data%ptr, itbl-1, jtbl-1)
      end do
    end do

    ! This is the size of the F field of the prim_data structure
    n_prim_f = 4*(libint1_get_libint_max_am() - 1) + 1

    ! Retrieve all elements from F, translating indices from C to Fortran
    if ( .not. associated(prim_data%f) ) then
      allocate(prim_data%f(n_prim_f))
    end if
    do itbl=1,n_prim_f
      prim_data%f(itbl) = libint1_get_prim_f(libint_data%ptr, itbl-1)
    end do

  end subroutine libint1f_get_prim

  ! **************************************************************************
  !> \brief Propagates a Fortran data structure to the PrimQuartet field of
  !!        a Libint_t data structure
  !!
  !! \param[inout] libint_data: the wrapped Libint_t structure to update (C)
  !! \param[in] prim_data: the prim_data structure to load from (Fortran)
  ! **************************************************************************
  subroutine libint1f_set_prim(libint_data, prim_data)

    implicit none

    ! Arguments
    type(libint1f_t), intent(inout) :: libint_data
    type(libint1f_prim_t), intent(in) :: prim_data

    ! Local variables
    integer :: itbl, jtbl, n_prim_f
    real(kind=c_double) :: setval

    ! ------------------------------------------------------------------------

    ! Propagate scalar fields (one-to-one correspondence)
    setval = prim_data%oo2n
    call libint1_set_prim_oo2n(libint_data%ptr, setval)
    setval = prim_data%oo2p
    call libint1_set_prim_oo2p(libint_data%ptr, setval)
    setval = prim_data%oo2z
    call libint1_set_prim_oo2z(libint_data%ptr, setval)
    setval = prim_data%oo2zn
    call libint1_set_prim_oo2zn(libint_data%ptr, setval)
    setval = prim_data%pon
    call libint1_set_prim_pon(libint_data%ptr, setval)
    setval = prim_data%poz
    call libint1_set_prim_poz(libint_data%ptr, setval)
    setval = prim_data%ss_r12_ss
    call libint1_set_prim_ss_r12_ss(libint_data%ptr, setval)
    setval = prim_data%twozeta_a
    call libint1_set_prim_twozeta_a(libint_data%ptr, setval)
    setval = prim_data%twozeta_b
    call libint1_set_prim_twozeta_b(libint_data%ptr, setval)
    setval = prim_data%twozeta_c
    call libint1_set_prim_twozeta_c(libint_data%ptr, setval)
    setval = prim_data%twozeta_d
    call libint1_set_prim_twozeta_d(libint_data%ptr, setval)

    ! Propagate all elements to U, translating indices from Fortran to C
    do itbl=1,6
      do jtbl=1,3
        setval = prim_data%u(jtbl,itbl)
        call libint1_set_prim_u(libint_data%ptr, itbl-1, jtbl-1, setval)
      end do
    end do

    ! This is the size of the F field of the prim_data structure
    n_prim_f = 4*(libint1_get_libint_max_am() - 1) + 1

    ! Propagate all elements to F, translating indices from Fortran to C
    do itbl=1,n_prim_f
      if ( associated(prim_data%f) ) then
        setval = prim_data%f(itbl)
        call libint1_set_prim_f(libint_data%ptr, itbl-1, setval)
      else
        call libint1_set_prim_f(libint_data%ptr, itbl-1, 0.0_dp)
      end if
    end do

  end subroutine libint1f_set_prim

  ! **************************************************************************
  !> \brief Calculates the ERIs corresponding to the specified indices and
  !!        Libint_t data structure
  !!
  !! \param[inout] libint_data: the wrapped Libint_t structure to update (C)
  !! \param[in] n_a: index of A
  !! \param[in] n_b: index of B
  !! \param[in] n_c: index of C
  !! \param[in] n_d: index of D
  !! \param[in] eri_size: the number of ERIs to store
  !! \param[out] eri_store: array to store the calculated ERIs
  ! **************************************************************************
  subroutine libint1f_build_eris(libint_data, n_d, n_c, n_b, n_a, &
&                eri_size, eri_store)

    implicit none

    ! Arguments
    type(libint1f_t), intent(inout) :: libint_data
    integer, intent(in) :: n_a, n_b, n_c, n_d
    integer, intent(in) :: eri_size
    real(dp), intent(out) :: eri_store(eri_size)

    ! Local variables
    integer :: itbl
    real(kind=c_double), pointer :: eri_buffer(:) => null()
    type(c_ptr) :: eri_addr = C_NULL_PTR

    ! ------------------------------------------------------------------------

    ! Prepare the buffer
    allocate(eri_buffer(0:eri_size-1))
    eri_addr = c_loc(eri_buffer)

    ! Calculate the ERIs
    call libint1_get_eris(libint_data%ptr, n_a, n_b, n_c, n_d, eri_addr)

    ! Retrieve the ERI values
    do itbl = 1,eri_size
      eri_store(itbl) = eri_buffer(itbl-1)
    end do

    ! Clean-up the mess
    deallocate(eri_buffer)
    eri_buffer => null()
    eri_addr = C_NULL_PTR

  end subroutine libint1f_build_eris

end module libint1f
