! *** Module: nao2gto_libint ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Interface to the Libint and Libderiv libraries
!!
!! \author Manuel Guidon
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 11.2007 Created [Manuel Guidon]
!!      - 10.2009 Refactored [Manuel Guidon]
!!      - 01.2016 Imported to SIESTA and edited [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_libint

  use, intrinsic :: iso_c_binding

  implicit none

  private

  public :: &
&   Libderiv_t, &
&   Libint_t, &
&   build_deriv1_eri_size, &
&   build_eri_size, &
&   get_derivs, &
&   get_eris, &
&   initialize_libderiv, &
&   initialize_libint, &
&   libint_max_am, &
&   libderiv_max_am1, &
&   prim_data, &
&   prim_data_f_size, &
&   show_libint, &
&   terminate_libderiv, &
&   terminate_libint

  !> \brief Default Libint parameters
  !!
  !! \warning libint_max_am and libderiv_max_am1 are equal to the values
  !!          of Libint's configure script plus one.
  integer, parameter :: libint_max_am = 5
  integer, parameter :: libderiv_max_am1 = 4
  integer, parameter :: build_eri_size = libint_max_am-1
  integer, parameter :: build_deriv1_eri_size = libderiv_max_am1-1
  integer, parameter :: libint_vrr_classes_size = 2*(libint_max_am-1)+1
  integer, parameter :: libint_dvrr_classes_size = 2*(libderiv_max_am1-1)+1
  integer, parameter :: prim_data_f_size = 4*(libint_max_am-1)+1

  !> \brief Internal coefficients for Libint and Libderiv
  type, bind(c) :: prim_data
    real(c_double) :: F(prim_data_f_size)
    real(c_double) :: U(3,6)
    real(c_double) :: twozeta_a
    real(c_double) :: twozeta_b
    real(c_double) :: twozeta_c
    real(c_double) :: twozeta_d
    real(c_double) :: oo2z
    real(c_double) :: oo2n
    real(c_double) :: oo2zn
    real(c_double) :: poz
    real(c_double) :: pon
    real(c_double) :: oo2p
    real(c_double) :: ss_r12_ss
  end type prim_data

  !> \brief Libint data structure
  type, bind(c) :: Libint_t
    type(c_ptr)    :: int_stack
    type(c_ptr)    :: PrimQuartet
    real(c_double) :: AB(3)
    real(c_double) :: CD(3)
    type(c_ptr)    :: &
      vrr_classes(libint_vrr_classes_size,libint_vrr_classes_size)
    type(c_ptr)    :: vrr_stack
  end type Libint_t

  !> \brief Libderiv data structure
  type, bind(c) :: Libderiv_t
    type(c_ptr)    :: int_stack
    type(c_ptr)    :: primquartet
    type(c_ptr)    :: zero_stack
    type(c_ptr)    :: abcd(156)
    real(c_double) :: ab(3)
    real(c_double) :: cd(3)
    type(c_ptr)    :: &
&     deriv_classes(12,libint_dvrr_classes_size,libint_dvrr_classes_size)
    type(c_ptr)    :: &
&     deriv2_classes(144,libint_dvrr_classes_size,libint_dvrr_classes_size)
    type(c_ptr)    :: &
&     dvrr_classes(libint_dvrr_classes_size,libint_dvrr_classes_size)
    type(c_ptr)    :: dvrr_stack
  end type Libderiv_t

  !> Libint is built around a matrix of function pointers
  type(c_funptr), dimension( &
&   0:build_eri_size, 0:build_eri_size, &
&   0:build_eri_size, 0:build_eri_size), &
&   bind(c) :: build_eri

  !> Libderiv is built around a matrix of function pointers
  type(c_funptr), dimension( &
&   0:build_deriv1_eri_size, 0:build_deriv1_eri_size, &
    0:build_deriv1_eri_size, 0:build_deriv1_eri_size), &
&   bind(c) :: build_deriv1_eri

  interface

    function build(libint_data, max_num_prim_comb) bind(c)
      import
      type(c_ptr)                :: build
      type(Libint_t)             :: libint_data
      integer(kind=c_int), value :: max_num_prim_comb
    end function build

    subroutine build_deriv1(libderiv_data, max_num_prim_comb) bind(c)
      import
      type(Libderiv_t)           :: libderiv_data
      integer(kind=c_int), value :: max_num_prim_comb
    end subroutine build_deriv1

    subroutine free_libderiv(libderiv_data) bind(c, name="free_libderiv")
      import
      type(Libderiv_t) :: libderiv_data
    end subroutine free_libderiv

    subroutine free_libint(libint_data) bind(c, name="free_libint")
      import
      type(Libint_t) :: libint_data
    end subroutine free_libint

    function init_libderiv1(libderiv_data, max_am, max_num_prim_comb, ccs) &
&                bind(c, name="init_libderiv1")
      import
      integer(kind=c_int)        :: init_libderiv1
      type(Libderiv_t)           :: libderiv_data
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_comb
      integer(kind=c_int), value :: ccs
    end function init_libderiv1

    subroutine init_libderiv_base() bind(c, name="init_libderiv_base")
      import
    end subroutine init_libderiv_base

    function init_libint(libint_data, max_am, max_num_prim_comb) &
&                bind(c, name="init_libint")
      import
      integer(kind=c_int)        :: init_libint
      type(Libint_t)             :: libint_data
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_comb
    end function init_libint

    subroutine init_libint_base() bind(c, name="init_libint_base")
      import
    end subroutine init_libint_base

    function libderiv1_storage_required(max_am, max_num_prim_comb, ccs) &
&                bind(c, name="libderiv1_storage_required")
      import
      integer(kind=c_int)        :: libderiv1_storage_required
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_comb
      integer(kind=c_int), value :: ccs
    end function libderiv1_storage_required

    function libint_storage_required(max_am, max_num_prim_comb) &
&                bind(c, name="libint_storage_required")
      import
      integer(kind=c_int)        :: libint_storage_required
      integer(kind=c_int), value :: max_am
      integer(kind=c_int), value :: max_num_prim_comb
    end function libint_storage_required

  end interface

contains

! *****************************************************************************
!> \brief ...
!! \param n_d ...
!! \param n_c ...
!! \param n_b ...
!! \param n_a ...
!! \param libderiv_data ...
!! \param prim ...
!! \param work_forces ...
!! \param a_mysize: array shape to pass to c_f_pointer (must be a vector)
! *****************************************************************************
  subroutine get_derivs(n_d, n_c, n_b, n_a, libderiv_data, prim, work_forces, a_mysize)

    use precision, only: dp
    use atm_types, only: nco

    implicit none

    ! Arguments
    integer, intent(in) :: n_d, n_c, n_b, n_a
    type(Libderiv_t), intent(inout) :: libderiv_data
    type(prim_data), target, intent(inout) :: prim
    real(dp), dimension(nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d),12), &
&     intent(out) :: work_forces
    integer, intent(in) :: a_mysize(1)

    ! Local variables
    integer :: i, k
    type(c_ptr) :: pc_result
    procedure(build_deriv1), pointer :: pbuild_deriv1
    real(c_double), dimension(:), pointer :: tmp_data => null()

    ! -------------------------------------------------------------------------

    libderiv_data%primquartet = c_loc(prim)
    call c_f_procpointer(build_deriv1_eri(n_d, n_c, n_b, n_a), pbuild_deriv1)
    call pbuild_deriv1(libderiv_data, 1)

    do k=1,12
      if ( (k == 4) .or. (k == 5) .or. (k == 6) ) cycle
      pc_result = libderiv_data%abcd(k)
      call c_f_pointer(pc_result, tmp_data, a_mysize)
      do i=1,a_mysize(1)
        work_forces(i,k) = tmp_data(i)
      enddo
   end do

 end subroutine get_derivs

! *****************************************************************************
!> \brief ...
!!
!! \param n_d ...
!! \param n_c ...
!! \param n_b ...
!! \param n_a ...
!! \param libint_data ...
!! \param prim ...
!! \param p_work ...
!! \param[in] a_mysize: array shape to pass to c_f_pointer (must be a vector)
! *****************************************************************************
  subroutine get_eris(n_d, n_c, n_b, n_a, libint_data, prim, p_work, a_mysize)

    use precision, only: dp

    implicit none

    ! Arguments
    integer, intent(in) :: n_d, n_c, n_b, n_a
    type(Libint_t), intent(inout) :: libint_data
    type(prim_data), target, intent(inout) :: prim
    real(dp), dimension(:), pointer, intent(out) :: p_work
    integer, intent(in) :: a_mysize(1)

    ! Local variables
    integer :: i
    procedure(build), pointer :: pbuild
    type(c_ptr) :: pc_result = C_NULL_PTR

    ! -------------------------------------------------------------------------

    p_work => null()
    libint_data%primquartet = c_loc(prim)
    call c_f_procpointer(build_eri(n_d, n_c, n_b, n_a), pbuild)
    pc_result = pbuild(libint_data, 1)
    call c_f_pointer(pc_result, p_work, a_mysize)

   end subroutine get_eris

! *****************************************************************************
!> \brief Initializes a Libderiv data structure
!! \param[in,out] libderiv_data: Libderiv data structure
!! \param[in] max_am: maximum angular momentum to consider
! *****************************************************************************
  subroutine initialize_libderiv(libderiv_data, max_am)

    use sys, only: die
    use atm_types, only: nco

    implicit none

    ! Arguments
    type(Libderiv_t), intent(inout) :: libderiv_data
    integer, intent(in) :: max_am

    ! Local variables
    integer(kind=c_int) :: deriv_storage, max_am_local, max_classes, max_prim

    ! -------------------------------------------------------------------------

    max_am_local = max_am
    max_prim = 1
    max_classes = nco(max_am)**4

    call init_libderiv_base()
    deriv_storage = init_libderiv1(libderiv_data, max_am_local, max_prim, max_classes)

    if ( deriv_storage < 0 ) then
      call die("The angular momentum needed exceeds the value configured in libint")
    endif

  end subroutine initialize_libderiv

! *****************************************************************************
!> \brief Initializes a Libint data structure
!! \param[in,out] libint_data: Libint data structure
!! \param[in] max_am: maximum angular momentum to consider
! *****************************************************************************
  subroutine initialize_libint(libint_data, max_am)

    use sys, only: die

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data
    integer, intent(in)           :: max_am

    ! Local variables
    integer(kind=c_int) :: lib_storage, max_am_local, max_prim

    ! -------------------------------------------------------------------------

    max_am_local = max_am
    max_prim = 1

    call init_libint_base()
    lib_storage = init_libint(libint_data, max_am_local, max_prim)

    if ( lib_storage < 0 ) then
      call die("The angular momentum needed exceeds the value configured in libint")
    endif

  end subroutine initialize_libint

! *****************************************************************************
!> \brief Displays the current status of NAO2GTO data structures associated
!! to Libint
!!
!! \author Yann Pouillon
!!
!! \param[in] libint_data: Libint data structure
! *****************************************************************************
  subroutine show_libint(libint_data)

    use, intrinsic :: iso_c_binding

    implicit none

    ! Arguments
    type(Libint_t), intent(in) :: libint_data

    write(6, '(/,a,/)') "nao2gto_libint_status: Libint data structure ----------------------------------"
    write(6,'(a,/,"---",/,a,":")') "# *** YAML START ***", "libint_data"
    write(6,'(2x,a,":",1x,a)') "int_stack", status_string(libint_data%int_stack)
    write(6,'(2x,a,":",1x,a)') "primquartet", status_string(libint_data%primquartet)
    write(6,'(2x,a,":",1x,a)') "vrr_classes", status_string(libint_data%vrr_classes(1,1))
    write(6,'(2x,a,":",1x,a)') "vrr_stack", status_string(libint_data%vrr_stack)
    write(6,'(2x,a,":",1x,"[",3(1x,e12.3)," ]")') "ab", dble(libint_data%ab(:))
    write(6,'(2x,a,":",1x,"[",3(1x,e12.3)," ]")') "cd", dble(libint_data%cd(:))
    write(6,'(a)') "# *** YAML STOP ***"
    write(6, '(/,a,/)') "nao2gto_libint_status: END ----------------------------------------------------"

  contains

    function status_string(var)

      use, intrinsic :: iso_c_binding, only: c_associated, c_ptr

      implicit none

      type(c_ptr), intent(in) :: var
      character(len=4) :: status_string

      if ( c_associated(var) ) then
        status_string = "addr"
      else
        status_string = "null"
      endif

    end function status_string

  end subroutine show_libint

! *****************************************************************************
!> \brief Properly terminates a Libderiv data structure
!! \param[in,out] deriv: data structure to terminate
! *****************************************************************************
  subroutine terminate_libderiv(libderiv_data)

    implicit none

    ! Arguments
    type(Libderiv_t), intent(inout) :: libderiv_data

    ! -------------------------------------------------------------------------

    call free_libderiv(libderiv_data)

  end subroutine terminate_libderiv

! *****************************************************************************
!> \brief Properly terminates a Libint data structure
!! \param[in,out] libint_data: data structure to terminate
! *****************************************************************************
  subroutine terminate_libint(libint_data)

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data

    ! -------------------------------------------------------------------------

    call free_libint(libint_data)

  end subroutine terminate_libint

end module nao2gto_libint
