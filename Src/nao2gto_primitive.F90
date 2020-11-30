! *** Module: nao2gto_primitive ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Interface to the Libint and Libderiv libraries
!!
!!  This module was originally extracted from CP2K: A general program to
!!  perform molecular dynamics simulations. It was then adapted for use within
!!  SIESTA.
!!
!! \note
!!      This file currently works with a version of Libint configured for
!!      LIBINT_MAX_AM=5 and LIBINT_MAX_AM1=4 only.
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \copyright
!!      - 2000-2009 CP2K Developers Group
!!      - 2010-2018 SIESTA Developers Group
!!
!! \par History
!!      - 11.2006 Created [Manuel Guidon]
!!      - 01.2018 Refactored for SIESTA [Yann Pouillon]
! *****************************************************************************
module nao2gto_primitive

  use nao2gto_common

  implicit none

  private

  integer, parameter :: full_perm1(12) = [1,2,3,4,5,6,7,8,9,10,11,12]
  integer, parameter :: full_perm2(12) = [4,5,6,1,2,3,7,8,9,10,11,12]
  integer, parameter :: full_perm3(12) = [1,2,3,4,5,6,10,11,12,7,8,9]
  integer, parameter :: full_perm4(12) = [4,5,6,1,2,3,10,11,12,7,8,9]
  integer, parameter :: full_perm5(12) = [7,8,9,10,11,12,1,2,3,4,5,6]
  integer, parameter :: full_perm6(12) = [7,8,9,10,11,12,4,5,6,1,2,3]
  integer, parameter :: full_perm7(12) = [10,11,12,7,8,9,1,2,3,4,5,6]
  integer, parameter :: full_perm8(12) = [10,11,12,7,8,9,4,5,6,1,2,3]

  public :: &
&   calc_primitive_deriv_eri, &
&   calc_primitive_eri, &
&   calc_primitive_eri2, &
&   calc_primitive_screen

contains

! *****************************************************************************
!> \brief Evaluates derivatives of electron repulsion integrals for a primitive
!!        quartet.
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \par History
!!      - 03.2007 Created [Manuel Guidon]
!!      - 08.2007 Refactured permutation part [Manuel Guidon]
!!      - 01.2018 Changed interface for SIESTA [Yann Pouillon]
!!
!! \param[in,out] deriv_data: Libderiv data structure
!! \param[in] A: ...
!! \param[in] B: ...
!! \param[in] C: ...
!! \param[in] D: ...
!! \param[in] ZetaA: ...
!! \param[in] ZetaB: ...
!! \param[in] ZetaC: ...
!! \param[in] ZetaD: ...
!! \param[in] n_a: ...
!! \param[in] n_b: ...
!! \param[in] n_c: ...
!! \param[in] n_d: ...
!! \param[in,out] work_forces: ...
!! \param[in] ncoa: ...
!! \param[in] ncob: ...
!! \param[in] ncoc: ...
!! \param[in] ncod: ...
!! \param[out] primitive_forces: ...
!! \param[in] hfx_opts: ...
!! \param[in] offset_a: ...
!! \param[in] offset_b: ...
!! \param[in] offset_c: ...
!! \param[in] offset_d: ...
!! \param[in] ZetaInv: ...
!! \param[in] EtaInv: ...
!! \param[in] ZetapEtaInv: ...
!! \param[in] Rho: ...
!! \param[in] RhoInv: ...
!! \param[in] S1234: ...
!! \param[in] P: ...
!! \param[in] Q: ...
!! \param[in] W: ...
!! \param[in] rpq2: ...
!! \param[in] AB: ...
!! \param[in] CD: ...
! *****************************************************************************
  subroutine calc_primitive_deriv_eri( &
&                deriv_data, A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
&                n_a, n_b, n_c, n_d, work_forces, ncoa, ncob, ncoc, ncod, &
&                primitive_forces, hfx_opts, max_contraction, tmp_max_all, &
&                offset_a, offset_b, offset_c,offset_d, &
&                ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
&                P, Q, W, rpq2, AB, CD)

    use atm_types,      only: nco
    use nao2gto_libint, only: Libderiv_t, get_derivs, prim_data
    use nao2gto_types,  only: hfx_options_type

    implicit none

    ! Arguments
    type(Libderiv_t), intent(inout)       :: deriv_data
    real(dp), intent(in)                  :: A(3), B(3), C(3), D(3), &
&                                            Zeta_A, Zeta_B, Zeta_C, Zeta_D
    integer, intent(in)                   :: n_a, n_b, n_c, n_d
    real(dp),&
&     dimension(nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d), 12), &
&     intent(inout)                       :: work_forces
    integer, intent(in)                   :: ncoa, ncob, ncoc, ncod
    real(dp), &
&     dimension(ncoa, ncob, ncoc, ncod, 12), & 
&     intent(out)                         :: primitive_forces
    real(dp), intent(in) :: max_contraction
    real(dp), intent(out) :: tmp_max_all
    integer, intent(in)                   :: offset_a, offset_b, offset_c, &
&                                            offset_d
    type(hfx_options_type), intent(in) :: hfx_opts
    real(dp), intent(in)                  :: ZetaInv, EtaInv, ZetapEtaInv, &
&                                            Rho, RhoInv, S1234, P(3), Q(3),&
&                                            W(3), rpq2, AB(3), CD(3)

    ! Local variables
    integer                 :: a_mysize(1), i, j, k, l, m_max, mysize, n, &
&                              p1, p2, p3, perm_case
    logical                 :: do_it
    real(dp) :: tmp_max
    type(prim_data), target :: prim

    ! -------------------------------------------------------------------------

    ! Permutation of configurations
    perm_case = 1
    if(n_a<n_b) then
      perm_case = perm_case + 1
    end if
    if(n_c<n_d) then
      perm_case = perm_case + 2
    end if
    if( n_a+n_b > n_c+n_d) then
      perm_case = perm_case + 4
    end if

    do_it = .true.
    m_max = n_a+n_b+n_c+n_d
    m_max = m_max + 1
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize=mysize

    select case(perm_case)

      case(1)
        call build_deriv_data(A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
&         m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, P, Q, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB = AB ! A-B
        deriv_data%CD = CD ! C-D
        call get_derivs(n_d, n_c, n_b, n_a, deriv_data, prim, work_forces, a_mysize)
        do k=4,6
          do i=1,mysize
             work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                          work_forces(i,k+3) + &
&                                          work_forces(i,k+6))
          end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do i = 1,nco(n_a)
            p1 = (i-1) *nco(n_b)
            do j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_c)
              do k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                do l = 1,nco(n_d)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm1(n)) = work_forces(p3+l,n)
                end do
              end do
            end do
          end do
        end do

      case(2)
        call build_deriv_data(B, A, C, D, Zeta_B, Zeta_A, Zeta_C, Zeta_D, &
&         m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, P, Q, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB = -AB ! B-A
        deriv_data%CD =  CD ! C-D
        call get_derivs(n_d, n_c, n_a, n_b, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
          do i=1,mysize
             work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                          work_forces(i,k+3) + &
&                                          work_forces(i,k+6))
          end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            do i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_c)
              do k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                do l = 1,nco(n_d)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm2(n)) = work_forces(p3+l,n)
                end do
              end do
            end do
          end do
        end do

     case(3)
        call build_deriv_data(A, B, D, C, Zeta_A, Zeta_B, Zeta_D, Zeta_C, &
&         m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, P, Q, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB =  AB ! A-B
        deriv_data%CD = -CD ! D-C
        call get_derivs(n_c, n_d, n_b, n_a, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp* (work_forces(i,k-3) + &
&                                          work_forces(i,k+3) + &
&                                          work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do i = 1,nco(n_a)
            p1 = (i-1)*nco(n_b)
            do j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_d)
              do l = 1,nco(n_d)
                p3 = (p2 + l-1) * nco(n_c)
                do k = 1,nco(n_c)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm3(n)) = work_forces(p3+k,n)
                end do
              end do
            end do
          end do
        end do

      case(4)
        call build_deriv_data(B, A, D, C, Zeta_B, Zeta_A, Zeta_D, Zeta_C, &
&         m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, P, Q, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB = -AB ! B-A
        deriv_data%CD = -CD ! D-C
        call get_derivs(n_c, n_d, n_a, n_b, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                           work_forces(i,k+3) + &
&                                           work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            do i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_d)
              do l = 1,nco(n_d)
                p3 = (p2 + l-1)*nco(n_c)
                do k = 1,nco(n_c)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm4(n)) = work_forces(p3+k,n)
                end do
              end do
            end do
          end do
        end do

      case(5)
        call build_deriv_data(C, D, A, B, Zeta_C, Zeta_D, Zeta_A, Zeta_B, &
&          m_max, hfx_opts, prim, do_it, EtaInv, ZetaInv, ZetapEtaInv, &
&          Rho, RhoInv, S1234, Q, P, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB = CD ! C-D
        deriv_data%CD = AB ! A-B
        call get_derivs(n_b, n_a, n_d, n_c, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                           work_forces(i,k+3) + &
&                                           work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            do l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_a)
              do i = 1,nco(n_a)
                p3 = (p2 + i-1)*nco(n_b)
                do j = 1,nco(n_b)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm5(n)) = work_forces(p3+j,n)
                end do
              end do
            end do
          end do
        end do

      case(6)
        call build_deriv_data(C, D, B, A, Zeta_C, Zeta_D, Zeta_B, Zeta_A, &
&         m_max, hfx_opts, prim, do_it, EtaInv, ZetaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, Q, P, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB =  CD ! C-D
        deriv_data%CD = -AB ! B-A
        call get_derivs(n_a, n_b, n_d, n_c, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
          do i=1,mysize
             work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                          work_forces(i,k+3) + &
&                                          work_forces(i,k+6))
          end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            do l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_b)
              do j = 1,nco(n_b)
                p3 = (p2 + j-1)*nco(n_a)
                do i = 1,nco(n_a)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm6(n)) = work_forces(p3+i,n)
                end do
              end do
            end do
          end do
        end do

      case(7)
        call build_deriv_data(D, C, A, B, Zeta_D, Zeta_C, Zeta_A, Zeta_B, &
&         m_max, hfx_opts, prim, do_it, EtaInv, ZetaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, Q, P, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB=-CD ! D-C
        deriv_data%CD=AB ! A-B
        call get_derivs(n_b, n_a, n_c, n_d, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                           work_forces(i,k+3) + &
&                                           work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            do k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_a)
              do i = 1,nco(n_a)
                p3 = (p2 + i-1) *nco(n_b)
                do j = 1,nco(n_b)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm7(n)) = work_forces(p3+j,n)
                end do
              end do
            end do
          end do
        end do

      case(8)
        call build_deriv_data(D, C, B, A, Zeta_D, Zeta_C, Zeta_B, Zeta_A, &
&         m_max, hfx_opts, prim, do_it, EtaInv, ZetaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, Q, P, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB=-CD ! D-C
        deriv_data%CD=-AB ! B-A
        call get_derivs(n_a, n_b, n_c, n_d, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                           work_forces(i,k+3) + &
&                                           work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            do k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_b)
              do j = 1,nco(n_b)
                p3 = (p2 + j-1) * nco(n_a)
                do i = 1,nco(n_a)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm8(n)) = work_forces(p3+i,n)
                end do
              end do
            end do
          end do
        end do
    end select

  end subroutine calc_primitive_deriv_eri

! *****************************************************************************
!> \brief Evaluates electron repulsion integrals for a primitive quartet.
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \par History
!!      - 11.2006 Created [Manuel Guidon]
!!      - 08.2007 Refactured permutation part [Manuel Guidon]
!!      - 01.2018 Changed interface for SIESTA [Yann Pouillon]
!!
!! \param[in,out] lib: Libint data structure
!! \param[in] A: ...
!! \param[in] B: ...
!! \param[in] C: ...
!! \param[in] D: ...
!! \param[in] n_a: ...
!! \param[in] n_b: ...
!! \param[in] n_c: ...
!! \param[in] n_d: ...
!! \param[in] ncoa: ...
!! \param[in] ncob: ...
!! \param[in] ncoc: ...
!! \param[in] ncod: ...
!! \param[in] offset_a: ...
!! \param[in] offset_b: ...
!! \param[in] offset_c: ...
!! \param[in] offset_d: ...
!! \param[out] primitives: ...
!! \param[in] hfx_opts: ...
!! \param[in] ZetaInv: ...
!! \param[in] EtaInv: ...
!! \param[in] ZetapEtaInv: ...
!! \param[in] Rho: ...
!! \param[in] RhoInv: ...
!! \param[in] S1234: ...
!! \param[in] P: ...
!! \param[in] Q: ...
!! \param[in] W: ...
!! \param[in] rpq2: ...
!! \param[in] AB: ...
!! \param[in] CD: ...
! *****************************************************************************
  subroutine calc_primitive_eri(lib, A, B, C, D, n_a, n_b, n_c, n_d, ncoa, &
&                ncob, ncoc, ncod, offset_a, offset_b, offset_c, offset_d, &
&                primitives, hfx_opts, max_contraction, tmp_max, neris, &
&                ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, &
&                rpq2, AB, CD, R1, R2)

    use atm_types,      only: nco
    use alloc,          only: de_alloc, re_alloc
    use nao2gto_libint, only: Libint_t, get_eris, prim_data
    use nao2gto_types,  only: hfx_options_type

    implicit none

    ! Arguments
    type(Libint_t), intent(inout)         :: lib
    real(dp), intent(in)                  :: A(3), B(3), C(3), D(3)
    integer, intent(in)                   :: n_a, n_b, n_c, n_d, ncoa, &
&                                            ncob, ncoc, ncod, offset_a, &
&                                            offset_b, offset_c, offset_d
    real(dp), &
&     dimension(ncoa, ncob, ncoc, ncod), &
&     intent(out)                         :: primitives
    type(hfx_options_type), intent(in) :: hfx_opts
    real(dp), intent(in) :: max_contraction
    real(dp), intent(out) :: tmp_max
    integer(int_8), intent(inout) :: neris
    real(dp), intent(in)                  :: ZetaInv, EtaInv, &
&                                            ZetapEtaInv, Rho, RhoInv, &
&                                            S1234, P(3), Q(3), W(3), &
&                                            rpq2, AB(3), CD(3), R1, R2

    ! Local variables
    integer                         :: a_mysize(1), i, j, k, l, &
&                                      m_max, mysize, p1, p2, p3, &
&                                      p4, perm_case, q1, q2, q3, q4
    logical                         :: do_it
    real(dp), dimension(:), pointer :: p_work => null()
    type(prim_data), target         :: prim

    ! -------------------------------------------------------------------------

    m_max = n_a + n_b + n_c + n_d
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize(:) = mysize

    do_it = .true.
    if(m_max/=0) then
      perm_case = 1
      if(n_a<n_b) then
        perm_case = perm_case + 1
      end if
      if(n_c<n_d) then
        perm_case = perm_case + 2
      end if
      if(n_a+n_b > n_c+n_d) then
        perm_case = perm_case + 4
      end if

      select case(perm_case)

        case(1)
          call build_quartet_data(A, B, C, D, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
          if( .not. do_it ) return
          lib%AB=AB ! A-B
          lib%CD=CD ! C-D
          call get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do i = 1,nco(n_a)
            p1 = (i-1) *nco(n_b)
            q1 = offset_a + i
            do j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_c)
              q2 = offset_b + j
              do k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                q3 = offset_c + k
                do l = 1,nco(n_d)
                  q4 = offset_d + l
                  p4 = p3 + l
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(2)
          call build_quartet_data(B, A, C, D, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
          if( .not. do_it ) return
          lib%AB=-AB ! B-A
          lib%CD=CD ! C-D
          call get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            q2 = offset_b+j
            do i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_c)
              q1 = offset_a+i
              do k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                q3 = offset_c+k
                do l = 1,nco(n_d)
                  q4 = offset_d+l
                  p4 = p3 + l
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(3)
          call build_quartet_data(A, B, D, C, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
          if( .not. do_it ) return
          lib%AB=AB ! A-B
          lib%CD=-CD ! D-C
          call get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do i = 1,nco(n_a)
            p1 = (i-1)*nco(n_b)
            q1 = offset_a + i
            do j = 1,nco(n_b)
              q2 = offset_b + j
              p2 = (p1 + j-1)*nco(n_d)
              do l = 1,nco(n_d)
                q4 = offset_d + l
                p3 = (p2 + l-1) * nco(n_c)
                do k = 1,nco(n_c)
                  q3 = offset_c + k
                  p4 = p3+k
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(4)
          call build_quartet_data(B, A, D, C, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
          if( .not. do_it ) return
          lib%AB=-AB ! B-A
          lib%CD=-CD ! D-C
          call get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            q2 = offset_b + j
            do i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_d)
              q1 = offset_a + i
              do l = 1,nco(n_d)
                p3 = (p2 + l-1)*nco(n_c)
                q4 = offset_d + l
                do k = 1,nco(n_c)
                  q3 = offset_c + k
                  p4 = p3 + k
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(5)
          call build_quartet_data(C, D, A, B, m_max, &
                            hfx_opts, prim, do_it, &
                            EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            Q, P, W, rpq2)
          if( .not. do_it ) return
          lib%AB=CD ! C-D
          lib%CD=AB ! A-B
          call get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do k = 1,nco(n_c)
            q3 = offset_c + k
            p1 = (k-1)*nco(n_d)
            do l = 1,nco(n_d)
              q4 = offset_d + l
              p2 = (p1 + l-1)*nco(n_a)
              do i = 1,nco(n_a)
                q1 = offset_a + i
                p3 = (p2 + i-1)*nco(n_b)
                do j = 1,nco(n_b)
                  q2 = offset_b + j
                  p4 = p3+j
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(6)
          call build_quartet_data(C, D, B, A, m_max, &
                            hfx_opts, prim, do_it, &
                            EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            Q, P, W, rpq2)
          if( .not. do_it ) return
          lib%AB=CD ! C-D
          lib%CD=-AB ! B-A
          call get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            q3 = offset_c + k
            do l = 1,nco(n_d)
              q4 = offset_d + l
              p2 = (p1 + l-1)*nco(n_b)
              do j = 1,nco(n_b)
                q2 = offset_b + j
                p3 = (p2 + j-1)*nco(n_a)
                do i = 1,nco(n_a)
                  p4 = p3+i
                  q1 = offset_a + i
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(7)
          call build_quartet_data(D, C, A, B, m_max, &
                            hfx_opts, prim, do_it, &
                            EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            Q, P, W, rpq2)
          if( .not. do_it ) return
          lib%AB=-CD ! D-C
          lib%CD=AB ! A-B
          call get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            q4 = offset_d + l
            do k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_a)
              q3 = offset_c + k
              do i = 1,nco(n_a)
                p3 = (p2 + i-1) *nco(n_b)
                q1 = offset_a + i
                do j = 1,nco(n_b)
                  q2 = offset_b + j
                  p4 = p3+j
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(8)
          call build_quartet_data(D, C, B, A, m_max, &
                            hfx_opts, prim, do_it, &
                            EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            Q, P, W, rpq2)
          if( .not. do_it ) return
          lib%AB=-CD ! D-C
          lib%CD=-AB ! B-A
          call get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do l = 1,nco(n_d)
            q4 = offset_d + l
            p1 = (l-1)*nco(n_c)
            do k = 1,nco(n_c)
              q3 = offset_c + k
              p2 = (p1 + k-1) * nco(n_b)
              do j = 1,nco(n_b)
                q2 = offset_b + j
                p3 = (p2 + j-1) * nco(n_a)
                do i = 1,nco(n_a)
                  q1 = offset_a + i
                  p4 = p3 + i
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

      end select

    else

      call build_quartet_data(A, B, C, D, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
      neris = neris + 1

      if( .not. do_it ) return

      primitives(offset_a+1,offset_b+1,offset_c+1,offset_d+1) = prim%F(1)
      tmp_max = max_contraction*abs(prim%F(1))

    end if

  end subroutine calc_primitive_eri

! *****************************************************************************
!> \brief sets the threshold for truncated calculations
!> \par History
!>      12.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
!  SUBROUTINE set_eps_cutoff(eps_schwarz)
!    REAL(dp), INTENT(IN)                     :: eps_schwarz
!
!    eps_cutoff = SQRT(-LOG(eps_schwarz))
!  END SUBROUTINE set_eps_cutoff

! *****************************************************************************
  SUBROUTINE calc_primitive_eri2(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                          n_a,n_b,n_c,n_d,&
                          ncoa,ncob,ncoc,ncod,&
                          offset_a,offset_b,offset_c,offset_d,&
                          primitives, hfx_opts,    &
                          Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                          P,Q,W,rpq2,AB,CD)

    use atm_types, only: nco
    use nao2gto_libint
    use nao2gto_types

    implicit none
    
    TYPE(Libint_t)                            :: lib
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d, ncoa, &
                                                ncob, ncoc, ncod, offset_a, &
                                                offset_b, offset_c, offset_d
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod)      :: primitives
    TYPE(hfx_options_type)                :: hfx_opts

    REAL(dp), INTENT(IN)                     :: Zeta, ZetaInv, Eta, EtaInv, &
                                                ZetapEtaInv, Rho, RhoInv, &
                                                S1234, P(3), Q(3), W(3), &
                                                rpq2, AB(3), CD(3)

    INTEGER                                  :: a_mysize(1), i, j, k, l, &
                                                m_max, mysize, p1, p2, p3, &
                                                p4, perm_case, q1, q2, q3, q4
    LOGICAL                                  :: do_it
    REAL(dp), DIMENSION(:), POINTER          :: p_work => null()
    TYPE(prim_data), TARGET                  :: prim
  
    m_max = n_a+n_b+n_c+n_d
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize = mysize

    do_it = .TRUE.
    IF(m_max/=0) THEN
      perm_case = 1
      IF(n_a<n_b) THEN
        perm_case = perm_case + 1
      END IF
      IF(n_c<n_d) THEN
        perm_case = perm_case + 2
      END IF
      IF(n_a+n_b > n_c+n_d) THEN
        perm_case = perm_case + 4
      END IF

      SELECT CASE(perm_case)
        CASE(1)
          CALL build_quartet_data(A,B,C,D, m_max,&
                            hfx_opts, prim, do_it,  &
                            ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=AB!A-B
          lib%CD=CD!C-D
          CALL get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)

          DO i = 1,nco(n_a)
            p1 = (i-1) *nco(n_b)
            q1 = offset_a + i
            DO j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_c)
              q2 = offset_b + j
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                q3 = offset_c + k
                DO l = 1,nco(n_d)
                  q4 = offset_d + l
                  p4 = p3 + l
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(2)
          CALL build_quartet_data(B,A,C,D, m_max,&
                            hfx_opts, prim, do_it, &
                            ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=-AB!B-A
          lib%CD=CD!C-D
          CALL get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            q2 = offset_b+j
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_c)
              q1 = offset_a+i
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                q3 = offset_c+k
                DO l = 1,nco(n_d)
                  q4 = offset_d+l
                  p4 = p3 + l
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(3)
          CALL build_quartet_data(A,B,D,C, m_max,&
                            hfx_opts, prim, do_it, &
                            ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=AB!A-B
          lib%CD=-CD!D-C
          CALL get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)

          DO i = 1,nco(n_a)
            p1 = (i-1)*nco(n_b)
            q1 = offset_a + i
            DO j = 1,nco(n_b)
              q2 = offset_b + j
              p2 = (p1 + j-1)*nco(n_d)
              DO l = 1,nco(n_d)
                q4 = offset_d + l
                p3 = (p2 + l-1) * nco(n_c)
                DO k = 1,nco(n_c)
                  q3 = offset_c + k
                  p4 = p3+k
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(4)
          CALL build_quartet_data(B,A,D,C, m_max,&
                            hfx_opts, prim, do_it, &
                            ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=-AB!B-A
          lib%CD=-CD!D-C
          CALL get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            q2 = offset_b + j
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_d)
              q1 = offset_a + i
              DO l = 1,nco(n_d)
                p3 = (p2 + l-1)*nco(n_c)
                q4 = offset_d + l
                DO k = 1,nco(n_c)
                  q3 = offset_c + k
                  p4 = p3 + k
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(5)
          CALL build_quartet_data(C,D,A,B, m_max,&
                            hfx_opts, prim, do_it ,&
                            EtaInv,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            Q,P,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=CD!C-D
          lib%CD=AB!A-B
          CALL get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)

          DO k = 1,nco(n_c)
            q3 = offset_c + k
            p1 = (k-1)*nco(n_d)
            DO l = 1,nco(n_d)
              q4 = offset_d + l
              p2 = (p1 + l-1)*nco(n_a)
              DO i = 1,nco(n_a)
                q1 = offset_a + i
                p3 = (p2 + i-1)*nco(n_b)
                DO j = 1,nco(n_b)
                  q2 = offset_b + j
                  p4 = p3+j
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(6)
          CALL build_quartet_data(C,D,B,A, m_max,&
                            hfx_opts, prim, do_it ,&
                            EtaInv,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            Q,P,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=CD!C-D
          lib%CD=-AB!B-A
          CALL get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)

          DO k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            q3 = offset_c + k
            DO l = 1,nco(n_d)
              q4 = offset_d + l
              p2 = (p1 + l-1)*nco(n_b)
              DO j = 1,nco(n_b)
                q2 = offset_b + j
                p3 = (p2 + j-1)*nco(n_a)
                DO i = 1,nco(n_a)
                  p4 = p3+i
                  q1 = offset_a + i
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(7)
          CALL build_quartet_data(D,C,A,B, m_max,&
                            hfx_opts, prim, do_it ,&
                            EtaInv,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            Q,P,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=-CD!D-C
          lib%CD=AB!A-B
          CALL get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)

          DO l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            q4 = offset_d + l
            DO k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_a)
              q3 = offset_c + k
              DO i = 1,nco(n_a)
                p3 = (p2 + i-1) *nco(n_b)
                q1 = offset_a + i
                DO j = 1,nco(n_b)
                  q2 = offset_b + j
                  p4 = p3+j
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(8)
          CALL build_quartet_data(D,C,B,A, m_max,&
                            hfx_opts, prim, do_it ,&
                            EtaInv,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            Q,P,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=-CD!D-C
          lib%CD=-AB!B-A
          CALL get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)

          DO l = 1,nco(n_d)
            q4 = offset_d + l
            p1 = (l-1)*nco(n_c)
            DO k = 1,nco(n_c)
              q3 = offset_c + k
              p2 = (p1 + k-1) * nco(n_b)
              DO j = 1,nco(n_b)
                q2 = offset_b + j
                p3 = (p2 + j-1) * nco(n_a)
                DO i = 1,nco(n_a)
                  q1 = offset_a + i
                  p4 = p3 + i
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
      END SELECT
    ELSE
      CALL build_quartet_data(A,B,C,D, m_max,&
                            hfx_opts, prim, do_it, &
                            ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
      IF( .NOT. do_it ) RETURN
      primitives(offset_a+1,offset_b+1,offset_c+1,offset_d+1) = prim%F(1)
    END IF

  END SUBROUTINE calc_primitive_eri2

! *****************************************************************************
!> \brief Evaluate electron repulsion integrals for a primitive quartet
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon

! *****************************************************************************
  SUBROUTINE calc_primitive_screen(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                                   n_a,n_b,n_c,n_d,&
                                   max_val, hfx_opts,R1, R2)

    use atm_types, only: nco
    use nao2gto_libint
    use nao2gto_types

    implicit none

    TYPE(Libint_t)                            :: lib
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d
    REAL(dp), INTENT(INOUT)                  :: max_val
    TYPE(hfx_options_type)                :: hfx_opts
    REAL(dp)                                 :: R1, R2

    INTEGER                                  :: a_mysize(1), i, m_max, &
                                                mysize, perm_case
    LOGICAL                                  :: do_it
    REAL(dp), DIMENSION(:), POINTER          :: p_work => null()
    TYPE(prim_data), TARGET                  :: prim

!permutation of configuration

    m_max = n_a+n_b+n_c+n_d
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize = mysize

    do_it = .TRUE.
    IF(m_max/=0) THEN
      perm_case = 1
      IF(n_a<n_b) THEN
        perm_case = perm_case + 1
      END IF
      IF(n_c<n_d) THEN
        perm_case = perm_case + 2
      END IF
      IF(n_a+n_b > n_c+n_d) THEN
        perm_case = perm_case + 4
      END IF

      SELECT CASE(perm_case)
        CASE(1)
          CALL build_quartet_data_screen(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                            hfx_opts, prim, do_it, n_a+n_b, n_c+n_d, R1, R2)
          lib%AB=A-B
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(2)
          CALL build_quartet_data_screen(B,A,C,D,Zeta_B, Zeta_A, Zeta_C, Zeta_D, m_max,&
                            hfx_opts, prim, do_it, n_b+n_a, n_c+n_d, R1,R2)
          lib%AB=B-A
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)

          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(3)
          CALL build_quartet_data_screen(A,B,D,C,Zeta_A, Zeta_B, Zeta_D, Zeta_C, m_max,&
                            hfx_opts, prim, do_it, n_a+n_b, n_d+n_c, R1,R2)
          lib%AB=A-B
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(4)
          CALL build_quartet_data_screen(B,A,D,C,Zeta_B, Zeta_A, Zeta_D, Zeta_C, m_max,&
                            hfx_opts, prim, do_it, n_b+n_a, n_d+n_c, R1,R2)
          lib%AB=B-A
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(5)
          CALL build_quartet_data_screen(C,D,A,B,Zeta_C, Zeta_D, Zeta_A, Zeta_B, m_max,&
                            hfx_opts, prim, do_it, n_c+n_d, n_a+n_b, R1,R2)
          lib%AB=C-D
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(6)
          CALL build_quartet_data_screen(C,D,B,A,Zeta_C, Zeta_D, Zeta_B, Zeta_A, m_max,&
                            hfx_opts, prim, do_it, n_c+n_d, n_b+n_a, R1,R2)
          lib%AB=C-D
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(7)
          CALL build_quartet_data_screen(D,C,A,B,Zeta_D, Zeta_C, Zeta_A, Zeta_B, m_max,&
                            hfx_opts, prim, do_it, n_d+n_c, n_a+n_b, R1,R2)
          lib%AB=D-C
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(8)
          CALL build_quartet_data_screen(D,C,B,A,Zeta_D, Zeta_C, Zeta_B, Zeta_A, m_max,&
                            hfx_opts, prim, do_it, n_d+n_c, n_b+n_a, R1,R2)
          lib%AB=D-C
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
      END SELECT

    ELSE
      CALL build_quartet_data_screen(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                              hfx_opts, prim, do_it, 0, 0, R1,R2)
      max_val = ABS(prim%F(1))
    END IF
  END SUBROUTINE calc_primitive_screen

! *****************************************************************************
! *** Private routines                                                      ***
! *****************************************************************************

! *****************************************************************************
!> \brief Fills the data structure used in Libderiv.
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \par History
!!      - 03.2007 Created [Manuel Guidon]
!!      - 01.2018 Changed interface for SIESTA [Yann Pouillon]
!!
!! \param[in,out] lib: Libint data structure
!! \param[in] A: ...
!! \param[in] B: ...
!! \param[in] C: ...
!! \param[in] D: ...
!! \param[in] ZetaA: ...
!! \param[in] ZetaB: ...
!! \param[in] ZetaC: ...
!! \param[in] ZetaD: ...
!! \param[in] m_max: ...
!! \param[in] hfx_opts: ...
!! \param[in,out] prim: ...
!! \param[in,out] do_it: ...
!! \param[in] ZetaInv: ...
!! \param[in] EtaInv: ...
!! \param[in] ZetapEtaInv: ...
!! \param[in] Rho: ...
!! \param[in] RhoInv: ...
!! \param[in] S1234: ...
!! \param[in] P: ...
!! \param[in] Q: ...
!! \param[in] W: ...
!! \param[in] PQ2: ...
! *****************************************************************************
  subroutine build_deriv_data(A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
&                m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, &
&                ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, PQ2)

    use units,          only: pi
    use nao2gto_common
    use nao2gto_libint, only: prim_data
    use nao2gto_utils,  only: fgamma => fgamma_0
    use nao2gto_types,  only: hfx_options_type

    implicit none

    ! Arguments
    real(dp), intent(in)               :: A(3), B(3), C(3), D(3)
    real(dp), intent(in)               :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    integer, intent(in)                :: m_max
    type(hfx_options_type), intent(in) :: hfx_opts
    type(prim_data), intent(out)       :: prim
    logical, intent(inout)             :: do_it
    real(dp), intent(in)               :: ZetaInv, EtaInv, ZetapEtaInv, &
&                                         Rho, RhoInv, S1234, P(3), Q(3), &
&                                         W(3), PQ2

    ! Local variables
    integer                 :: i
    real(dp)                :: factor, omega2, omega_corr, omega_corr2, T, tmp
    real(dp), dimension(17) :: Fm

    ! -------------------------------------------------------------------------

    T = Rho*PQ2
    do_it = .true.

    select case (hfx_opts%potential_type)

      case(do_hfx_potential_coulomb)
        call fgamma(m_max, T, prim%F(1))
        factor = 2.0_dp*Pi*RhoInv

      case(do_hfx_potential_short)
        call fgamma(m_max, T, prim%F)
        omega2 = hfx_opts%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        call fgamma(m_max,T,Fm)
        tmp = - omega_corr
        do i=1,m_max+1
          prim%F(i) = prim%F(i) + Fm(i)*tmp
          tmp = tmp * omega_corr2
        end do
        factor = 2.0_dp*Pi*RhoInv

      case default
        call die('Unknown Hartree-Fock potential type')

    end select

    tmp = (Pi*ZetapEtaInv)**3
    factor = factor*S1234*SQRT(tmp)

    do i=1,m_max+1
      prim%F(i) = prim%F(i)*factor
    end do

    prim%U(:,1)    = P-A
    prim%U(:,2)    = P-B
    prim%U(:,3)    = Q-C
    prim%U(:,4)    = Q-D
    prim%U(:,5)    = W-P
    prim%U(:,6)    = W-Q
    prim%twozeta_a = 2.0_dp*Zeta_A
    prim%twozeta_b = 2.0_dp*Zeta_B
    prim%twozeta_c = 2.0_dp*Zeta_C
    prim%twozeta_d = 2.0_dp*Zeta_D
    prim%oo2z      = 0.5_dp*ZetaInv
    prim%oo2n      = 0.5_dp*EtaInv
    prim%oo2zn     = 0.5_dp*ZetapEtaInv
    prim%poz       = Rho*ZetaInv
    prim%pon       = Rho*EtaInv
    prim%oo2p      = 0.5_dp*RhoInv

  end subroutine build_deriv_data

! *****************************************************************************
!> \brief Fills the data structure used in Libint.
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \par History
!!      - 03.2007 Created [Manuel Guidon]
!!      - 01.2018 Changed interface for SIESTA [Yann Pouillon]
!!
!! \param[in] A: ...
!! \param[in] B: ...
!! \param[in] C: ...
!! \param[in] D: ...
!! \param[in] m_max: ...
!! \param[in] hfx_opts: ...
!! \param[in,out] prim: ...
!! \param[in,out] do_it: ...
!! \param[in] ZetaInv: ...
!! \param[in] EtaInv: ...
!! \param[in] ZetapEtaInv: ...
!! \param[in] Rho: ...
!! \param[in] RhoInv: ...
!! \param[in] S1234: ...
!! \param[in] P: ...
!! \param[in] Q: ...
!! \param[in] W: ...
!! \param[in] PQ2: ...
! *****************************************************************************
  subroutine build_quartet_data(A, B, C, D, m_max, hfx_opts, &
&                prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, &
&                S1234, P, Q, W, PQ2)

    use units,          only: pi
    use nao2gto_common
    use nao2gto_libint, only: prim_data
    use nao2gto_utils,  only: fgamma => fgamma_0
    use nao2gto_types,  only: hfx_options_type

    implicit none

    ! Arguments
    real(dp), intent(in)               :: A(3), B(3), C(3), D(3)
    integer, intent(in)                :: m_max
    type(prim_data), intent(inout)     :: prim
    type(hfx_options_type), intent(in) :: hfx_opts
    logical, intent(inout)             :: do_it
    real(dp), intent(in)               :: ZetaInv, EtaInv, ZetapEtaInv, &
&                                         Rho, RhoInv, S1234, P(3), Q(3), &
&                                         W(3), PQ2

    ! Local variables
    integer                 :: i
    real(dp)                :: factor, T, tmp
    real(dp)                :: omega2, omega_corr, omega_corr2
    real(dp), dimension(17) :: Fm

    ! -------------------------------------------------------------------------

    T = Rho*PQ2
    do_it = .true.

    select case (hfx_opts%potential_type)

      case(do_hfx_potential_coulomb)
        call fgamma(m_max, T, prim%F(1))
        factor = 2.0_dp*Pi*RhoInv

      case(do_hfx_potential_short)
        call fgamma(m_max, T, prim%F(1))
        omega2 = hfx_opts%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        call fgamma(m_max,T,Fm)
        tmp = - omega_corr
        do i=1,m_max+1
          prim%F(i) = prim%F(i) + Fm(i)*tmp
          tmp = tmp * omega_corr2
        end do
        factor = 2.0_dp*Pi*RhoInv

    end select

    tmp    = (Pi*ZetapEtaInv)**3
    factor = factor*S1234*SQRT(tmp)

    do i=1,m_max+1
      prim%F(i) = prim%F(i)*factor
    end do

    prim%U(1:3,1) = P-A
    prim%U(1:3,3) = Q-C
    prim%U(1:3,5) = W-P
    prim%U(1:3,6) = W-Q
    prim%oo2z     = 0.5_dp*ZetaInv
    prim%oo2n     = 0.5_dp*EtaInv
    prim%oo2zn    = 0.5_dp*ZetapEtaInv
    prim%poz      = Rho*ZetaInv
    prim%pon      = Rho*EtaInv
    prim%oo2p     = 0.5_dp*RhoInv

  end subroutine build_quartet_data

! *****************************************************************************
  SUBROUTINE build_quartet_data_screen(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max, &
                                       hfx_opts, prim, do_it, np, nq,R11,R22)

    use units, only: pi
    use nao2gto_libint
    use nao2gto_types
    use nao2gto_utils, only: fgamma => fgamma_0

    implicit none

    REAL(KIND=dp)                            :: A(3), B(3), C(3), D(3)
    REAL(KIND=dp), INTENT(IN)                :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: m_max
    TYPE(hfx_options_type)                :: hfx_opts
    TYPE(prim_data)                          :: prim
    LOGICAL, INTENT(INOUT)                   :: do_it
    INTEGER                                  :: np, nq
    REAL(dp)                                 :: R11, R22

    INTEGER                                  :: i
    LOGICAL                                  :: use_gamma
    REAL(KIND=dp) :: AB(3), AB2, CD(3), CD2, Eta, EtaInv, factor, omega2, &
      omega_corr, omega_corr2, P(3), PQ(3), PQ2, Q(3), R, R1, R2, Rho, &
      RhoInv, S1234, T, tmp, W(3), Zeta, ZetaInv, ZetapEtaInv
    REAL(KIND=dp), DIMENSION(17)             :: Fm
    Zeta = Zeta_A + Zeta_B
    ZetaInv = 1.0_dp/Zeta
    Eta  = Zeta_C + Zeta_D
    EtaInv = 1.0_dp/Eta
    ZetapEtaInv = Zeta+Eta
    ZetapEtaInv = 1.0_dp/ZetapEtaInv
    Rho  = Zeta*Eta*ZetapEtaInv
    RhoInv = 1.0_dp/Rho

    DO i=1,3
      P(i) = (Zeta_A*A(i) + Zeta_B*B(i))*ZetaInv
      Q(i) = (Zeta_C*C(i) + Zeta_D*D(i))*EtaInv
      AB(i) = A(i)-B(i)
      CD(i) = C(i)-D(i)
      PQ(i) = P(i)-Q(i)
      W(i) = (Zeta*P(i) + Eta*Q(i))*ZetapEtaInv
    END DO

    AB2 = DOT_PRODUCT(AB,AB)
    CD2 = DOT_PRODUCT(CD,CD)
    PQ2 = DOT_PRODUCT(PQ,PQ)

    S1234= EXP((-Zeta_A*Zeta_B*ZetaInv*AB2)+(-Zeta_C*Zeta_D*EtaInv*CD2))
    T = Rho*PQ2

    do_it = .TRUE.

    SELECT CASE(hfx_opts%potential_type)
      CASE(do_hfx_potential_coulomb)
        CALL fgamma(m_max,T,prim%F)
        factor = 2.0_dp*Pi*RhoInv
      CASE(do_hfx_potential_short)
         R = hfx_opts%cutoff_radius*SQRT(rho)
         R1 = R11
         R2 = R22
    !     IF (PQ2 > (R1+R2+hfx_opts%cutoff_radius)**2) THEN
    !        RETURN
    !     END IF
        CALL fgamma(m_max,T,prim%F)
        omega2 = hfx_opts%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,Fm)
        tmp = - omega_corr
        DO i=1,m_max+1
          prim%F(i)=prim%F(i) + Fm(i)*tmp
          tmp = tmp * omega_corr2
        END DO
        factor = 2.0_dp*Pi*RhoInv
    END SELECT

    tmp    = (Pi*ZetapEtaInv)**3
    factor = factor*S1234*SQRT(tmp)

    DO i=1,m_max+1
       prim%F(i)=prim%F(i)*factor
    ENDDO
    prim%U(1:3,1) = P-A
    prim%U(1:3,3) = Q-C
    prim%U(1:3,5) = W-P
    prim%U(1:3,6) = W-Q
    prim%oo2z      = 0.5_dp*ZetaInv
    prim%oo2n      = 0.5_dp*EtaInv
    prim%oo2zn     = 0.5_dp*ZetapEtaInv
    prim%poz       = Rho*ZetaInv
    prim%pon       = Rho*EtaInv
    prim%oo2p      = 0.5_dp*RhoInv
  END SUBROUTINE build_quartet_data_screen



end module nao2gto_primitive
