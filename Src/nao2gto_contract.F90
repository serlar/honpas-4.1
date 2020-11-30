! *** Module: nao2gto_contract ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief NAO2GTO contraction routines
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.2010 Edited [Xinming Qin]
!!      - 01.2016 Edited [Honghui Shang]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_contract

  use nao2gto_common

  implicit none

  private

  real(dp), parameter :: erfc_epsilon = 0.99998871620832942100_dp

  public :: calc_contract_eri, calc_contract_deriv_eri, calc_contract_eri2

contains

! *****************************************************************************
!> \brief Computes contractions of spherical Gaussians
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] libint_data: Libint data structure
!! \param[in] cell: supercell in real space
!! \param[in] rcell: supercell in reciprocal space
!! \param[in] ra: ...
!! \param[in] rb: ...
!! \param[in] rc: ...
!! \param[in] rd: ...
!! \param[in] npgfa: ...
!! \param[in] npgfb: ...
!! \param[in] npgfc: ...
!! \param[in] npgfd: ...
!! \param[in] la: ...
!! \param[in] lb: ...
!! \param[in] lc: ...
!! \param[in] ld: ...
!! \param[in] ncoa: ...
!! \param[in] ncob: ...
!! \param[in] ncoc: ...
!! \param[in] ncod: ...
!! \param[in] zeta: ...
!! \param[in] zetb: ...
!! \param[in] zetc: ...
!! \param[in] zetd: ...
!! \param[in] sphia: ...
!! \param[in] sphib: ...
!! \param[in] sphic: ...
!! \param[in] sphid: ...
!! \param[in] hfx_opts: data structure containing Hartree-Fok exchange
!!                           parameters
!! \param[out] eri: ...
! *****************************************************************************
  subroutine calc_contract_eri(libint_data, cell, rcell, ra, rb, rc, rd,  &
&                npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, &
&                ncoa, ncob, ncoc, ncod, zeta, zetb, zetc, zetd, &
&                sphia, sphib, sphic, sphid, neris_tmp, max_contraction, &
&                max_val2_set, log10_pmax, R1_pgf, R2_pgf, pgf1, &
&                pgf2, hfx_opts, eri)

    use units,              only: pi
    use atm_types,          only: nso, nco
    use alloc,              only: de_alloc, re_alloc
    use nao2gto_data
    use nao2gto_libint,     only: Libint_t
    use nao2gto_pbc,        only: trans_pbc
    use nao2gto_primitive,  only: calc_primitive_eri
    use nao2gto_types

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data
    integer, intent(in)   :: npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, &
&                            ncoa, ncob, ncoc, ncod
    integer(int_8), intent(out) :: neris_tmp
    real(dp), intent(in)  :: cell(3,3), rcell(3,3), &
&     ra(3), rb(3), rc(3), rd(3), &
&     zeta(npgfa), zetb(npgfb), zetc(npgfc), zetd(npgfd), &
&     sphia(ncoa,nso(la)), sphib(ncob,nso(lb)), &
&     sphic(ncoc,nso(lc)), sphid(ncod,nso(ld))
    real(dp), intent(in) :: max_contraction, max_val2_set, log10_pmax
    real(dp), intent(out) :: eri(nso(la),nso(lb),nso(lc),nso(ld))
    type(hfx_screen_coeff_type), dimension(:,:), intent(in) :: &
&     R1_pgf, R2_pgf, pgf1, pgf2
    type(hfx_options_type), intent(in) :: hfx_opts

    ! Local variables
    integer :: ipgf, jpgf, kpgf, lpgf, offset_a, offset_b, &
&              offset_c, offset_d, am, bm, cm, dm, i, &
&              index_primitive_integrals
    integer :: ieri, j, k, l
    real(dp) :: Eta, EtaInv, P(3), Q(3), Rho, RhoInv, rpq2, S1234, S1234a, &
&               tmp_max, W(3), Zeta1, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
&               ZetaInv, ZetapEtaInv, rab2, rcd2, r_temp(3), r_pbc_temp(3), &
&               alpha_P, alpha_Q, R_P, R_Q, Kab, Kcd, theta_w, &
&               rc_trans(3), rd_trans(3), rab(3), rcd(3), R1, R2
    real(dp) :: far_eri, pgf_max_1, pgf_max_2, cart_estimate

    real(dp), dimension(:), pointer :: primitive_integrals => null()
    real(dp), dimension(:), pointer :: T1 => null()

    external :: dgemm

    ! -------------------------------------------------------------------------

    call re_alloc(primitive_integrals, 1, ncoa*ncob*ncoc*ncod, &
&     name='primitive_integrals', routine='calc_contract_eri')
    call re_alloc(T1, 1, ncoa*ncob*ncoc*ncod, &
&     name='T1', routine='calc_contract_eri')
    primitive_integrals(:) = 0.0_dp
    T1(:) = 0.0_dp

    neris_tmp = 0
    cart_estimate = 0.0_dp
    eri(:,:,:,:) = 0.0_dp

    rab(1:3) = ra(1:3) - rb(1:3)
    rcd(1:3) = rc(1:3) - rd(1:3)
    rab2 = (rab(1)**2) + (rab(2)**2) + (rab(3)**2)
    rcd2 = (rcd(1)**2) + (rcd(2)**2) + (rcd(3)**2)

    do ipgf = 1,npgfa
      offset_a = (ipgf-1)*nco(la)
      Zeta_A = zeta(ipgf)

      do jpgf = 1,npgfb
        offset_b = (jpgf-1)*nco(lb)
        pgf_max_1 = pgf1(jpgf,ipgf)%x(1)*rab2+pgf1(jpgf,ipgf)%x(2)
        if ( pgf_max_1+max_val2_set+log10_pmax < log10_eps_schwarz ) cycle

        Zeta_B = zetb(jpgf)
        Zeta1 = Zeta_A + Zeta_B
        ZetaInv = 1.0_dp/Zeta1
        S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
        P(1:3) = (Zeta_A*ra(1:3) + Zeta_B*rb(1:3))*ZetaInv

        do kpgf = 1,npgfc
          offset_c = (kpgf-1)*nco(lc)
          Zeta_C = zetc(kpgf)

          do lpgf = 1,npgfd
            offset_d = (lpgf - 1)*nco(ld)
            pgf_max_2 = pgf2(lpgf,kpgf)%x(1)*rcd2+pgf2(lpgf,kpgf)%x(2)
            if ( pgf_max_1+pgf_max_2+log10_pmax < log10_eps_schwarz ) cycle

            Zeta_D = zetd(lpgf)
            Eta = Zeta_C + Zeta_D
            EtaInv = 1.0_dp/Eta
            ZetapEtaInv = Zeta1+Eta
            ZetapEtaInv = 1.0_dp/ZetapEtaInv
            Rho = Zeta1*Eta*ZetapEtaInv
            RhoInv = 1.0_dp/Rho
            S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
            Q(1:3) = (Zeta_C*rc(1:3) + Zeta_D*rd(1:3))*EtaInv
            r_temp(1:3) = Q(1:3) - P(1:3)
            call trans_pbc(r_temp, cell,rcell, r_pbc_temp)

            Q(1:3) = P(1:3) + r_pbc_temp(1:3)
            rc_trans(1:3) = rc(1:3) + r_pbc_temp(1:3) - r_temp(1:3)
            rd_trans(1:3) = rd(1:3) + r_pbc_temp(1:3) - r_temp(1:3)

            rpq2 = ((P(1)-Q(1))**2) + ((P(2)-Q(2))**2) + ((P(3)-Q(3))**2)
            W(1:3) = (Zeta1*P(1:3) + Eta*Q(1:3))*ZetapEtaInv

            tmp_max = 0.0_dp

            if ( .not. hfx_opts%far_field ) then

              call calc_primitive_eri(libint_data, ra, rb, rc_trans, rd_trans, &
&               la, lb, lc ,ld, ncoa, ncob, ncoc, ncod, &
&               offset_a, offset_b, offset_c, offset_d, &
&               primitive_integrals, hfx_opts, &
&               max_contraction, tmp_max, neris_tmp, &
&               ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
&               P, Q, W, rpq2, rab, rcd, R1, R2)

            else ! farfield screening for GTO ERIs

              alpha_P = Zeta_A*Zeta_B*ZetaInv
              R_P = aint(1.0_dp/(sqrt(2.0_dp*alpha_P)*erfc_epsilon)) + 1.0_dp
              alpha_Q = Zeta_C*Zeta_D*EtaInv
              R_Q = aint(1.0_dp/(sqrt(2.0_dp*alpha_Q)*erfc_epsilon)) + 1.0_dp

              if ( sqrt(rpq2) .gt. (R_P + R_Q) ) then

                Kab = sqrt(2.0_dp)*pi**1.25_dp*ZetaInv*exp(-alpha_P*rab2)
                Kcd = sqrt(2.0_dp)*pi**1.25_dp*EtaInv*exp(-alpha_Q*rcd2)
                theta_w = 1.0_dp/ &
&                 (1.0_dp/alpha_P+1.0_dp/alpha_Q+82.644628099173553719)

                far_eri = max_contraction * &
&                 (Kab*Kcd*erfc((theta_w**0.5_dp)*sqrt(rpq2))/sqrt(rpq2))

                if ( far_eri .gt. hfx_opts%eps_far ) then

                  call calc_primitive_eri(libint_data, ra, rb, rc_trans, rd_trans, &
&                   la, lb, lc ,ld, ncoa, ncob, ncoc, ncod, &
&                   offset_a, offset_b, offset_c, offset_d, &
&                   primitive_integrals, hfx_opts, &
&                   max_contraction, tmp_max, neris_tmp, &
&                   ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
&                   P, Q, W, rpq2, rab, rcd, R1, R2)

                endif

              else ! near field

                call calc_primitive_eri(libint_data, ra, rb, rc_trans, rd_trans, &
&                 la, lb, lc ,ld, ncoa, ncob, ncoc, ncod, &
&                 offset_a, offset_b, offset_c, offset_d, &
&                 primitive_integrals, hfx_opts, &
&                 max_contraction, tmp_max, neris_tmp, &
&                 ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
&                 P, Q, W, rpq2, rab, rcd, R1, R2)

              endif ! sqrt(rpq2) .gt. (R_P + R_Q)

            endif ! .not. hfx_opts%far_field

            cart_estimate = max(tmp_max, cart_estimate)

          enddo ! lpgf
        enddo ! kpgf
      enddo ! jpgf
    enddo ! ipgf

    if ( cart_estimate >= hfx_opts%eps_schwarz ) then

      call dgemm("T", "N", ncob*ncoc*ncod, nso(la), ncoa, &
&       1.0_dp, primitive_integrals(1), ncoa, &
&       sphia(1,1), ncoa, 0.0_dp, T1(1), ncob*ncoc*ncod)

      call dgemm("T", "N", nso(la)*ncoc*ncod, nso(lb), ncob, &
&       1.0_dp, T1(1), ncob, sphib(1,1), ncob, &
&       0.0_dp, primitive_integrals(1), nso(la)*ncoc*ncod)

      call dgemm("T", "N", nso(la)*nso(lb)*ncod, nso(lc), ncoc, &
&       1.0_dp, primitive_integrals(1), ncoc, &
&       sphic(1,1), ncoc, 0.0_dp, T1(1), nso(la)*nso(lb)*ncod)

      call dgemm("T", "N", nso(la)*nso(lb)*nso(lc), nso(ld), ncod, &
&       1.0_dp, T1(1), ncod, sphid(1,1), ncod, &
&       0.0_dp, primitive_integrals(1), nso(la)*nso(lb)*nso(lc))

    end if

    do dm=1,nso(ld)
      do cm=1,nso(lc)
        do bm=1,nso(lb)
          do am=1,nso(la)
            index_primitive_integrals = &
&             (dm-1)*nso(lc)*nso(lb)*nso(la) + (cm-1)*nso(lb)*nso(la) + &
&             (bm-1)*nso(la) + am
            eri(am,bm,cm,dm) = primitive_integrals(index_primitive_integrals)
          enddo
        enddo
      enddo
    enddo

    call de_alloc(primitive_integrals, &
&     name='primitive_integrals', routine='calc_contract_eri')
    call de_alloc(T1, &
&     name='primitive_integrals', routine='calc_contract_eri')

  end subroutine calc_contract_eri

! *****************************************************************************
!> \brief Computes contracted derivatives of spherical Gaussians
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] deriv_data: Libderiv data structure
!! \param[in] cell: supercell in real space
!! \param[in] rcell: supercell in reciprocal space
!! \param[in] ra: ...
!! \param[in] rb: ...
!! \param[in] rc: ...
!! \param[in] rd: ...
!! \param[in] npgfa: ...
!! \param[in] npgfb: ...
!! \param[in] npgfc: ...
!! \param[in] npgfd: ...
!! \param[in] la: ...
!! \param[in] lb: ...
!! \param[in] lc: ...
!! \param[in] ld: ...
!! \param[in] ncoa: ...
!! \param[in] ncob: ...
!! \param[in] ncoc: ...
!! \param[in] ncod: ...
!! \param[in] zeta: ...
!! \param[in] zetb: ...
!! \param[in] zetc: ...
!! \param[in] zetd: ...
!! \param[in] sphia: ...
!! \param[in] sphib: ...
!! \param[in] sphic: ...
!! \param[in] sphid: ...
!! \param[in] hfx_opts: data structure containing Hartree-Fok exchange
!!                           parameters
!! \param[out] eri_force: ...
! *****************************************************************************
  subroutine calc_contract_deriv_eri(deriv_data, cell, rcell, ra, rb, rc, rd, &
&                npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, &
&                ncoa, ncob, ncoc, ncod, zeta, zetb, zetc, zetd, &
&                sphia, sphib, sphic, sphid, hfx_opts, eri_force, &
&                max_contraction, max_val2_set, log10_pmax, &
&                R1_pgf, R2_pgf, pgf1, pgf2)

    use units,             only: pi
    use alloc,             only: de_alloc, re_alloc
    use atm_types,         only: nso, nco
    use nao2gto_data
    use nao2gto_libint,    only: Libderiv_t
    use nao2gto_pbc,       only: cell_pbc
    use nao2gto_primitive, only: calc_primitive_deriv_eri
    use nao2gto_types

    implicit none

    ! Arguments
    type(Libderiv_t), intent(inout) :: deriv_data
    integer, intent(in)  :: npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, &
&                           ncoa, ncob, ncoc, ncod
    real(dp), intent(in)  :: cell(3,3), rcell(3,3), ra(3), rb(3), &
&     rc(3), rd(3), zeta(npgfa), zetb(npgfb), zetc(npgfc), zetd(npgfd), &
&     sphia(ncoa,nso(la)), sphib(ncob,nso(lb)), &
&     sphic(ncoc,nso(lc)), sphid(ncod,nso(ld))
    type(hfx_options_type), intent(in) :: hfx_opts
    real(dp), intent(out) :: eri_force(nso(la),nso(lb),nso(lc),nso(ld),12)
    real(dp), intent(in) :: max_contraction, max_val2_set, log10_pmax
    type(hfx_screen_coeff_type), dimension(:, :), pointer :: &
&     R1_pgf, R2_pgf, pgf1, pgf2

    ! Local variables
    integer  :: ipgf, jpgf, kpgf, lpgf, offset_a, offset_b, &
&     offset_c,offset_d, am, bm, cm, dm, i, index_primitive_force, coord
    real(dp) :: Eta, EtaInv, P(3), Q(3), Rho, &
&     RhoInv, rpq2, S1234, S1234a, tmp_max, W(3), Zeta1, Zeta_A, Zeta_B, &
&     Zeta_C, Zeta_D, ZetaInv, ZetapEtaInv, rab(3), rcd(3), rab2, rcd2, &
&     r_temp(3), r_pbc_temp(3), alpha_P, alpha_Q, R_P, R_Q, Kab, Kcd, &
&     theta_w, rc_trans(3), rd_trans(3)
    real(dp) :: far_eri, R1, R2, pgf_max_1, pgf_max_2, cart_estimate

    real(dp), dimension(:), pointer :: primitive_force => null()
    real(dp), dimension(:), pointer :: work_forces => null()
    real(dp), dimension(:), pointer :: T1 => null(), T2 => null()

    external :: dgemm

    ! -------------------------------------------------------------------------

    call re_alloc(primitive_force, 1, ncoa*ncob*ncoc*ncod*12, &
&     name='primitive_force', routine='calc_contract_deriv_eri')
    primitive_force(:) = 0.0_dp
    call re_alloc(T1, 1, ncoa*ncob*ncoc*ncod, &
&     name='work_forces', routine='calc_contract_deriv_eri')
    T1(:) = 0.0_dp
    call re_alloc(work_forces, 1, nco(la)*nco(lb)*nco(lc)*nco(ld)*12, &
&     name='work_forces', routine='calc_contract_deriv_eri')

    cart_estimate = 0.0_dp
    rab(1:3) = ra(1:3) - rb(1:3)
    rcd(1:3) = rc(1:3) - rd(1:3)
    rab2 = rab(1)**2 + rab(2)**2 + rab(3)**2
    rcd2 = rcd(1)**2 + rcd(2)**2 + rcd(3)**2

    do ipgf = 1,npgfa
      offset_a = (ipgf-1)*nco(la)
      Zeta_A = zeta(ipgf)

      do jpgf = 1,npgfb
        offset_b = (jpgf-1)*nco(lb)
        pgf_max_1 = pgf1(jpgf,ipgf)%x(1)*rab2+pgf1(jpgf,ipgf)%x(2)
        if ( pgf_max_1+max_val2_set+log10_pmax < log10_eps_schwarz ) cycle

        Zeta_B = zetb(jpgf)
        Zeta1 = Zeta_A + Zeta_B
        ZetaInv = 1.0_dp/Zeta1
        S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
        P(1:3) = (Zeta_A*ra(1:3) + Zeta_B*rb(1:3))*ZetaInv

        do kpgf = 1,npgfc
          offset_c = (kpgf-1)*nco(lc)
          Zeta_C = zetc(kpgf)

          do lpgf = 1,npgfd
            offset_d = (lpgf-1)*nco(ld)
            pgf_max_2 = pgf2(lpgf,kpgf)%x(1)*rcd2+pgf2(lpgf,kpgf)%x(2)
            if (pgf_max_1+pgf_max_2+log10_pmax < log10_eps_schwarz ) cycle

            Zeta_D = zetd(lpgf)
            Eta = Zeta_C + Zeta_D
            EtaInv = 1.0_dp/Eta
            ZetapEtaInv = Zeta1+Eta
            ZetapEtaInv = 1.0_dp/ZetapEtaInv
            Rho = Zeta1*Eta*ZetapEtaInv
            RhoInv = 1.0_dp/Rho
            S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
            Q(1:3) = (Zeta_C*rc(1:3) + Zeta_D*rd(1:3))*EtaInv
            r_temp(1:3) = Q(1:3)-P(1:3)
            r_pbc_temp = cell_pbc(r_temp,cell,rcell)

            Q(1:3) = P(1:3) + r_pbc_temp(1:3)
            rc_trans(1:3) = rc(1:3) + r_pbc_temp(1:3) - r_temp(1:3)
            rd_trans(1:3) = rd(1:3) + r_pbc_temp(1:3) - r_temp(1:3)

            rpq2 = (P(1)-Q(1))**2+(P(2)-Q(2))**2+(P(3)-Q(3))**2
            W(1:3) = (Zeta1*P(1:3)+Eta*Q(1:3))*ZetapEtaInv

            tmp_max = 0.0_dp

            if ( .not. hfx_opts%far_field ) then

              call calc_primitive_deriv_eri(deriv_data, ra, rb, &
&               rc_trans, rd_trans, zeta(ipgf), zetb(jpgf), &
&               zetc(kpgf), zetd(lpgf), la, lb, lc ,ld, work_forces, &
&               ncoa, ncob, ncoc, ncod, primitive_force, hfx_opts, &
&               max_contraction, tmp_max, &
&               offset_a, offset_b, offset_c, offset_d, ZetaInv, EtaInv, &
&               ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, rpq2, rab, rcd)

            else

              alpha_P = Zeta_A*Zeta_B*ZetaInv
              R_P = aint(1.0_dp/(sqrt(2.0_dp*alpha_P)*erfc_epsilon)) + 1.0_dp
              alpha_Q = Zeta_C*Zeta_D*EtaInv
              R_Q = aint(1.0_dp/(sqrt(2.0_dp*alpha_Q)*erfc_epsilon)) + 1.0_dp

              if ( sqrt(rpq2) > (R_P + R_Q) ) then ! far field

                Kab = sqrt(2.0_dp)*pi**1.25_dp*ZetaInv*exp(-alpha_P*rab2)
                Kcd = sqrt(2.0_dp)*pi**1.25_dp*EtaInv*exp(-alpha_Q*rcd2)

                ! Note: 1.0_dp/(0.11*0.11)
                theta_w = 1.0_dp/ &
&                 (1.0_dp/alpha_P+1.0_dp/alpha_Q+82.644628099173553719)

                far_eri = max_contraction * &
&                 (Kab*Kcd*erfc((theta_w**0.5_dp)*sqrt(rpq2))/sqrt(rpq2))

                if ( far_eri > hfx_opts%eps_far ) then

                  call calc_primitive_deriv_eri(deriv_data, ra, rb, &
&                   rc_trans, rd_trans, zeta(ipgf), zetb(jpgf), &
&                   zetc(kpgf), zetd(lpgf), la, lb, lc ,ld, work_forces, &
&                   ncoa, ncob, ncoc, ncod, primitive_force, hfx_opts, &
&                   max_contraction, tmp_max, &
&                   offset_a, offset_b, offset_c, offset_d, ZetaInv, &
&                   EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, &
&                   rpq2, rab, rcd)

                endif

              else !near field

                call calc_primitive_deriv_eri(deriv_data, ra, rb, &
&                 rc_trans, rd_trans, zeta(ipgf), zetb(jpgf), &
&                 zetc(kpgf), zetd(lpgf), la, lb, lc ,ld, work_forces, &
&                 ncoa, ncob, ncoc, ncod, primitive_force, hfx_opts, &
&                 max_contraction, tmp_max, &
&                 offset_a, offset_b, offset_c, offset_d, ZetaInv, &
&                 EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, &
&                 rpq2, rab, rcd)

              endif ! sqrt(rpq2) .gt. (R_P + R_Q)

            endif ! .not. hfx_opts%far_field

            cart_estimate = max(tmp_max, cart_estimate)

          enddo ! lpgf
        enddo ! kpgf
      enddo ! jpgf
    enddo ! ipgf

    if ( cart_estimate < hfx_opts%eps_schwarz ) then

      do coord = 1,12

        T2 => primitive_force( &
&         (coord-1)*ncoa*ncob*ncoc*ncod+1:coord*ncoa*ncob*ncoc*ncod)

        call dgemm("T", "N", ncob*ncoc*ncod, nso(la), ncoa, &
&         1.0_dp, T2(1), ncoa, sphia(1,1), ncoa, &
&         0.0_dp, T1(1), ncob*ncoc*ncod)

        call dgemm("T", "N", nso(la)*ncoc*ncod, nso(lb), ncob, &
&         1.0_dp, T1(1), ncob, sphib(1,1), ncob, &
&         0.0_dp, T2(1), nso(la)*ncoc*ncod)

        call dgemm("T", "N", nso(la)*nso(lb)*ncod, nso(lc), ncoc, &
&         1.0_dp, T2(1), ncoc, sphic(1,1), ncoc, &
&         0.0_dp, T1(1), nso(la)*nso(lb)*ncod)

        call dgemm("T", "N", nso(la)*nso(lb)*nso(lc), nso(ld), ncod, &
&         1.0_dp, T1(1), ncod, sphid(1,1), ncod, &
&         0.0_dp, T2(1), nso(la)*nso(lb)*nso(lc))

        do dm=1,nso(ld)
          do cm=1,nso(lc)
            do bm=1,nso(lb)
              do am=1,nso(la)
                index_primitive_force = (dm-1)*nso(lc)*nso(lb)*nso(la) + &
&                                       (cm-1)*nso(lb)*nso(la) + &
&                                       (bm-1)*nso(la) + am
                eri_force(am,bm,cm,dm,coord) = T2(index_primitive_force)
              enddo ! am
            enddo ! bm
          enddo ! cm
        enddo ! dm

      enddo ! coord

    end if

    nullify(T2)
    call de_alloc(primitive_force, &
&     name='primitive_force', routine='calc_contract_deriv_eri')
    call de_alloc(work_forces, &
&     name='work_forces', routine='calc_contract_deriv_eri')
    call de_alloc(T1, &
&     name='T1', routine='calc_contract_deriv_eri')

  end subroutine calc_contract_deriv_eri

! *****************************************************************************
!> \brief Computes contractions of spherical Gaussians
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] libint_data: Libint data structure
!! \param[in] cell: supercell in real space
!! \param[in] rcell: supercell in reciprocal space
!! \param[in] ra: ...
!! \param[in] rb: ...
!! \param[in] rc: ...
!! \param[in] rd: ...
!! \param[in] npgfa: ...
!! \param[in] npgfb: ...
!! \param[in] npgfc: ...
!! \param[in] npgfd: ...
!! \param[in] la: ...
!! \param[in] lb: ...
!! \param[in] lc: ...
!! \param[in] ld: ...
!! \param[in] ncoa: ...
!! \param[in] ncob: ...
!! \param[in] ncoc: ...
!! \param[in] ncod: ...
!! \param[in] zeta: ...
!! \param[in] zetb: ...
!! \param[in] zetc: ...
!! \param[in] zetd: ...
!! \param[in] sphia: ...
!! \param[in] sphib: ...
!! \param[in] sphic: ...
!! \param[in] sphid: ...
!! \param[in] hfx_opts: data structure containing Hartree-Fok exchange
!!                           parameters
!! \param[out] eri: ...
! *****************************************************************************
  subroutine calc_contract_eri2(lib, cell, rcell, ra, rb, rc, rd, &
      npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, ncoa, ncob, ncoc, ncod, &
      zeta, zetb, zetc, zetd, sphia, sphib, sphic, sphid, hfx_opts, eri)

    use alloc
    use atm_types, only: nso, nco
    use nao2gto_libint
    use nao2gto_pbc
    use nao2gto_primitive
    use nao2gto_types

    ! FIXME: for debugging
    use parallel

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: lib
    integer, intent(in) :: npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, &
&     ncoa, ncob, ncoc, ncod
    real(dp) :: cell(3,3), rcell(3,3), ra(3), rb(3), rc(3), rd(3), &
&     rab(3), rcd(3), zeta(npgfa), zetb(npgfb), zetc(npgfc), zetd(npgfd), &
&     sphia(ncoa,nso(la)), sphib(ncob,nso(lb)), sphic(ncoc,nso(lc)), &
&     sphid(ncod,nso(ld))
    type(hfx_options_type) :: hfx_opts
    real(dp) :: eri(nso(la),nso(lb),nso(lc),nso(ld))

    ! Local variables
    real(dp), parameter ::  pi=3.14159265358979323846264338_dp
    real(dp), parameter ::  erfc_epsilon = 0.99998871620832942100_dp
    integer :: ipgf, jpgf, kpgf, lpgf, offset_a, offset_b, offset_c, &
&     offset_d, am, bm, cm, dm, ieri, i, j, k, l, index_primitive_integrals
    real(dp) :: Eta, EtaInv, P(3), Q(3), Rho, RhoInv, rpq2, S1234, S1234a, &
&     tmp_max, W(3), Zeta1, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
&     ZetaInv, ZetapEtaInv, rab2, rcd2, r_temp(3), r_pbc_temp(3), &
&     alpha_P, alpha_Q, R_P, R_Q, Kab, Kcd, theta_w, rc_trans(3), rd_trans(3)
    real(dp), pointer :: primitive_integrals(:) => null(), T1(:) => null()

    ! -------------------------------------------------------------------------

    call re_alloc(primitive_integrals, 1, ncoa*ncob*ncoc*ncod, &
&     name='primitive_integrals', routine='calc_contract_eri2')
    primitive_integrals(:) = 0.0_dp
    call re_alloc(T1, 1, ncoa*ncob*ncoc*ncod, &
&     name='T1', routine='calc_contract_eri2')
    T1(:) = 0.0_dp

    rab(1:3) = ra(1:3) - rb(1:3)
    rcd(1:3) = rc(1:3) - rd(1:3)
    rab2 = rab(1)**2 + rab(2)**2 + rab(3)**2
    rcd2 = rcd(1)**2 + rcd(2)**2 + rcd(3)**2

    do ipgf=1,npgfa
      offset_a = (ipgf-1)*nco(la)
      Zeta_A = zeta(ipgf)

      do jpgf=1,npgfb
        offset_b = (jpgf-1)*nco(lb)
        Zeta_B = zetb(jpgf)
        Zeta1 = Zeta_A + Zeta_B
        ZetaInv = 1.0d0/Zeta1
        S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
        P(1:3) = (Zeta_A*ra(1:3) + Zeta_B*rb(1:3))*ZetaInv

        do kpgf=1,npgfc
          offset_c = (kpgf-1)*nco(lc)
          Zeta_C = zetc(kpgf)

          do lpgf=1,npgfd
            offset_d = (lpgf-1)*nco(ld)
            Zeta_D = zetd(lpgf)
            Eta = Zeta_C + Zeta_D
            EtaInv = 1.0d0/Eta
            ZetapEtaInv = Zeta1+Eta
            ZetapEtaInv = 1.0d0/ZetapEtaInv
            Rho = Zeta1*Eta*ZetapEtaInv
            RhoInv = 1.0d0/Rho
            S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
            Q(1:3) = (Zeta_C*rc(1:3) + Zeta_D*rd(1:3))*EtaInv
            r_temp(1:3) = Q(1:3)-P(1:3)
            call trans_pbc(r_temp, cell, rcell, r_pbc_temp)

            Q(1:3) = P(1:3)+r_pbc_temp(1:3)
            rc_trans(1:3) = rc(1:3) + r_pbc_temp(1:3) - r_temp(1:3)
            rd_trans(1:3) = rd(1:3) + r_pbc_temp(1:3) - r_temp(1:3)
            rpq2 = (P(1)-Q(1))**2 + (P(2)-Q(2))**2 + (P(3)-Q(3))**2
            W(1:3) = (Zeta1*P(1:3) + Eta*Q(1:3))*ZetapEtaInv

            call calc_primitive_eri2(lib, ra, rb, rc_trans, rd_trans, &
&             zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf), la, lb, lc ,ld, &
&             ncoa, ncob, ncoc, ncod, offset_a, offset_b, offset_c, offset_d, &
&             primitive_integrals, hfx_opts, Zeta1, ZetaInv, Eta, EtaInv, &
&             ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, rpq2, rab, rcd)

          enddo
        enddo
      enddo
    enddo

    call dgemm("T", "N", ncob*ncoc*ncod, nso(la), ncoa, &
      1.0_dp, primitive_integrals(1), ncoa, sphia(1,1), ncoa, &
      0.0_dp, T1(1), ncob*ncoc*ncod)

    call dgemm("T", "N", nso(la)*ncoc*ncod, nso(lb), ncob, &
      1.0_dp, T1(1), ncob, sphib(1,1), ncob, &
      0.0_dp, primitive_integrals(1), nso(la)*ncoc*ncod)

    call dgemm("T", "N", nso(la)*nso(lb)*ncod, nso(lc), ncoc, &
      1.0_dp, primitive_integrals(1), ncoc, sphic(1,1), ncoc, &
      0.0_dp, T1(1), nso(la)*nso(lb)*ncod)

    call dgemm("T", "N", nso(la)*nso(lb)*nso(lc), nso(ld), ncod, &
      1.0_dp, T1(1), ncod, sphid(1,1), ncod, &
      0.0_dp, primitive_integrals(1), nso(la)*nso(lb)*nso(lc))

    do dm=1,nso(ld)
      do cm=1,nso(lc)
        do bm=1,nso(lb)
          do am=1,nso(la)
            index_primitive_integrals = (dm-1)*nso(lc)*nso(lb)*nso(la) + &
              (cm-1)*nso(lb)*nso(la) + (bm-1)*nso(la) + am
            eri(am,bm,cm,dm) = primitive_integrals(index_primitive_integrals)
          enddo
        enddo
      enddo
    enddo

    call de_alloc(primitive_integrals, name='primitive_integrals', &
      routine='calc_contract_eri2')
    call de_alloc(T1, name='T1', routine='calc_contract_eri2')

  end subroutine calc_contract_eri2

end module nao2gto_contract
