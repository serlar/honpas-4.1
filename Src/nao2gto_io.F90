! *** Module: nao2gto_io ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief I/O module for Gaussian-based Hartree-Fock exchange
!!
!!  This module bridges the input file of SIESTA with the NAO2GTO routines,
!!  which calculate the Hartree-Fock exchange interaction using Gaussians.
!!
!! \note
!!      This file currently works with a version of Libint configured for
!!      LIBINT_MAX_AM=5 and LIBINT_MAX_AM1=4 only.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \copyright
!!      - 2010-2018 SIESTA Developers Group
!!
!! \par History
!!      - 11.2017 Reviewed for inclusion in SIESTA [Xinming Qin]
!!      - 01.2018 Brought together from separate files [Yann Pouillon]
! *****************************************************************************
module nao2gto_io

  use nao2gto_common

  implicit none

  private

  public :: &
&   nao2gto_print_info, &
&   nao2gto_read_info, &
&   nao2gto_transfer

contains

! *****************************************************************************
!> \brief Displays the values of the specified Hartree-Fock excahnge data
!!        structure
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] hfx_opts: data structure storing Hartree-Fock exchange
!!                           parameters
! *****************************************************************************
  subroutine nao2gto_print_info(hfx_opts)

    use parallel, only: node
#ifdef BSC_CELLXC
    use bsc_xcmod, only: nXCfunc, XCfunc, XCauth
#else
    use siestaXC, only: getXC
#endif /* BSC_CELLXC */
    use nao2gto_common
    use nao2gto_types, only: hfx_options_type

    implicit none

    ! Arguments
    type(hfx_options_type), intent(in) :: hfx_opts

    ! Local variables
    integer :: nf
#ifndef BSC_CELLXC
    integer :: nXCfunc
    character(len=20) :: XCfunc(10), XCauth(10)
#endif /* BSC_CELLXC */

    ! -------------------------------------------------------------------------

#ifndef BSC_CELLXC
    call getXC(nXCfunc, XCfunc, XCauth)
#endif /* BSC_CELLXC */

    if ( node.eq.0 ) then

      write(*,'(/,a,/)') "nao2gto_print_info: Hartree-Fock information ----------------------------------"

      write(*,'(a,/,a,/,2x,a)') "# *** YAML START ***", "nao2gto_options:", "functionals:"
      do nf=1,nXCfunc
        write(*,'(4x,"- {type: ",a,", name: ",a,", omega: ",E10.3,", potential_type: ",I2,"}")') &
&         trim(XCfunc(nf)), trim(XCauth(nf)), hfx_opts%omega, &
&         hfx_opts%potential_type
      enddo

      write(*,'(2x,a)') "tolerances:"
      write(*,'(4x,"eps_schwarz: ",E10.3)') hfx_opts%eps_schwarz
      write(*,'(4x,"eps_pairlist: ",E10.3)') hfx_opts%eps_pairlist
      write(*,'(4x,"eps_stored: ",E10.3)') hfx_opts%eps_stored
      write(*,'(4x,"eps_farfield: ",E10.3)') hfx_opts%eps_far

      write(*,'(2x,a)') "tweaks:"
      write(*,'(4x,"on_the_fly: ",L1)') hfx_opts%on_the_fly
      write(*,'(4x,"dm_truncate: ",L1)') hfx_opts%DM_trunc
      write(*,'(4x,"far_field: ",L1)') hfx_opts%far_field

      write(*,'(a)') "# *** YAML STOP ***"
      write(*, '(/,a,/)') "nao2gto_print_info: END -------------------------------------------------------"

    endif

  end subroutine  nao2gto_print_info

! *****************************************************************************
!> \brief Reads Hartree-Fock exchange parameters from the FDF input
!!
!! This routine reads the values of the Hartree-Fock exchange parameters
!! from the FDF input of SIESTA and strores them into the specified data
!! structure. It also sets the default values for parameters not present
!! in the input file.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[out] hfx_opts: data structure storing Hartree-Fock exchange
!!                            parameters
! *****************************************************************************
  subroutine nao2gto_read_info(hfx_opts)

    use fdf, only: fdf_get
#ifdef BSC_CELLXC
    use bsc_xcmod, only: nXCfunc, XCfunc, XCauth
#else
    use siestaXC, only: getXC
#endif /* BSC_CELLXC */
    use nao2gto_common
    use nao2gto_types, only: hfx_options_type

    implicit none

    ! Arguments
    type(hfx_options_type), intent(out) :: hfx_opts

    ! Local variables
    integer :: nf
#ifndef BSC_CELLXC
    integer :: nXCfunc
    character(len=20) :: XCfunc(10), XCauth(10)
#endif /* BSC_CELLXC */

    ! -------------------------------------------------------------------------

#ifndef BSC_CELLXC
    call getXC(nXCfunc, XCfunc, XCauth)
#endif /* BSC_CELLXC */

    do nf = 1,nXCfunc
      if ( ((XCauth(nf).eq.'pbe0') .or. (XCauth(nf).eq.'PBE0')) &
&          .and. (XCfunc(nf).eq.'GGA') ) then
        hfx_opts%potential_type = do_hfx_potential_coulomb
      else
        if ( ((XCauth(nf).eq.'hse06') .or. (XCauth(nf).eq.'HSE06')) &
&            .and. (XCfunc(nf).eq.'GGA') ) then
          hfx_opts%potential_type = do_hfx_potential_short
        else
          hfx_opts%potential_type = do_hfx_potential_truncated
        endif
      endif
    enddo
    hfx_opts%omega = fdf_get("Omega", 0.11d0)

    hfx_opts%on_the_fly = fdf_get("On_the_fly", .false.)
    hfx_opts%eps_schwarz = fdf_get('Schwarz.Tolerance', 1.0d-6)
    hfx_opts%eps_stored = fdf_get('ERIs.store', 1.0d-6)
    hfx_opts%eps_pairlist = fdf_get('Pairlist.Tolerance', 1.0d-6)
    hfx_opts%eps_far = fdf_get('Farfield.Tolerance', 1.0d-6)
    hfx_opts%far_field = fdf_get("Farfield", .true.)
    hfx_opts%eri_on_disk = fdf_get("ERIsOnDisk", .false.)
    hfx_opts%DM_trunc = fdf_get("DM_trunc", .true.)
    hfx_opts%dump_fit_data = fdf_get("HFX.DumpFitData", .true.)

  end subroutine nao2gto_read_info

! *****************************************************************************
!> \brief Reads coefficients and builds the corresponding spherical Gaussians
!!
!! This routine reads the coefficients of Natural Atomic Orbitals (NAO) onto
!! Gaussian-Type Orbitals (GTO) from the FDF input of SIESTA and builds the
!! corresponding spherical Gaussians.
!!
!! \author Honghui Shang
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2010 Imported and connected to atm_transfer [Honghui Shang]
!!      - 12.2013 Modified for input from FDF format [Xinming Qin]
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] l_max: ...
!! \param[in] co: ...
!! \param[out] orbtramat: ...
! *****************************************************************************
  subroutine nao2gto_transfer(hfx_opts)

    use units,     only: pi
    use alloc,     only: re_alloc, de_alloc
    use atm_types, only: maxnorbs, nspecies
    use atm_types, only: species, species_info
    use atm_types, only: ncon_max,l_max,nco,nso,ncosum,co,coset,indco, &
                         maxn_orbnl,maxn_contract
    use atmfuncs,  only: lofio,mofio,labelfis
    use atomlist,  only: rmaxo
    use basis_specs, only: label2species
    use chemical,  only: species_label
    use m_io,      only: io_assign
    use radial,    only: rad_get
    use sys,       only: die
    use fdf
    use parsing
    use nao2gto_common
    use nao2gto_types, only: hfx_options_type
    use nao2gto_utils, only: dfac, exp_radius, exp_radius_very_extended
    use nao2gto_transform, only: calc_c2s_matrix, cphi2sphi, orbtramat_type
    use gaufre, only: gaufre_driver_t, gaufre_driver_init, &
&                     gaufre_driver_find, gaufre_driver_free, &
&                     GAUFRE_LM_MINPACK

    implicit none

    ! Arguments
    type(hfx_options_type), intent(in) :: hfx_opts

    ! Local variables
    character(len=8) :: fit_str
    character(len=132) :: msg
    logical :: do_fit
    integer  :: fit_fd, isp, io, iorbrd, inlz, l, m, n, jnlz, &
&     ipgf, irad, itry, jpgf, npts, lx, ly, lz, nco_max, nso_max, &
&     ico, ifdfblk, i_cphi, ren, rel, rezeta, recoeffs
    integer  :: num_cphi(maxnorbs)
    real(dp) :: expzet, fnorm, zeta, zetb, prefac, gcca, gccb, &
&     gauss_resid, dummy, fit_rad, fit_orb, fit_amin, fit_amax, fit_sig3, &
&     new_sig3
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline
    type(species_info), pointer :: spp

    logical, dimension(:), pointer    :: chk_coeffs => null()
    integer, dimension(:,:), pointer :: hfx_contract => null()
    real(dp), dimension(:,:), pointer :: gauss_fit => null()
    real(dp), dimension(:,:,:), pointer :: nao2gto_zeta => null()
    real(dp), dimension(:,:,:), pointer :: nao2gto_coefficient => null()
    real(dp), dimension(:), pointer   :: rad_pts => null()
    real(dp), dimension(:), pointer   :: orb_pts => null()
    type(orbtramat_type), dimension(:), pointer :: orbtramat => null()
    type(gaufre_driver_t) :: fit_driver

    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    !> Step 1: Initialization of nco, nso, co, coset, orbtramat, etc
    ! -------------------------------------------------------------------------
    ncosum(-1:l_max)=0
    do l=0,l_max
      nco(l) = (l + 1)*(l + 2)/2
      nso(l) = 2*l + 1
      ncosum(l) = ncosum(l-1) + nco(l)
    end do

    do lx=0,l_max
      do ly=0,l_max
        do lz=0,l_max
          l = lx + ly + lz
          if ( l > l_max ) cycle
          co(lx,ly,lz) = 1 + (l - lx)*(l - lx + 1)/2 + lz
          coset(lx,ly,lz) = ncosum(l-1) + co(lx,ly,lz)
        end do
      end do
    end do

    indco(:,:) = 0
    do l=0,l_max
      do lx=0,l
        do ly=0,l-lx
          lz = l - lx - ly
          indco(1:3,coset(lx,ly,lz)) = (/lx,ly,lz/)
        end do
      end do
    end do

    ! One cannot use re_alloc with orbtramat, because it is a vector
    ! of structured types
    allocate(orbtramat(0:l_max))
    do l=0,l_max
      nco_max = nco(l)
      nso_max = nso(l)
      write(msg,'("orbtramat(",I1,")%c2s")') l
      nullify(orbtramat(l)%c2s)
      call re_alloc(orbtramat(l)%c2s, 1, nso_max, 1, nco_max, &
&       name=trim(msg), routine='nao2gto_transfer')
    enddo
    call calc_c2s_matrix(l_max, co, orbtramat)

    call re_alloc(hfx_contract, 1, maxn_orbnl, 1, nspecies, &
&     "nao2gto_transfer")
    call re_alloc(nao2gto_zeta, 1, maxn_contract, 1, maxn_orbnl, &
&     1, nspecies, "nao2gto_transfer")
    call re_alloc(nao2gto_coefficient, 1, maxn_contract, 1, maxn_orbnl, &
&     1, nspecies, "nao2gto_transfer")

    hfx_contract(:,:) = 0
    nao2gto_zeta(:,:,:) = 0.0_dp
    nao2gto_coefficient(:,:,:) = 0.0_dp

    ! -------------------------------------------------------------------------
    !> Step 2: Read Hartree-Fock exchange parameters from FDF input file
    !!
    !! \note See ldau_specs.f to understand how to read blocks
    ! -------------------------------------------------------------------------
    do_fit = .true.
    if (fdf_block('NAO2GTO', bfdf)) then

      write(*,'("[DEBUG] NAO2GTO: ",A)') "Will not fit the orbitals"
      do_fit = .false.

      ifdfblk = 1
      do while(fdf_bline(bfdf,pline))

        !> Step 2.a(nofit): Read species
        if ( .not. fdf_bmatch(pline,'ni') ) then
          write(msg,'(A," (line ",I4,")")') &
&           'Wrong format in NAO2GTO', ifdfblk
          call die(trim(msg))
        endif
        isp = label2species(fdf_bnames(pline,1))
        if (isp .eq. 0) then
          write(*,'(a,1x,a)') &
&           'WRONG species symbol in NAO2GTO:', &
&         trim(fdf_bnames(pline,1))
          call die()
        endif
        spp => species(isp)
        ifdfblk = ifdfblk + 1

        !> Step 2.b(nofit): Prepare variables that will tell whether
        !! all coefficients for all orbitals have been read
        call re_alloc(chk_coeffs, 1, spp%n_orbnl, "nao2gto_transfer")
        chk_coeffs(:) = .false.

        !> Step 2.c(nofit): Read data for each (n, l, zeta) triplet
        !! \note Coefficients can be provided in any order.
        do iorbrd=1,spp%n_orbnl

          !> Step 2.c.1(nofit): Read information about which orbital the
          !! coefficients correspond to
          if (.not. fdf_bline(bfdf, pline)) &
&           call die('Not enough information on the Gaussian expansion')
          if (fdf_bmatch(pline,'iiii')) then
            ren = fdf_bintegers(pline,1)
            rel = fdf_bintegers(pline,2)
            rezeta = fdf_bintegers(pline,3)
            recoeffs = fdf_bintegers(pline,4)
            write(*,fmt='(4(A,I2))') "[DEBUG] NAO2GTO: expecting ", recoeffs, &
&             " coefficients for n=", ren, ", l=", rel, ", zeta=", rezeta
          else
            write(msg,'(A," (line ",I4,")")') &
&             'Wrong format in NAO2GTO', ifdfblk
            call die(trim(msg))
          endif

          !> Step 2.c.2(nofit): Translate n, l, zeta into an orbital index
          inlz = orb_nlz_to_index(spp, ren, rel, rezeta)
          if ( inlz == -1 ) then
            write(msg,'(A,3(1X,I2),A," (line ",I4,")")') &
&             'Could not find orbital(', ren, rel, rezeta, ') in SIESTA', &
&             ifdfblk
            call die(trim(msg))
          end if
          hfx_contract(inlz,isp) = recoeffs
          ifdfblk = ifdfblk + 1

          !> Step 2.c.3(nofit): Check for duplicates
          if ( .not. chk_coeffs(inlz) ) then
            chk_coeffs(inlz) = .true.
          else
            write(msg,'(3(A),3(1X,I2),A," (line ",I4,")")') &
&             'Duplicate coefficients for ', spp%symbol, ' orbital(', &
&             ren, rel, rezeta,') in NAO2GTO block', ifdfblk
            call die(trim(msg))
          end if

          !> Step 2.c.3(nofit): Read Gaussian coefficients
          do jnlz=1,hfx_contract(inlz,isp)
            if (.not. fdf_bline(bfdf, pline)) &
&               call die('Not enough information on the Gaussian expansion')
            if (fdf_bmatch(pline,'vv')) then
              nao2gto_zeta(jnlz,inlz,isp) = fdf_bvalues(pline,1)
              nao2gto_coefficient(jnlz,inlz,isp) = fdf_bvalues(pline,2)
            else
              write(msg,'(A," (line ",I4,")")') &
&                 'Wrong format in NAO2GTO', ifdfblk
              call die(trim(msg))
            endif
            ifdfblk = ifdfblk + 1
          enddo

        enddo ! iorbrd=1,spp%n_orbnl

        !> Step 2.d(nofit): Check that all coefficients for the
        !! current species have actually been read
        do iorbrd=1,spp%n_orbnl
          if ( .not. chk_coeffs(iorbrd) ) then
            write(msg,'(A,3(1X,I2),A," (line ",I4,")")') &
&             'Missing coefficients for orbital(', ren, rel, rezeta, &
&             ') in NAO2GTO block', ifdfblk
            call die(trim(msg))
          end if
        end do
        call de_alloc(chk_coeffs, 'chk_coeffs', 'nao2gto_transfer')

      enddo ! while fdf_bline(bfdf,pline)

    else

      write(*,'("[DEBUG] NAO2GTO: ",A)') "Will fit the orbitals"

      do isp=1,nspecies

        spp => species(isp)
        spp%label = species_label(isp)

        do inlz=1,spp%n_orbnl

          !> Step 2.a(fit): Translate rad_func representation of the orbital
          !! \bug Hard-coded the array size!
          call re_alloc(gauss_fit, 1, 2, 1, maxn_contract, &
&                 'gauss_fit', 'nao2gto_transfer')
          call re_alloc(rad_pts, 1, NPTS_GTO, 'rad_pts', 'nao2gto_transfer')
          call re_alloc(orb_pts, 1, NPTS_GTO, 'orb_pts', 'nao2gto_transfer')
          do irad=1,NPTS_GTO
            rad_pts(irad) = spp%orbnl(inlz)%cutoff * real(irad-1, dp) / &
&                   real(NPTS_GTO, dp)
            call rad_get(spp%orbnl(inlz), rad_pts(irad), orb_pts(irad), dummy)
          enddo

          !> Step 2.c(fit): Compute acceptance criteria
          !!
          !! \note
          !!   - The minimum exponent must correspond to at least
          !!     2 sigmas within the cutoof radius.
          !!   - The maximum exponent must satisfy the
          !!     Shannon-Nyquist theorem.
          !!
          fit_amin = 1.0_dp / (spp%orbnl(inlz)%cutoff**2)
          fit_amax = real((NPTS_GTO - 1)**2, dp) / &
&           (18.0_dp*(spp%orbnl(inlz)%cutoff**2))
          write(*,'("[DEBUG]",2(1X,A," = ",E15.5))') &
&           "alpha_min", fit_amin, "alpha_max", fit_amax
          fit_sig3 = 3.0_dp
          jpgf = 0

          !> Step 2.c(fit): Perform Gaussian fitting and filter results
          write(*,'("[DEBUG] CONTRACT:",I4)') &
&                 hfx_contract(inlz,isp)
          write(*,'("[DEBUG] GAUSS_FIT:",2(I4))') &
&                 shape(gauss_fit)
          write(*,'("[DEBUG] RAD:",I4," ORB:",I4)') &
&                 size(rad_pts,1), size(orb_pts,1)
          do ipgf=3,6
            do itry=1,20
              call gaufre_driver_init(fit_driver, GAUFRE_LM_MINPACK, &
&               ipgf, NPTS_GTO, rad_pts, orb_pts)
              call gaufre_driver_find(fit_driver)
              if ( fit_driver%fit_data%converged ) then
                if ( minval(fit_driver%fit_data%coeffs(1,:)) > fit_amin .and. &
&                    maxval(fit_driver%fit_data%coeffs(1,:)) < fit_amax ) then
                  new_sig3 = 3.0_dp - spp%orbnl(inlz)%cutoff * &
&                   sqrt(2.0_dp*minval(fit_driver%fit_data%coeffs(1,:)))
                  if ( abs(new_sig3) < abs(fit_sig3) ) then
                    do jnlz=1,ipgf
                      nao2gto_zeta(jnlz,inlz,isp) = &
&                       fit_driver%fit_data%coeffs(1,jnlz)
                      nao2gto_coefficient(jnlz,inlz,isp) = &
&                       fit_driver%fit_data%coeffs(2,jnlz)
                    end do
                    fit_sig3 = new_sig3
                    jpgf = ipgf
                  end if
                end if
              end if
              call gaufre_driver_free(fit_driver)
            end do
          end do
          hfx_contract(inlz,isp) = jpgf
          write(*,'("[DEBUG] Final distance from 3 sigmas: ",E15.5)') fit_sig3

          !> Step 2.d(fit): Clean-up the mess
          call de_alloc(gauss_fit, 'gauss_fit', 'nao2gto_transfer')
          call de_alloc(rad_pts, 'rad_pts', 'nao2gto_transfer')
          call de_alloc(orb_pts, 'orb_pts', 'nao2gto_transfer')

          !> Step 2.e(fit): Check the validity of the obtained coefficients
          if ( maxval(nao2gto_zeta(:,inlz,isp)) < fit_amin ) then
            write(msg, '(A," (n=",I2,", l=",I2,", zeta=",I2,")")') &
&             "Gaussian fitting of the orbitals failed for", &
&             spp%orbnl_n(inlz), spp%orbnl_l(inlz), spp%orbnl_z(inlz)
            call die(trim(msg))
          end if

        end do   ! inlz

      end do   ! isp

    endif   ! NAO2GTO block

    ! -------------------------------------------------------------------------
    !> Step 3: Process input data for each species
    ! -------------------------------------------------------------------------
    do isp=1,nspecies
      spp => species(isp)
      inlz = 0

      write(*,*) "[DEBUG] ", spp%symbol, " norbs=", spp%norbs
      write(*,*) "[DEBUG]   hfx_contract ", shape(hfx_contract)

      !> Step 3.a: Store GTOs for normal orbitals (m == -l)
      !!
      !! \bug Did not check the use of negative Z numbers
      !!      in original implementation
      do io=1,spp%norbs
        l = spp%orb_l(io)
        m = spp%orb_m(io)
        if ( m .ne. -l ) cycle   ! not a normal orbital
        inlz = inlz + 1
        spp%orbnl_contract(inlz) = hfx_contract(inlz,isp)
        do ipgf=1,spp%orbnl_contract(inlz)
          spp%orbnl_zeta(ipgf,inlz) = nao2gto_zeta(ipgf,inlz,isp)
          spp%orbnl_coefficient(ipgf,inlz) = nao2gto_coefficient(ipgf,inlz,isp)
        enddo
        spp%orbnl_adjoined_zeta(inlz) = &
&         minval(spp%orbnl_zeta(1:spp%orbnl_contract(inlz),inlz))
      enddo

      !> Step 3.b: Add cutoff radius for primitive GTOs and contracted
      !! GTOs, actually PAO. PAO itself has cutoff radius in siesta, we
      !! use the shell cutoff just for comparison
      !! (added by Xinming Qin, Oct. 2018).
      spp%pgf_radius(1:maxn_contract,1:maxn_orbnl)=0.0d0
      spp%shell_radius(1:maxn_orbnl)=0.0d0
      inlz = 0
      do io=1,spp%norbs
        l = spp%orb_l(io)
        m = spp%orb_m(io)
        if ( m /= -l ) cycle   ! Not a normal orbital
        inlz = inlz+1
        do ipgf=1,spp%orbnl_contract(inlz)
          gcca = nao2gto_coefficient(ipgf,inlz,isp)
          zeta = nao2gto_zeta(ipgf,inlz,isp)
          spp%pgf_radius(ipgf,inlz)= exp_radius(l, zeta, 1.0E-5_dp, gcca)
        enddo
        spp%shell_radius(inlz) = &
&         maxval(spp%pgf_radius(1:spp%orbnl_contract(inlz),inlz))
      enddo
      spp%kind_radius = maxval(spp%shell_radius(1:maxn_orbnl))

      !> Step 3.c: Compute the number of cphi coefficients for each orbital
      !! and store their total number in i_cphi
      num_cphi(:) = 0
      i_cphi = 1
      do io=1,spp%norbs
        if ( spp%orb_m(io) .lt. 2 ) then
          num_cphi(i_cphi) = io
          i_cphi = i_cphi + 1
        else if ( spp%orb_m(io) .eq. 2 ) then
          num_cphi(i_cphi) = io
          num_cphi(i_cphi+1) = io
          i_cphi = i_cphi + 2
        else if ( spp%orb_m(io) .eq. 3 ) then
          num_cphi(i_cphi) = io
          num_cphi(i_cphi+1) = io
          num_cphi(i_cphi+2) = io
          num_cphi(i_cphi+3) = io
          i_cphi = i_cphi + 4
        endif
      enddo
      i_cphi = i_cphi - 1

      !> Step 3.d: Initialize indices for cartesian coefficients
      !!
      !! \note Here, we only consider the lmax <= 3 case (s, p, d, and f)
      if ( spp%lmax_basis .lt. 2 ) then
        spp%norbs_cphi = spp%norbs   !< Coefficients for s and p orbitals
      else
        spp%norbs_cphi = i_cphi !< Coefficients for d and f orbitals
      endif

      write(*,'("num_cphi ",6(1X,I4))') num_cphi(:)
      do i_cphi=1,spp%norbs_cphi
        io = num_cphi(i_cphi)
        spp%orb_n_cphi(i_cphi) = spp%orb_n(io)
        spp%orb_l_cphi(i_cphi) = spp%orb_l(io)
        spp%orb_index_cphi(i_cphi) = spp%orb_index(io)
        write(*,'(A,A,5(1X,A,"=",I4))') "[DEBUG][orb_index] ", spp%symbol, &
&         "norbs_cphi", spp%norbs_cphi, &
&         "io", io, &
&         "n", spp%orb_n(io), &
&         "l", spp%orb_l(io), &
&         "orb_index", spp%orb_index(io)
      enddo

      io = 0
      do inlz=1,spp%n_orbnl
        l = spp%orbnl_l(inlz)

        if ( inlz .eq. 1 ) then
          spp%orbnl_index_cphi(inlz) = 1
          spp%orbnl_index_sphi(inlz) = 1
        else
          spp%orbnl_index_cphi(inlz) = spp%orbnl_index_cphi(inlz-1) &
&                                    + nco(spp%orbnl_l(inlz-1))
          spp%orbnl_index_sphi(inlz) = spp%orbnl_index_sphi(inlz-1) &
&                                    + nso(spp%orbnl_l(inlz-1))
        endif

        do ico=ncosum(l-1)+1,ncosum(l)
          io = io + 1
          spp%orb_cphi_lx(io) = indco(1,ico)
          spp%orb_cphi_ly(io) = indco(2,ico)
          spp%orb_cphi_lz(io) = indco(3,ico)
        enddo
      enddo

      !> Step 3.e: Compute norms of cartesian GTOs
      do io=1,spp%norbs_cphi
        l = spp%orb_l_cphi(io)
        expzet = 0.5_dp*REAL(2*l + 3,dp)
        fnorm = 0.0_dp
        inlz = spp%orb_index_cphi(io)
        do ipgf=1,spp%orbnl_contract(inlz)
          gcca = spp%orbnl_coefficient(ipgf,inlz)
          zeta = spp%orbnl_zeta(ipgf,inlz)
          do jpgf=1,spp%orbnl_contract(inlz)
            gccb = spp%orbnl_coefficient(jpgf,inlz)
            zetb = spp%orbnl_zeta(jpgf,inlz)
            fnorm = fnorm + gcca*gccb/((zeta + zetb)**expzet)
          end do
        end do

        fnorm = (0.5_dp**l)*(pi**1.5_dp)*fnorm
        lx = spp%orb_cphi_lx(io)
        ly = spp%orb_cphi_ly(io)
        lz = spp%orb_cphi_lz(io)
        prefac = dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1)
        spp%norm_cphi(io) = 1.0_dp/SQRT(prefac*fnorm)
      enddo

      !> Step 3.f: Compute cphi and sphi
      spp%cphi(1:ncon_max,1:spp%norbs_cphi) = 0.0_dp
      spp%sphi(1:ncon_max,1:spp%norbs) = 0.0_dp

      do io=1,spp%norbs_cphi
        inlz = spp%orb_index_cphi(io)
        lx = spp%orb_cphi_lx(io)
        ly = spp%orb_cphi_ly(io)
        lz = spp%orb_cphi_lz(io)
        do jnlz=1,spp%orbnl_contract(inlz)
          spp%cphi((jnlz-1)*nco(spp%orb_l_cphi(io))+co(lx,ly,lz),io) =  &
&           spp%orbnl_coefficient(jnlz,inlz) * spp%norm_cphi(io)
        enddo
      enddo

      call cphi2sphi(ncon_max, spp%norbs_cphi, spp%norbs, l_max, &
&                    spp%n_orbnl, spp%orbnl_l, nco, nso, &
&                     spp%orbnl_index_cphi, spp%orbnl_index_sphi, &
&                     spp%cphi, spp%sphi, orbtramat)

      inlz = 0
      do io=1,spp%norbs
        l = spp%orb_l(io)
        m = spp%orb_m(io)
        if ( m .ne. -l ) cycle     ! not a normal orbital
        inlz = inlz+1
        spp%orbnl_contraction_coeff(inlz) =  &
&         maxval([(sum(abs(spp%sphi(1:nco(l)*spp%orbnl_contract(inlz),jnlz))), &
&           jnlz=io,io+nso(l)-1)])
      enddo

    enddo   ! is=1,nspecies

    ! -------------------------------------------------------------------------
    !> Step 4: Report about fitting parameters
    ! -------------------------------------------------------------------------
    write(*,'(/,a,/)') "nao2gto_transfer: NAO2GTO fitting information ---------------------------------"
    do isp=1,nspecies
      spp => species(isp)
      spp%label = species_label(isp)

      write(*,*) 'Write NAO2GTO fitting information:'

      write(*,'(2x,a)') &
&       "#(species label, n, l, z, ngto, is_polarized, popul, index)"
      write(*,'(2x,2(a16))') 'zeta', 'coefficent'

      do inlz=1,spp%n_orbnl
        write(*,'(2a,4i3,l5,e15.5,i3)') "# ", trim(spp%label), &
&         spp%orbnl_n(inlz), spp%orbnl_l(inlz), spp%orbnl_z(inlz),&
&         spp%orbnl_contract(inlz), spp%orbnl_ispol(inlz), &
&         spp%orbnl_pop(inlz), inlz

        do n=1,spp%orbnl_contract(inlz)
          write(*,'(2x,2f16.9)') spp%orbnl_zeta(n,inlz),&
&           spp%orbnl_coefficient(n,inlz)
        enddo
      enddo
    enddo

    write(*,'(/,a,/)') "nao2gto_transfer: -------------------------------------------------------------"

    if ( hfx_opts%dump_fit_data ) then

      call io_assign(fit_fd)
      open(unit=fit_fd, file="nao2gto_fit.yml", form="formatted", &
&       position="rewind", action="write", status="unknown")
      write(unit=fit_fd, fmt='(A/"---")') "%YAML 1.2"

      write(unit=fit_fd, fmt='(/A,":")') "hfx_options"
      write(unit=fit_fd, fmt='(2X,A,": ",I2)') "potential_type", &
&       hfx_opts%potential_type
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "omega", &
&       hfx_opts%omega
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_far", &
&       hfx_opts%eps_far
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_pairlist", &
&       hfx_opts%eps_pairlist
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_schwarz", &
&       hfx_opts%eps_schwarz
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_stored", &
&       hfx_opts%eps_stored
      if ( hfx_opts%DM_trunc ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "DM_trunc", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "DM_trunc", "False"
      end if
      write(unit=fit_fd, fmt='(2X,A,": ",A)') "dump_fit_data", "True"
      if ( hfx_opts%eri_on_disk ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "eri_on_disk", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "eri_on_disk", "False"
      end if
      if ( hfx_opts%far_field ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "far_field", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "far_field", "False"
      end if
      if ( hfx_opts%on_the_fly ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "on_the_fly", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "on_the_fly", "False"
      end if

      write(unit=fit_fd, fmt='(/A,":")') "orbitals"
      do isp=1,nspecies
        spp => species(isp)
        spp%label = species_label(isp)

        do inlz=1,spp%n_orbnl
          write(unit=fit_fd, fmt='(2X,"- ",A,": ",A)') "species", &
&           trim(adjustl(spp%symbol))
          write(unit=fit_fd, fmt='(4X,A,": ",I4)') "qn_n", spp%orbnl_n(inlz)
          write(unit=fit_fd, fmt='(4X,A,": ",I4)') "qn_l", spp%orbnl_l(inlz)
          write(unit=fit_fd, fmt='(4X,A,": ",I4)') "zeta", spp%orbnl_z(inlz)
          if ( spp%orbnl_ispol(inlz) ) then
            write(unit=fit_fd, fmt='(4X,A,": ",A)') "polarized", "True"
          else
            write(unit=fit_fd, fmt='(4X,A,": ",A)') "polarized", "False"
          end if
          write(unit=fit_fd, fmt='(4X,A,": ",E15.5)') &
&           "population", spp%orbnl_pop(inlz)

          if ( do_fit ) then
            write(unit=fit_fd, fmt='(4X,A,": ",A)') "do_fit", "True"
          else
            write(unit=fit_fd, fmt='(4X,A,": ",A)') "do_fit", "False"
          end if

          write(unit=fit_fd, fmt='(4X,A,":")') "points"
          do irad=1,NPTS_GTO
            fit_rad = spp%orbnl(inlz)%cutoff * real(irad-1, dp) / &
&               real(NPTS_GTO, dp)
            call rad_get(spp%orbnl(inlz), fit_rad, fit_orb, dummy)
            write(unit=fit_fd, fmt='(6X,"- [",E20.8,", ",E20.8,"]")') &
&               fit_rad, fit_orb
          enddo

          write(unit=fit_fd, fmt='(4X,A,":")') "fit"
          do jnlz=1,hfx_contract(inlz,isp)
            write(unit=fit_fd, fmt='(6X,"- [",E20.8,", ",E20.8,"]")') &
&             nao2gto_zeta(jnlz,inlz,isp), nao2gto_coefficient(jnlz,inlz,isp)
          end do
        end do
      end do

      write(unit=fit_fd, fmt='(/A)') "..."
      close(unit=fit_fd)

    end if   ! hfx_opts%dump_fit_data

    ! -------------------------------------------------------------------------
    !> Step 5: Properly free memory
    ! -------------------------------------------------------------------------
    do l=0,l_max
      write(msg,'("orbtramat(",I1,")%c2s")') l
      call de_alloc(orbtramat(l)%c2s, name=trim(msg), &
&       routine='nao2gto_transfer')
      nullify(orbtramat(l)%c2s)
    enddo
    deallocate(orbtramat)
    nullify(orbtramat)

    call de_alloc(hfx_contract, name="hfx_contract", &
&     routine="nao2gto_transfer")
    call de_alloc(nao2gto_zeta, name="nao2gto_zeta", &
&     routine="nao2gto_transfer")
    call de_alloc(nao2gto_coefficient, name="nao2gto_coefficient", &
&     routine="nao2gto_transfer")

  end subroutine nao2gto_transfer

  function orb_nlz_to_index(spp, orb_n, orb_l, orb_z) result(orb_index)

    use atm_types, only: species_info

    implicit none

    ! Arguments
    type(species_info), intent(in) :: spp
    integer, intent(in) :: orb_n
    integer, intent(in) :: orb_l
    integer, intent(in) :: orb_z

    ! Local variables
    integer :: orb_index
    integer :: iorb

    orb_index = -1

    do iorb=1,spp%n_orbnl
      if ( spp%orbnl_n(iorb) == orb_n ) then
        if ( spp%orbnl_l(iorb) == orb_l ) then
          if ( spp%orbnl_z(iorb) == orb_z ) then
            orb_index = iorb
            exit
          end if
        end if
      end if
    end do

  end function orb_nlz_to_index

end module nao2gto_io
