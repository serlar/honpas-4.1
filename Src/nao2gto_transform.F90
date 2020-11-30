! *** Module: nao2gto_transform ***

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
module nao2gto_transform

  use precision, only: dp

  implicit none

  private

  public :: &
&   calc_c2s_matrix, &
&   cphi2sphi, &
&   orbtramat_type

    !> \brief Spherical harmonics and orbital transformation matrices
    !!
    !! This type stores spherical harmonics and the corresponding orbital
    !! transformation matrices.
    !!
    !! \par Literature
    !!      H. B. Schlegel, M. J. Frisch, Int. J. Quantum Chem. 54, 83 (1995)
    type :: orbtramat_type
      real(dp), dimension(:,:), pointer :: c2s
    end type orbtramat_type

contains

! *****************************************************************************
!> \brief Computes the cartesian-to-spherical transformation matrix
!!
!! This routine builds the orbital transformation matrix for the
!! conversion from cartesian to spherical orbitals (c2s, formula 15)
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] l_max: ...
!! \param[in] co: ...
!! \param[out] orbtramat: ...
! *****************************************************************************
  subroutine calc_c2s_matrix(l_max, co, orbtramat)

    use precision,     only: dp
    use nao2gto_utils, only: fac

    implicit none

    ! Arguments
    integer, intent(in) :: l_max, co(0:l_max,0:l_max,0:l_max)
    type(orbtramat_type), dimension(0:l_max), intent(out) :: orbtramat

    ! Local variables
    integer :: expo, i, ic, isp, j, k, l, lx, ly, lz, m, ma
    real(dp) :: s, s1, s2

    ! -------------------------------------------------------------------------

    do l=0,l_max
      do lx=0,l
        do ly=0,l-lx
          lz = l - lx - ly
          ic = co(lx,ly,lz)
          do m=-l,l
            isp = l + m + 1
            ma = abs(m)
            j = lx + ly - ma
            if ( (j >= 0) .and. (modulo(j,2) == 0) ) then
              j = j/2
              s1 = 0.0_dp
              do i=0,(l-ma)/2
                s2 = 0.0_dp
                do k=0,j
                  if ( ((m < 0) .and. (modulo(abs(ma-lx),2) == 1)) .or. &
&                      ((m > 0) .and. (modulo(abs(ma-lx),2) == 0)) ) then
                    expo = (ma - lx + 2*k)/2
                    s = (-1.0_dp)**expo*sqrt(2.0_dp)
                  else if ( (m == 0) .and. (modulo(lx,2) == 0) ) then
                    expo = k - lx/2
                    s = (-1.0_dp)**expo
                  else
                    s = 0.0_dp
                  end if
                  s2 = s2 + binomial(j,k)*binomial(ma,lx-2*k)*s
                end do
                s1 = s1 + binomial(l,i)*binomial(i,j)* &
&                    (-1.0_dp)**i*fac(2*l-2*i)/fac(l-ma-2*i)*s2
              end do
              orbtramat(l)%c2s(isp,ic) = &
&               sqrt((fac(2*lx)*fac(2*ly)*fac(2*lz)*fac(l)*fac(l-ma))/ &
&                    (fac(lx)*fac(ly)*fac(lz)*fac(2*l)*fac(l+ma)))*s1/ &
&                    ((2.0_dp**l)*fac(l))
            else
              orbtramat(l)%c2s(isp,ic) = 0.0_dp
            end if
          end do
        end do
      end do
    end do

  end subroutine calc_c2s_matrix

! *****************************************************************************
!> \brief Computes cartesian to spherical coordinates transforms
!!
!! This routine converts Gaussian coefficients expressed in a cartesian
!! coordinates system to spherical coordinates.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] ncon_max: ...
!! \param[in] norb_cphi: ...
!! \param[in] norb_sphi: ...
!! \param[in] l_max: ...
!! \param[in] norb_nl: ...
!! \param[in] orbnl_l: ...
!! \param[in] nco: ...
!! \param[in] nso: ...
!! \param[in] index_cphi: ...
!! \param[in] index_sphi: ...
!! \param[in] cphi: ...
!! \param[in,out] sphi: ...
!! \param[in] orbtramat: ...
!!
!! \note The sphi variable must be declared inout because of LAPACK's DGEMM.
! *****************************************************************************
  subroutine cphi2sphi(ncon_max, norb_cphi, norb_sphi, l_max, norb_nl, &
&                orbnl_l, nco, nso, index_cphi, index_sphi, cphi, sphi, &
&                orbtramat)

    use precision, only: dp

    implicit none

    ! Arguments
    integer, intent(in)     :: ncon_max, norb_cphi, norb_sphi, l_max, norb_nl
    integer, intent(in)     :: orbnl_l(norb_nl), nco(0:l_max), nso(0:l_max), &
&                            index_cphi(norb_nl), index_sphi(norb_nl)
    real(dp), intent(in)    :: cphi(ncon_max,norb_cphi)
    real(dp), intent(inout) :: sphi(ncon_max,norb_sphi)
    type(orbtramat_type), dimension(0:l_max), intent(in) :: orbtramat

    ! Local variables
    integer :: i, l, ncgf, nsgf, first_cgf, first_sgf

    ! -------------------------------------------------------------------------

    do i=1,norb_nl
      l = orbnl_l(i)
      ncgf = nco(l)
      nsgf = nso(l)
      first_cgf = index_cphi(i)
      first_sgf = index_sphi(i)
      call dgemm("N", "T", ncon_max, nsgf, ncgf, 1.0_dp, &
&       cphi(1,first_cgf), ncon_max, orbtramat(l)%c2s(1,1), nsgf, 0.0_dp, &
&       sphi(1,first_sgf), ncon_max)
    enddo

  end subroutine cphi2sphi

! *****************************************************************************
! *** Private routines                                                      ***
! *****************************************************************************

! *****************************************************************************
!> \brief Returns the number of binomial combinations of two integers
!!
!! This function computes the binomial combinations of two integers k and n
!! as (n!/(k!(n-k)!)).
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] n: sample size
!! \param[in] k: number of sample elements selected
!! \return number of binomial combinations of k elements among n
! *****************************************************************************
  function binomial(n,k) result(n_over_k)

    use precision,     only: dp
    use nao2gto_utils, only: fac

    implicit none

    ! Arguments
    integer, intent(in) :: n, k

    ! Local variables
    real(dp) :: n_over_k

    ! -------------------------------------------------------------------------

    if ( (k >= 0) .and. (k <= n) ) then
      n_over_k = fac(n)/(fac(n-k)*fac(k))
    else
      n_over_k = 0.0_dp
    endif

  end function binomial

end module nao2gto_transform
