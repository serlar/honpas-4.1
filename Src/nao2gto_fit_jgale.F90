#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! ****************************************************************************
!> \brief This module contains the subroutines required to fit Gaussians to
!!        a SIESTA orbital.
!!
!! \author
!!     Julian Gale, NRI, Curtin University of Technology, April 2005
! ****************************************************************************
module nao2gto_fit_jgale

  use precision, only: dp

  implicit none

  private

  integer,               save :: ll
  real(dp),              save :: rcut
  real(dp), allocatable, save :: orbital(:)
  real(dp), allocatable, save :: radius(:)

  public :: fitgauss

contains

  ! **************************************************************************
  !> \brief Main fitting driver that calls Minpack subroutine
  !!
  !! \author
  !!     Julian Gale, NRI, Curtin University of Technology, [2005.04.30]
  !!
  !! \par History
  !!     2005.04.30: Creation [Julian Gale]
  !!     2018.07.18: Upgrade to Fortran 2003 and Doxygen [Yann Pouillon]
  !!
  !! \param[in] nradialpoints: number of points in the radial mesh
  !! \param[in] radialpoints: radial mesh points
  !! \param[in] orbitalonradial: values of the orbital on the radial mesh points
  !! \param[in] ngaussians: number of Gaussians in the expansion
  !! \param[in] l: angular momentum of the orbital
  !! \param[out] gaussianprefactor: coefficients of the individual Gaussians
  !! \param[out] gaussianzeta: exponents of the individual Gaussians
  !! \param[out] sumofsquares: residual norm of the Gaussian expansion
  !! \param[out] info: return value of MINPACK's LMDIF
  ! **************************************************************************
  subroutine fitgauss(nradialpoints, radialpoints, orbitalonradial, &
&     ngaussians, l, gaussianprefactor, gaussianzeta, sumofsquares, info)

    implicit none

    ! Arguments
    integer,    intent(in)    :: nradialpoints
    integer,    intent(in)    :: ngaussians
    integer,    intent(in)    :: l
    real(dp),   intent(inout) :: gaussianprefactor(ngaussians)
    real(dp),   intent(inout) :: gaussianzeta(ngaussians)
    real(dp),   intent(in)    :: orbitalonradial(nradialpoints)
    real(dp),   intent(in)    :: radialpoints(nradialpoints)
    integer,    intent(out)   :: info
    real(dp),   intent(out)   :: sumofsquares

    ! Local variables
    integer                   :: i
    integer                   :: j
    integer                   :: ldfjac
    integer                   :: nfev
    integer                   :: ngvar
    integer                   :: maxfev
    integer                   :: mode
    integer                   :: n
    integer                   :: nprint
    integer,allocatable, save :: ipvt(:)
    real(dp), allocatable, save :: diag(:)
    real(dp), allocatable, save :: f(:)
    real(dp), allocatable, save :: fjac(:,:)
    real(dp), allocatable, save :: qtf(:)
    real(dp), allocatable, save :: x(:)
    real(dp), allocatable, save :: wa1(:)
    real(dp), allocatable, save :: wa2(:)
    real(dp), allocatable, save :: wa3(:)
    real(dp), allocatable, save :: wa4(:)
    real(dp)                    :: alpha
    real(dp)                    :: epsfcn
    real(dp)                    :: factor
    real(dp)                    :: fshift
    real(dp)                    :: ftol
    real(dp)                    :: gtol
    real(dp)                    :: xtol
    real(dp)                    :: zeta
    real(dp)                    :: sumfshift

    ! ------------------------------------------------------------------------

    ! Set size parameters
    ngvar = 2*ngaussians
    ldfjac = nradialpoints

    ! Allocate dependent arrays for minpack
    allocate(ipvt(ngvar))
    allocate(diag(ngvar))
    allocate(f(nradialpoints))
    allocate(fjac(ldfjac,ngvar))
    allocate(qtf(ngvar))
    allocate(x(ngvar))
    allocate(wa1(ngvar))
    allocate(wa2(ngvar))
    allocate(wa3(ngvar))
    allocate(wa4(nradialpoints))

    ! Allocate arrays for module
    allocate(orbital(nradialpoints))
    allocate(radius(nradialpoints))

    ! Set details of fitting data in local arrays
    orbital(1:nradialpoints) = orbitalonradial(1:nradialpoints)
    radius(1:nradialpoints) = radialpoints(1:nradialpoints)
    rcut = radius(nradialpoints)
    ll = l

    ! Set details of fitting functions in local arrays
    !
    !   ngvar = number of Gaussian variables
    !   x     = initial values of variables
    !
    n = 0
    do i = 1,ngaussians
      x(n+1) = gaussianprefactor(i)
      x(n+2) = gaussianzeta(i)
      n = n + 2
    enddo

    ! Initialise values for minpack parameters
    !
    !   ftol / gtol / xtol = convergence tolerances
    !   maxfev             = maximum number of function evaluations
    !   epsfcn             = forward difference step size
    !
    maxfev = 10000
    mode = 1
    nprint = - 1
    epsfcn = 1.0d-6
    factor = 100.0_dp
    ftol = 1.0d-3
    gtol = 1.0d-4
    xtol = 1.0d-4

    ! Call function to perform least squares from Minpack
    call lmdif(nradialpoints, ngvar, x, f, ftol, xtol, gtol, maxfev, &
&       epsfcn, diag, mode, factor, nprint, info, nfev, fjac, ldfjac, &
&       ipvt, qtf, wa1, wa2, wa3, wa4)

    ! Return values to passed arrays
    n = 0
    do i = 1,ngaussians
      gaussianprefactor(i) = x(n+1)
      gaussianzeta(i) = abs(x(n+2))
      n = n + 2
    enddo

    ! Compute final sum of squares
    sumofsquares = 0.0_dp
    do i = 1,nradialpoints
      sumofsquares = sumofsquares + f(i)*f(i)
    enddo
    sumofsquares = sqrt(sumofsquares)

    ! Deallocate memory
    deallocate(radius)
    deallocate(orbital)
    deallocate(wa4)
    deallocate(wa3)
    deallocate(wa2)
    deallocate(wa1)
    deallocate(x)
    deallocate(qtf)
    deallocate(fjac)
    deallocate(f)
    deallocate(diag)
    deallocate(ipvt)

  end subroutine fitgauss

  ! **************************************************************************
  !> \brief Routine called internally by LMDIF1
  !!
  !! Details at: https://www.math.utah.edu/software/minpack/minpack/lmdif1.html
  !!
  !! \par History
  !!      - 11.2017 Created [Yann Pouillon]
  !!      - 01.2018 Added Doxygen documentation [Yann Pouillon]
  !!
  !! \note The interface of this routine is imposed by MINPACK.
  !!
  !! \param[in] M: number of data points (npts in the rest of this module)
  !! \param[in] N: number of variables (2*ngfs in the rest of this module)
  !! \param[in] X: trial set of variables
  !! \param[in,out] FVEC: computed residuals
  !! \param[in,out] IFLAG: controls what to compute (1=residuals,
  !!                       2=Jacobian residuals)
  ! **************************************************************************
  subroutine fcn(nrpts, ngvar, x, f, iflag)

    implicit none

    ! Arguments
    integer,  intent(in)  :: ngvar
    integer,  intent(in)  :: nrpts
    integer,  intent(in)  :: iflag
    real(dp), intent(in)  :: x(ngvar)
    real(dp), intent(out) :: f(nrpts)

    ! Local variables
    integer  :: i
    integer  :: j
    integer  :: n
    integer  :: ngauss
    real(dp) :: alpha
    real(dp) :: zeta
    real(dp) :: fshift
    real(dp) :: fsum

    ! ------------------------------------------------------------------------

    !  Calculate function at points
    f(1:nrpts) = 0.0_dp
    ngauss = ngvar / 2
    n = 0
    do i = 1,ngauss
      alpha = x(n+1)
      zeta = x(n+2)
      fshift = alpha*(rcut**ll)*exp(-abs(zeta)*rcut*rcut)
      fshift=0
      n = n + 2
      do j = 1,nrpts
        f(j) = f(j) + alpha*(radius(j)**ll) * &
&              exp(-abs(zeta)*radius(j)*radius(j)) - fshift
      enddo
    enddo

    !  Compute error
    do j = 1,nrpts
      f(j) = orbital(j) - f(j)
    enddo

  end subroutine fcn

end module nao2gto_fit_jgale
