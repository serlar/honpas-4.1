! *** Module: nao2gto_dm ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Density-matrix utilities for NAO2GTO
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_dm

  use precision,     only: dp
  use atm_types,     only: nso, nco
  use atmfuncs,      only: lofio, mofio
  use atomlist,      only: indxuo
  use listsc_module, only: listsc

  implicit none

  private

  public :: sparse_dense, get_pmax_shell

contains

! *****************************************************************************
!> \brief Builds the dense density matrix from the sparse one
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] nuo: local number of orbitals within the unit cell in this node
!!                 (nuo <= nuotot)
!! \param[in] norb: local number of orbitals within the supercell
!!              (norb=nuotot*ncells)
!! \param[in] maxnd: maximum size of the sparse matrix
!! \param[in] numd(nuo): number of nonzero elements of each row
!! \param[in] listd(maxnd): nonzero hamiltonian-matrix element
!! \param[in] listdptr(nuo): control vector of density matrix
!! \param[in] M_sparse(mxnd): sparse matrix
!! \param[out] M_dense(mxnd): dense matrix
! *****************************************************************************
  subroutine sparse_dense(nspin, nuo, nuotot, norb, maxnd, numd, listdptr, &
&                listd, M_sparse, M_dense)

    use precision, only: dp

    implicit none

    integer , intent(in)  :: nspin, nuo, nuotot, norb, maxnd
    integer , intent(in)  :: numd(nuo), listdptr(nuo), listd(*)
    real(dp), intent(in)  :: M_sparse(maxnd, nspin)
    real(dp), intent(out) :: M_dense(norb,norb,nspin)

    ! Local variables
    integer :: ispin, io, jo, j, iu, ju, ind
    real(dp), dimension(:,:), pointer :: M_aux => null()

    ! -------------------------------------------------------------------------

    allocate(M_aux(nuotot,norb))

    M_aux(:,:) = 0.0_dp
    M_dense(:,:,:) = 0.0_dp

    do ispin=1,nspin
      do io=1,nuotot
        do j=1,numd(io)
          ind = listdptr(io) + j
          jo = listd(ind)
          M_aux(io,jo) = M_sparse(ind,ispin)
          M_dense(io,jo,ispin) = M_aux(io,jo)
        enddo
      enddo

      do io=nuotot+1,norb
        iu = indxuo(io)
        do j=1,numd(iu)
          ind=listdptr(iu) +j
          ju=listd(ind)
          jo=listsc(io,iu,ju)
          M_dense(io,jo,ispin)=M_aux(iu,ju)
        enddo
      enddo
    enddo

    deallocate(M_aux)

  end subroutine sparse_dense

! *****************************************************************************
!> \brief ...
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] nspin: number of spin states
!! \param[in] norb: local number of orbitals within the supercell
!! \param[in] iaorb: ...
!! \param[in] iphorb: ...
!! \param[in] nuo: local number of orbitals within the unit cell in this node
!! \param[in] na: ...
!! \param[in] isa: ...
!! \param[in] maxnd: maximum size of the sparse matrix
!! \param[in] numd: number of nonzero elements of each row
!! \param[in] listdptr: control vector of density matrix
!! \param[in] listd: nonzero hamiltonian-matrix element
!! \param[in] M_dense: dense matrix
!! \param[out] P_max: ...
! *****************************************************************************
  subroutine get_pmax_shell(nspin, norb, iaorb, iphorb, nuo, na, isa, &
&                maxnd, numd, listdptr, listd, M_dense, P_max)

    ! Arguments
    integer , intent(in)  :: nspin, norb, nuo, iaorb(norb), &
&     iphorb(norb), na, isa(na), maxnd, numd(nuo), listdptr(nuo), &
&     listd(*)
    real(dp), intent(in)  :: M_dense(norb, norb, nspin)
    real(dp), intent(out) :: P_max(norb, norb)

    ! Local variables
    integer :: ia, ja, io, ioa, is, jo, joa, js, j, iu, ju, ind, &
&     il, im, jl, jm

    ! -------------------------------------------------------------------------

    do io  = 1, norb
      ia = iaorb(io)
      ioa = iphorb(io)
      is = isa(ia)
      il = lofio(is,ioa)
      im = mofio(is,ioa)
      if ( im .ne. -il ) cycle

      iu = indxuo(io)
      do j=1, numd(iu)
        ind = listdptr(iu) + j
        ju = listd(ind)
        jo = listsc(io,iu,ju)
        ja = iaorb(jo)
        joa = iphorb(jo)
        js = isa(ja)
        jl = lofio(js,joa)
        jm = mofio(js,joa)
        if ( jm .ne. -jl ) cycle
        P_max(io,jo) = P_max(io,jo) &
&         + maxval(abs(M_dense(io:io-1+nso(il),jo:jo-1+nso(jl),1:nspin)))
      enddo
    enddo

  end subroutine get_pmax_shell

end module nao2gto_dm
