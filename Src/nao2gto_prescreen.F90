! *** Module: nao2gto_prescreen ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Prescreen ERIs using the Schwarz inequality to get list_ij,
!!        list_uv and list_mn.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_prescreen

  use nao2gto_common

  implicit none

  private

  public :: &
&   build_pair_list, &
&   calc_prescreen_eri, &
&   init_prescreen_eri, &
&   calc_pair_dist_radii, &
&   calc_screening_functions

contains

! *****************************************************************************
!> \brief ...
!!
!! This routine uses the Schwarz inequality to compute the requested values.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] libint_data: Libint data structure
!! \param[in] hfx_opts: data structure storing Hartree-Fock exchange options
!! \param[in] hfx_sys: data structure storing system information
!! \param[out] max_eri: maximum value of eri_prescreen
! *****************************************************************************
  subroutine calc_prescreen_eri(libint_data, hfx_opts, hfx_sys, &
&   maxnh, numh, listhptr, listh, max_eri)

    use atm_types,      only: species, species_info, &
                              l_max, nco, nso
    use parallel,       only: Node, Nodes
    use atmfuncs,       only: lofio, mofio, rcut
    use atomlist,       only: indxuo, rmaxo, lasto
    use atomlist,       only: rco, rmaxkb, lastkb
    use listsc_module,  only: listsc
    use alloc
    use sorting
    use neighbour,      only: jna=>jan, xij, r2ij, maxna=>maxnna
    use neighbour,      only: mneighb
    use nao2gto_contract
    use nao2gto_data
    use nao2gto_index,  only: indexsc
    use nao2gto_libint, only: Libint_t
    use nao2gto_pbc,    only: trans_pbc
    use nao2gto_powell
    use nao2gto_types

    ! Arguments
    type(Libint_t)        , intent(inout) :: libint_data
    type(hfx_options_type), intent(in)    :: hfx_opts
    type(hfx_system_type) , intent(in)    :: hfx_sys
    integer               , intent(in)    :: maxnh
    integer               , intent(in)    :: numh(hfx_sys%nuotot)
    integer               , intent(in)    :: listhptr(hfx_sys%nuotot)
    integer               , intent(in)    :: listh(maxnh)
    real(dp)              , intent(out)   :: max_eri

    ! Local variables
    integer :: ia, ioa, io, iu, is, j, ja, joa, jo, ju, js, ind, &
&     l_i, l_j, m_i, m_j, ncoi, ncoj, npgfi, npgfj
    real(dp) :: ri(3), rj(3), r_temp(3), r_pbc(3)
    type(species_info), pointer :: ispp => null(), jspp => null()
    real(dp), dimension(:,:,:,:), pointer :: eri => null()

    ! -------------------------------------------------------------------------

    do io=1,hfx_sys%nuotot
      ia   =  hfx_sys%iaorb(io)
      is   =  hfx_sys%isa(ia)
      ispp => species(is)
      ioa  =  hfx_sys%iphorb(io)
      l_i  =  lofio(is,ioa)
      m_i  =  mofio(is,ioa)

      if ( m_i .ne. -l_i ) cycle

      ri(:) = hfx_sys%xa(:,ia)
      npgfi = ispp%orbnl_contract(ispp%orb_index(ioa))
      ncoi  = nco(l_i)*npgfi

      do j=1,numh(io)
        ind = listhptr(io) + j
        jo  = listh(ind)
        ja  = hfx_sys%iaorb(jo)
        js  = hfx_sys%isa(ja)
        jspp => species(js)
        joa = hfx_sys%iphorb(jo)
        l_j = lofio(js,joa)
        m_j = mofio(js,joa)

        if ( m_j .ne. -l_j ) cycle

        rj(:)  = hfx_sys%xa(:,ja)
        r_temp = rj - ri
        call trans_pbc(r_temp, hfx_sys%cell, hfx_sys%cell_r, r_pbc)
        rj = ri + r_pbc

        npgfj = jspp%orbnl_contract(jspp%orb_index(joa))
        ncoj = nco(l_j)*npgfj

        call re_alloc(eri, 1, nso(l_i), 1, nso(l_j), 1, nso(l_i), &
&         1, nso(l_j), name='eri', routine='calc_prescreen_eri')
        eri(:,:,:,:) = 0.0_dp

        call calc_contract_eri2(libint_data, hfx_sys%cell, hfx_sys%cell_r, &
&         ri, rj, ri, rj, npgfi, npgfj, npgfi, npgfj, &
&         l_i, l_j, l_i, l_j, ncoi, ncoj, ncoi, ncoj, &
&         ispp%orbnl_zeta(1:npgfi, ispp%orb_index(ioa)), &
&         jspp%orbnl_zeta(1:npgfj, jspp%orb_index(joa)), &
&         ispp%orbnl_zeta(1:npgfi, ispp%orb_index(ioa)), &
&         jspp%orbnl_zeta(1:npgfj, jspp%orb_index(joa)), &
&         ispp%sphi(1:ncoi,ioa:ioa+nso(l_i)-1), &
&         jspp%sphi(1:ncoj,joa:joa+nso(l_j)-1), &
&         ispp%sphi(1:ncoi,ioa:ioa+nso(l_i)-1), &
&         jspp%sphi(1:ncoj,joa:joa+nso(l_j)-1), &
&         hfx_opts, eri)

        eri_prescreen(io,jo) = 2.0_dp * max(maxval(eri), -minval(eri))

        call de_alloc(eri, name='eri', routine='calc_prescreen_eri')
      end do
    end do

    do io=hfx_sys%nuotot+1,hfx_sys%norb
      is  = hfx_sys%isa(hfx_sys%iaorb(io))
      ioa = hfx_sys%iphorb(io)
      l_i = lofio(is,ioa)
      m_i = mofio(is,ioa)

      if( m_i .ne. -l_i ) cycle

      iu = indxuo(io)
      do j=1,numh(iu)
        ind = listhptr(iu) +j
        ju  = listh(ind)
        jo  = listsc(io,iu,ju)
        js  = hfx_sys%isa(hfx_sys%iaorb(jo))
        joa = hfx_sys%iphorb(jo)
        l_j = lofio(js, joa)
        m_j = mofio(js, joa)

        if ( m_j .ne. -l_j ) cycle

        eri_prescreen(io,jo) = eri_prescreen(iu,ju)
      end do
    end do

    max_eri = max(maxval(eri_prescreen),-minval(eri_prescreen))
    max_eri = sqrt(max_eri)

    if ( Node .eq. 0 ) then
      write(*,'("calc_prescreen_eri:",1X,A,E20.8)') &
&       "max_eri_prescreen", max_eri
    end if

  end subroutine calc_prescreen_eri

! *****************************************************************************
!> \brief ...
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] libint_data: Libint data structure
!! \param[in] hfx_opts: data structure storing Hartree-Fock exchange options
!! \param[in] hfx_sys: data structure storing system information
! *****************************************************************************
  subroutine init_prescreen_eri(libint_data, hfx_opts, hfx_sys)

    use alloc         , only: de_alloc, re_alloc
    use atmfuncs      , only: lofio, mofio
    use atm_types     , only: maxn_contract, maxn_orbnl, nspecies
    use atomlist      , only: indxuo
    use parallel      , only: IOnode, Node, Nodes
    use parallelsubs
    use listsc_module , only: listsc
    use nao2gto_common
    use nao2gto_data
    use nao2gto_libint, only: Libint_t, show_libint
    use nao2gto_types , only: hfx_options_type, hfx_system_type
#ifdef MPI
    use fdf
    use mpi_siesta
    use nao2gto_parallel
#endif

    implicit none

    ! Arguments
    type(Libint_t)        , intent(inout) :: libint_data
    type(hfx_options_type), intent(in)    :: hfx_opts
    type(hfx_system_type) , intent(in)    :: hfx_sys

    ! Local variables
    integer  :: ic, jc, ncells, nelem_ij, nelem_kl
    integer  :: ind, io, ioa, ipgf, is, ishell, iu, j, jo, joa, jpgf, &
&     js, jshell, ju, l, l_j, m, m_j
    real(dp) :: max_eri
#ifdef MPI
    integer :: MPIerror, bsizedefault
    integer :: BNode, nuog, maxnumh, maxnhg, iio
    integer,  dimension(:), pointer :: numhg => null(), &
&                                      listhptrg => null(), &
&                                      listhg => null()
#endif

    ! ---------------------------------------------------------------------

#ifdef MPI
    call re_alloc(numhg, 1, hfx_sys%nuotot, name='numhg', &
&     routine='setup_hfx')
    call re_alloc(listhptrg, 1, hfx_sys%nuotot, name='listhptrg', &
&     routine='setup_hfx')

    ! Globalize numh
    do io=1,hfx_sys%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node .eq. BNode ) then
        call GlobalToLocalOrb(io, Node, Nodes, iio)
        numhg(io) = hfx_sys%numh(iio)
      endif
      call MPI_Bcast(numhg(io), 1, MPI_integer, BNode, &
&       MPI_Comm_World, MPIerror)
    enddo

    ! Build global listhptr
    listhptrg(1) = 0
    do io=2,hfx_sys%nuotot
      listhptrg(io) = listhptrg(io-1) + numhg(io-1)
    enddo

    ! Globalise listh
    maxnhg = listhptrg(hfx_sys%nuotot) + numhg(hfx_sys%nuotot)
    call re_alloc(listhg, 1, maxnhg, name='listhg', routine='setup_hfx')

    do io=1,hfx_sys%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node .eq. BNode ) then
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        do jo=1,numhg(io)
          listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = &
&           hfx_sys%listh(hfx_sys%listhptr(iio)+1: &
&             hfx_sys%listhptr(iio)+hfx_sys%numh(iio))
        enddo
      endif

      call MPI_Bcast(listhg(listhptrg(io)+1), numhg(io), MPI_integer,&
&       BNode, MPI_Comm_World, MPIerror)
    enddo
#endif

    ! Init HFX cell parameters
    ncells   = hfx_sys%nsc(1)*hfx_sys%nsc(2)*hfx_sys%nsc(3)
    nelem_ij = ncells*hfx_sys%nuotot*(hfx_sys%nuotot + 1)/2
    nelem_kl = ncells*ncells*hfx_sys%nuotot*(hfx_sys%nuotot+1)/2

    ! Prepare data for Hartree-Fock exchange
    ! to calculate ERI prescreen matrix
    call re_alloc(eri_prescreen, 1, hfx_sys%norb, 1, hfx_sys%norb, &
&         name='eri_prescreen', routine='init_prescreen_eri')
    eri_prescreen(:,:) = 0.0_dp

#ifdef MPI
    call calc_prescreen_eri(libint_data, hfx_opts, hfx_sys, &
&     maxnhg, numhg, listhptrg, listhg, max_eri)
#else
    call calc_prescreen_eri(libint_data, hfx_opts, hfx_sys, &
&     hfx_sys%maxnh, hfx_sys%numh, hfx_sys%listhptr, hfx_sys%listh, max_eri)
#endif

    ! October 2018
    call re_alloc(um_cut, 1, hfx_sys%norb, 1, hfx_sys%norb, &
&     name='um_cut', routine='init_prescreen_eri')
    um_cut(:,:) = .true.

    do io=1,hfx_sys%norb
      is = hfx_sys%isa(hfx_sys%iaorb(io))
      ioa = hfx_sys%iphorb(io)
      l = lofio(is,ioa)
      m = mofio(is,ioa)
      if ( m /= -l ) cycle

      iu = indxuo(io)
#ifdef MPI
      do j=1,numhg(iu)
        ind = listhptrg(iu) + j
        ju = listhg(ind)
        jo = listsc(io,iu,ju)
        js = hfx_sys%isa(hfx_sys%iaorb(jo))
        joa = hfx_sys%iphorb(jo)
        l_j = lofio(js,joa)
        m_j = mofio(js,joa)
        if ( m_j /= -l_j ) cycle
        um_cut(io,jo) = .false.
      end do
#else
      do j=1,hfx_sys%numh(iu)
        ind = hfx_sys%listhptr(iu) + j
        ju = hfx_sys%listh(ind)
        jo = listsc(io,iu,ju)
        js = hfx_sys%isa(hfx_sys%iaorb(jo))
        joa = hfx_sys%iphorb(jo)
        l_j = lofio(js,joa)
        m_j = mofio(js,joa)
        if ( m_j /= -l_j ) cycle
        um_cut(io,jo) = .false.
      end do
#endif
    end do

    ! Add screenfunc method (Xinming Qin)
    if ( associated(pair_dist_radii_pgf) ) then
      deallocate(pair_dist_radii_pgf)
      pair_dist_radii_pgf => null()
    end if
    allocate(pair_dist_radii_pgf(maxn_contract, maxn_contract, maxn_orbnl, &
&     maxn_orbnl, nspecies, nspecies))

    if ( associated(sfc_kind) ) then
      deallocate(sfc_kind)
      sfc_kind => null()
    end if
    allocate(sfc_kind(nspecies, nspecies))

    if ( associated(sfc_pgf) ) then
      deallocate(sfc_pgf)
      sfc_pgf => null()
    end if
    allocate(sfc_pgf(maxn_contract, maxn_contract, maxn_orbnl, maxn_orbnl, &
&     nspecies, nspecies))

    if ( associated(sfc_shell) ) then
      deallocate(sfc_shell)
      sfc_shell => null()
    end if
    allocate(sfc_shell(maxn_orbnl,maxn_orbnl,nspecies,nspecies))

    do is=1,nspecies
      do js=1,nspecies
        sfc_kind(js,is)%x(:) = 0.0_dp
        do ishell=1,maxn_orbnl
          do jshell=1,maxn_orbnl
            sfc_shell(jshell,ishell,js,is)%x(:) = 0.0_dp
            do ipgf=1,maxn_contract
              do jpgf=1,maxn_contract
                 pair_dist_radii_pgf(jpgf,ipgf,jshell,ishell,js,is)%x(:) = &
&                  0.0_dp
                 sfc_pgf(jpgf,ipgf,jshell,ishell,js,is)%x(:) = 0.0_dp
              end do
            end do
          end do
        end do
      end do
    end do

    call calc_pair_dist_radii(hfx_opts, pair_dist_radii_pgf)
    call calc_screening_functions(libint_data, hfx_opts, hfx_sys%cell, &
&     hfx_sys%cell_r, sfc_pgf, sfc_shell, pair_dist_radii_pgf)

    coeffs_kind_max0 = maxval(sfc_shell(:,:,:,:)%x(2))
    if ( IONode ) write(6,*) "max0 = ", coeffs_kind_max0

    !       Build left and right orbital pair_list
    !         max_element1 = ncells*hfx_sys%nuotot*(hfx_sys%nuotot+1)/2
    !         max_element2 = ncells**2*hfx_sys%nuotot*(hfx_sys%nuotot+1)/2
    !
    ! Note: element is not a simple type, hence we cannot use re_alloc
    if ( associated(list_ij%element) ) then
      deallocate(list_ij%element)
    end if
    allocate(list_ij%element(nelem_ij))
    if ( associated(list_kl%element) ) then
      deallocate(list_kl%element)
    end if
    allocate(list_kl%element(nelem_kl))

    write(*,*) "[DEBUG] CALL BUILD_PAIR_LIST -> LIST_IJ", nelem_ij, nelem_kl
    call build_pair_list(hfx_sys%nua, hfx_sys%na, hfx_opts, hfx_sys, &
&     eri_prescreen, max_eri, list_ij)

    write(*,*) "[DEBUG] CALL BUILD_PAIR_LIST -> LIST_KL", nelem_ij, nelem_kl
    call build_pair_list(hfx_sys%na, hfx_sys%na, hfx_opts, hfx_sys, &
&     eri_prescreen, max_eri, list_kl)

    log10_eps_schwarz = log10(hfx_opts%eps_schwarz)

    ! FIXME: I/O should be done in nao2gto_io.F90
    ! FIXME: LibFDF will soon drop MPI support
    if ( IOnode ) then
      write(*,'(A,1X,I12)')  'nelem_ij =', nelem_ij
      write(*,'(A,1X,I12)')  'nelem_kl =', nelem_kl
      write(*,'(A,1X,I12)')  'list_ij  :',list_ij%nelement
      write(*,'(A,1X,I12)')  'list_kl  :',list_kl%nelement
      write(*,'(A,1X,E20.6)') 'max_eri =', max_eri
      write(*,'(A,1X,E20.6," (log =",E20.6,")")') &
&       'eps_schwarz =', hfx_opts%eps_schwarz, log10_eps_schwarz
    end if

#ifdef MPI
    if ( IOnode ) then
      call set_bsize_NAO2GTO(Nodes, list_kl%nelement, bsizedefault)
      nao2gto_bsize = fdf_integer('blocksize_NAO2GTO', bsizedefault)
      write(*,'(A,1X,I4)') "Blocksize for ERIs distribution:", nao2gto_bsize
    end if
    call MPI_Bcast(nao2gto_bsize, 1, MPI_integer, 0, MPI_Comm_World, MPIerror)

    call de_alloc(listhg, name='listhg', routine='setup_hfx')
    call de_alloc(listhptrg, name='listhptrg', routine='setup_hfx')
    call de_alloc(numhg, name='numhg', routine='setup_hfx')
#endif

  end subroutine init_prescreen_eri

! *****************************************************************************
! *** Private routines                                                      ***
! *****************************************************************************

! ---------------------------------------------------------------------------
!  New subroutine to find shell pair list (uv) and (mn)
!
!  Coded by xmqin, 12. 06, 2016,  different from shanghui's method.
!
!  (1) Find neighbours for all atoms in supercell using rmax ,
!      Number of neighbours  : nnia    max : maxnna = 200
!      Neighbours' index: jna(nnia)
!      PBC vectors to neighbours : xij(jna)
!      Squared distances to neighbours : r2ij(jna),
!
!      Actually, I call the original subritine "neighbour()" of SISTEA
!      program to find atomic pair_list without additional definition
!      for pair_atomic_list '(type)'.
!      Loop over atoms to give xij and rij in batch !
!
!  (2) Find all shell pairs over orbital loop
!     (a) Hermite and trans symmetry : ishell_unit > jshell_unit
!     (b) NAO overlap criteria: rcut(io)+rcut(jo) > rij, io overlap with jo,
!     (c) Schwarz inequality included denity matrix
!---------------------------------------------------------------------------
!
! Useful information:
!
! (1)  To deal with periodic boundary conditions, the unit cell is 'extended'
!      on each side.
!      There are two kinds of index for images( unit cell, atoms and obitals ) :
!      one is in the PBC supercell, the other is in the normal (or auxiliary) supercell.
!
!  PBC unit cell : Wigner-Seitz cell A(i) = -1/2a(i), 1/2a(i) , a, Ai : lattice vector.
!  PBC supercell : To extend Wigner cell on each side of direction: -nsc(i)*A/2, nsc(i)*A/2
!  Normal supercell : Direct supercell : A(i) = 1, nsc(i)*a !
!
!  The index of images should be shift between pbc supercell and normal supercell
!  because we usually use normal index.
!
!  (1) In order to take into account the 8-fold symmetry under pbc :
!
!   (u0v[R]|m[R']n[R"])        =   (v0u[-R]|m[R'-R]n[R"-R])
! = (u0v[R]|n[R"]m[R'])        =   (v0u[-R]|n[R"-R]m[R'-R])
! = (m0n[R"-R']|u[-R']v[R-R']) =   (n0m[R'-R"]|u[-R"]v[R-R"])
! = (m0n[R"-R']|v[R-R']u[-R']) =   (n0m[R'-R"]|v[R-R"]u[-R"])
!
!  The first index must be always in unit cell, so we should transfer R to 0.
!
!   pair_list by iterating in the following way.
!   we should transfer R to 0
!
!   DO iatom=1,natom
!      iatom_unit = mod(iatom-1, unit_atoms)+1
!
!      DO jatom=iatom, negbours(iatom)
!         jatom_pbc = ind1 (ind2(jatom)-ind2(iatom)+ind2(iatom_unit))
!         jatom_unit_pbc = mod(jatom-1, unit_atoms)+1
!         if (jatom_unit<iatom_unit) CYCLE
!
!         atom_ij = ncells*(iatom_unit)*(iatom_unit-1)/2
!                  +mod(jatom_pbc,unit_atoms)*iatom_unit+jatom_unit_pbc
!
!  atom_ij is the uv or mn pair list index gives u[0]v[R] or
!  m[0]n[R],  it is not equal to number of loop.
!
!  Here "_pbc" denotes images' index in pbc supercell.
!  v[0] means that we choose v[R] as reference unitcell.
!  then u, m, n must be shift to extended supercell related to v[R].
!
!!
!! if (katom+latom<=iatom+jatom)  then
!!   if ( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
! or
!----------------------------------------------------------------------------
! integer nia: Atom index of orbital (u)s into the loop .
!              To build (uv), nia = nua  for u
!                       (mn), nia = na   for m

! integer norb         : Number of orbitals in supercell
! real*8  scell(3,3)   : Supercell vectors SCELL(IXYZ,IVECT)
! integer nsc(3)       : Num. of unit cells in each supercell direction
! real*8  xa(3,na)     : Atomic positions in cartesian coordinates
! integer lasto(0:na)  : Last orbital of each atom in array iphorb
! integer iphorb(norb) : Orbital index of each orbital in its atom,
!                       where no=lasto(na)
! integer iphorb(norb)
! integer isa(na)      : Species index of each atom

!---------------------------------------------------------------------------
! *****************************************************************************
!> \brief ...
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] libint_data: Libint data structure
! *****************************************************************************
  subroutine build_pair_list(nia, na, hfx_opts, hfx_sys, &
                 eri_prescreen, max_eri, list_orb)

    use precision    , only: dp
    use alloc        , only: de_alloc, re_alloc
    use atmfuncs     , only: lofio, mofio, rcut
    use atomlist     , only: indxuo, rmaxo, lasto
    use neighbour    , only: jna=>jan, xij, r2ij, maxna=>maxnna, mneighb
    use sorting      , only: iorder, order, ordix
    use nao2gto_data , only: subshell
    use nao2gto_index, only: indexsc
    use nao2gto_pbc  , only: trans_pbc
    use nao2gto_types, only: hfx_options_type, hfx_system_type, pair_list_type

    implicit none

    ! Arguments
    integer               , intent(in)  :: nia, na
    type(hfx_options_type), intent(in)  :: hfx_opts
    type(hfx_system_type) , intent(in)  :: hfx_sys
    real(dp)              , intent(in)  :: &
&     eri_prescreen(hfx_sys%norb,hfx_sys%norb)
    real(dp)              , intent(in)  :: max_eri
    type(pair_list_type)  , intent(out) :: list_orb

    ! Local variables
    integer  :: ncells, n_shells, nnia
    integer  :: ia, io, ioa, is, i, ja, jn, jo, joa, js, iuo, juo, jos, iua
    integer  :: l_i, l_j, m_i, m_j, ishell, jshell, &
&     ishell_unit, jshell_unit, unit_shells
    real(dp) :: rij, r2, ri(3), rj(3), r_temp(3), r_pbc(3)
    integer, dimension(:), pointer :: index => null()

    ! -------------------------------------------------------------------------

    ! ......................
    !  rmaxo = max rcut of all orbitals
    !      if(node.eq.0) then
    !       write(6,'(a)') 'To find the neighbours of an atom in supercell.'
    !       write(6,'(a,f12.6,a)') "using rmaxo= ", rmaxo+rmaxkb, " bohr."
    !      end if

    ! Initialize neighb subroutine

    ! Allocate local memory
    unit_shells = subshell(hfx_sys%nuotot)
    ncells = hfx_sys%norb / hfx_sys%nuotot
    call mneighb(hfx_sys%cell, 2.0_dp*rmaxo, na, hfx_sys%xa, 0, 0, nnia )
    call re_alloc(index, 1, maxna, name='index', routine='build_pair_list')

    list_orb%nelement = 0
    do ia = 1,nia ! the last cell 5*5*5 scell
      iua = mod(ia-1, hfx_sys%nua) + 1
      ri(1:3) = hfx_sys%xa(1:3,ia)
      call mneighb(hfx_sys%cell, 2.0_dp*rmaxo, na, hfx_sys%xa, ia, 0, nnia)

      ! Sort by by distance
      call ordix(r2ij, 1, nnia, index)
      call iorder(jna, 1, nnia, index)
      call order(r2ij, 1, nnia, index)
      call order(xij, 3, nnia, index)

      do io=lasto(ia-1)+1,lasto(ia)
        ioa = hfx_sys%iphorb(io)
        is  = hfx_sys%isa(hfx_sys%iaorb(io))
        l_i = lofio(is,ioa)
        m_i = mofio(is,ioa)

        if ( m_i .ne. -l_i ) cycle

        ishell = subshell(io)
        iuo = mod(io-1, hfx_sys%nuotot) + 1
        ishell_unit = subshell(iuo)

        !  Find neighbor atoms using rmaxo = max(rcut(natom))
        do jn=1,nnia
          ja = jna(jn)
          r_temp(1:3) = hfx_sys%xa(1:3,ja) - hfx_sys%xa(1:3,ia)
          ! xmqin for debug
          call trans_pbc(r_temp, hfx_sys%cell, hfx_sys%cell_r, r_pbc)
          r2 = (r_pbc(1)**2) + (r_pbc(2)**2) + (r_pbc(3)**2)
          rj = ri + r_pbc
          rij = sqrt(r2)

          do jo=lasto(ja-1)+1,lasto(ja)
            joa = hfx_sys%iphorb(jo)
            js  = hfx_sys%isa(hfx_sys%iaorb(jo))
            l_j = lofio(js,joa)
            m_j = mofio(js,joa)

            if ( m_j .ne. -l_j ) cycle

            jos = indexsc(io, iuo, jo)
            juo = mod(jos-1, hfx_sys%nuotot) + 1
            jshell = subshell(jos)
            jshell_unit = subshell(juo)

            if ( jshell_unit .le. ishell_unit ) then
              if ( (rcut(is,ioa) + rcut(js,joa)) .le. rij ) cycle

              if ( sqrt(abs(eri_prescreen(io,jo)))*max_eri .ge. &
                   hfx_opts%eps_pairlist ) then
                n_shells = ncells*(ishell_unit)*(ishell_unit-1)/2 &
&                        + ((jshell-1)/unit_shells)*ishell_unit + jshell_unit

                list_orb%nelement = list_orb%nelement + 1
                i = list_orb%nelement
                list_orb%element(i)%pair(1) = io
                list_orb%element(i)%pair(2) = jo
                list_orb%element(i)%nl_index = n_shells
                list_orb%element(i)%r1(1:3) = ri(1:3)
                list_orb%element(i)%r2(1:3) = rj(1:3)
                list_orb%element(i)%dist2 = r2
              end if
            end if
          end do   ! jo
        end do   ! jn
      end do   ! io
    end do   ! ia

    call de_alloc(index, name='index', routine='build_pair_list')

  end subroutine build_pair_list


  ! ***************************************************************************
  !> \brief Calculates distances between pairs of orbitals
  ! ***************************************************************************
  subroutine calc_pair_dist_radii(hfx_opts, radii_pgf)

    use alloc
    use atm_types
    use nao2gto_types
    use nao2gto_utils

    implicit none

    ! Arguments
    type(hfx_options_type), intent(in) :: hfx_opts
    type(hfx_screen_coeff_type), dimension(:,:,:,:,:,:), pointer :: radii_pgf

    ! Local variables
    integer :: handle, i, ikind, ipgf, isave, iset, jkind, jpgf, jset, la, &
&     lb, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb, stat
    integer :: ia, ioa, io, iu, is, ja, joa, jo, ju, js, ind, l_i, l_j, &
&     m_i, m_j, ncoi, ncoj, npgfi, npgfj, ishell, jshell
    integer, dimension(:), pointer :: npgfa => null(), npgfb => null(), &
&     nsgfa => null(), nsgfb => null()
    real(dp) :: cutoff, table(2,0:100), ff, max_contraction_a, &
&     max_contraction_b, prefactor, R1, R_max, ra(3), rab(3), rab2, radius, &
&     rap(3), rb(3), rp(3), x(2), zetp
    real(dp), dimension(:), pointer :: set_radius_a, set_radius_b
    real(dp), dimension(:,:), pointer :: rpgfa => null(), rpgfb => null(), &
&     sphi_a => null(), sphi_b => null(), zeta => null(), zetb => null()
    type(species_info), pointer :: ispp => null(), jspp => null()

    ! -------------------------------------------------------------------------

    ra = 0.0_dp
    rb = 0.0_dp
    table = 0.0_dp

    do is=1,nspecies
      ispp => species(is)
      do js=1,nspecies
        jspp => species(js)

        ishell = 0
        do io=1,ispp%norbs
           l_i = ispp%orb_l(io)
           m_i = ispp%orb_m(io)
           if ( m_i .ne. -l_i ) cycle     ! Not a normal orbital
           ishell = ishell+1
           npgfi = ispp%orbnl_contract(ishell)
           ncoi = nco(l_i)*npgfi

          jshell = 0
          do jo=1,jspp%norbs
            l_j = jspp%orb_l(jo)
            m_j = jspp%orb_m(jo)
            if ( m_j .ne. -l_j ) cycle     ! Not a normal orbital
            jshell = jshell+1
            npgfj = jspp%orbnl_contract(jshell)
            ncoj = nco(l_j)*npgfj

            do ipgf=1,npgfi
              do jpgf=1,npgfj
                 radius = ispp%pgf_radius(ipgf,ishell) + &
&                         jspp%pgf_radius(jpgf,jshell)
                do i=0,100
                  rb(1) = 0.0_dp + 0.01_dp * radius*i
                  R_max = 0.0_dp
                  zetp = ispp%orbnl_zeta(ipgf,ishell) + &
&                        jspp%orbnl_zeta(jpgf,jshell)
                  ff = jspp%orbnl_zeta(jpgf,jshell)/zetp
                  rab = 0.0_dp
                  rab(1) = rb(1)
                  rab2 = rb(1)**2
                  prefactor = exp(-ispp%orbnl_zeta(ipgf,ishell)*ff*rab2)
                  rap(:) = ff*rab(:)
                  rp(:) = ra(:) + rap(:)
                  rb(:) = ra(:) + rab(:)
                  cutoff = 1.0_dp

                  R1 = exp_radius_very_extended(l_i, l_i, l_j, l_j, &
&                   ra=ra, rb=rb, rp=rp, zetp=zetp, eps=1.0e-6_dp, &
&                   prefactor=prefactor, cutoff=cutoff,epsin=1.0e-12_dp)
                  R_max = MAX(R_max,R1)

                  table(1,i) = rb(1)
                  table(2,i) = R_max
                end do

                call run_powell_optimize(table,x,0.0_dp)
                radii_pgf(jpgf,ipgf,jshell,ishell,js,is)%x = x
              end do !jpgf
            end do !ipgf
          end do
        end do
      end do
    end do

  end subroutine calc_pair_dist_radii

  ! ***************************************************************************
  !> \brief Calculates screening functions
  ! ***************************************************************************
  subroutine calc_screening_functions(libint_data, hfx_opts, cell, rcell,  &
&   coeffs_pgf, coeffs_set, radii_pgf)

    use alloc
    use atm_types
    use nao2gto_contract
    use nao2gto_libint, only: Libint_t
    use nao2gto_primitive
    use nao2gto_types

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data
    type(hfx_options_type), intent(in) :: hfx_opts
    real(dp), intent(in)  :: cell(3,3), rcell(3,3)
    type(hfx_screen_coeff_type), dimension(:,:,:,:), pointer :: coeffs_set
    type(hfx_screen_coeff_type), dimension(:,:,:,:,:,:), pointer :: &
&     radii_pgf, coeffs_pgf

    ! Local variables
    integer :: handle, i, ikind, ipgf, iset, jkind, jpgf, jset, &
&     la, lb, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb, stat
    integer :: ia, ioa, io, iu, is, ja, joa, jo, ju, js, ind, l_i, l_j, &
&     m_i, m_j, ncoi, ncoj, npgfi, npgfj, ishell, jshell
    integer, dimension(:), pointer :: la_max => null(), la_min => null(), &
&     lb_max => null(), lb_min => null(), npgfa => null(), npgfb => null(), &
&     nsgfa => null(), nsgfb => null()
    real(dp) :: zeta_i, zeta_j
    real(dp) :: table(2,0:100), kind_radius_a, kind_radius_b, &
&     max_contraction_a, max_contraction_b, max_val, max_val_temp, &
& R1, ra(3), radius, rb(3), x(2), max_val2, max_val_temp2
    real(dp), dimension(:,:,:,:), pointer :: eri => null()
    real(dp), dimension(:), pointer :: set_radius_a => null(), &
&     set_radius_b => null()
    real(dp), dimension(:, :), pointer :: rpgfa => null(), rpgfb => null(), &
&     sphi_a => null(), sphi_b => null(), zeta => null(), zetb => null()
    type(hfx_screen_coeff_type), dimension(:,:), pointer :: tmp_R_1 => null()
    type(species_info), pointer  :: ispp => null(), jspp => null()

    ! -------------------------------------------------------------------------

    ra = 0.0_dp
    rb = 0.0_dp
    table = 0.0_dp

    do is=1,nspecies
       ispp => species(is)
      do js=1,nspecies
         jspp => species(js)

        ishell = 0
        do io=1,ispp%norbs
          l_i = ispp%orb_l(io)
          m_i = ispp%orb_m(io)
          if ( m_i /= -l_i ) cycle     ! Not a normal orbital
          ishell = ishell+1
          npgfi = ispp%orbnl_contract(ishell)
          ncoi  = nco(l_i)*npgfi
          max_contraction_a = maxval((/(sum(abs(ispp%sphi(1:ncoi,i))), &
&           i=io,io+nso(l_i)-1)/))

          jshell = 0
          do jo = 1,jspp%norbs
            l_j = jspp%orb_l(jo)
            m_j = jspp%orb_m(jo)
            if ( m_j /= -l_j ) cycle     ! Not a normal orbital
            jshell = jshell+1
            npgfj = jspp%orbnl_contract(jshell)
            ncoj  = nco(l_j)*npgfj
            max_contraction_b =  maxval((/(sum(abs(jspp%sphi(1:ncoj,i))), &
&             i=jo,jo+nso(l_j)-1)/))
            do ipgf = 1, npgfi
              zeta_i = ispp%orbnl_zeta(ipgf,ishell)
              do jpgf = 1, npgfj
                zeta_j = jspp%orbnl_zeta(jpgf,jshell)
                radius=ispp%pgf_radius(ipgf,ishell)+jspp%pgf_radius(jpgf,jshell)

                do i=0,100
                  rb(1) = 0.0_dp + real(i,dp) * 0.01_dp * radius
                  max_val = 0.0_dp
                  max_val_temp = 0.0_dp

                  R1 = max(0.0_dp, &
&                   radii_pgf(jpgf,ipgf,jshell,ishell,js,is)%x(1)*rb(1)**2 &
                       + radii_pgf(jpgf,ipgf,jshell,ishell,js,is)%x(2))

                  call calc_primitive_screen(libint_data, ra, rb, ra, rb, &
&                   zeta_i, zeta_j, zeta_i, zeta_j, l_i, l_j, l_i ,l_j, &
&                   max_val_temp, hfx_opts, R1, R1)
                  max_val = max(max_val, max_val_temp)
                  max_val = sqrt(max_val)

                  max_val = max_val * max_contraction_a * max_contraction_b
                  table(1,i) = rb(1)
                  if ( max_val == 0.0_dp ) then
                    table(2,i) = powell_min_log
                  else
                    table(2,i) = LOG10((max_val))
                  end if
                end do

                call run_powell_optimize(table,x,powell_min_log)
                coeffs_pgf(jpgf,ipgf,jshell,ishell,js,is)%x = x

              end do
            end do
          end do
        end do
      end do
    end do

    ra = 0.0_dp
    rb = 0.0_dp
    table = 0.0_dp

    do is=1,nspecies
      ispp => species(is)
      do js = 1,nspecies
        jspp => species(js)

        ishell = 0
        do io = 1,ispp%norbs
           l_i = ispp%orb_l(io)
           m_i = ispp%orb_m(io)
           if (m_i .ne. -l_i) cycle     ! Not a normal orbital
           ishell = ishell+1
           npgfi = ispp%orbnl_contract(ishell)
           ncoi  = nco(l_i)*npgfi

           max_contraction_a = maxval((/(sum(abs(ispp%sphi(1:ncoi,i))), &
&            i=io,io+nso(l_i)-1)/))

          jshell = 0
          do jo = 1,jspp%norbs
            l_j = jspp%orb_l(jo)
            m_j = jspp%orb_m(jo)
            if ( m_j .ne. -l_j ) cycle     ! Not a normal orbital
             jshell = jshell+1
             npgfj = jspp%orbnl_contract(jshell)
             ncoj  = nco(l_j)*npgfj
             max_contraction_b = maxval((/(sum(abs(jspp%sphi(1:ncoj,i))), &
&              i=jo,jo+nso(l_j)-1)/))
             radius = ispp%shell_radius(ishell)+jspp%shell_radius(jshell)

             tmp_R_1 => radii_pgf(:,:,jshell,ishell,js,is)

            do i=0,100
              rb(1) = 0.0_dp + real(i,dp) * 0.01_dp * radius
              max_val = 0.0_dp
              max_val_temp = 0.0_dp

              call re_alloc(eri, 1, nso(l_i), 1, nso(l_j), 1, nso(l_i), &
&               1, nso(l_j), name='eri', routine='calc_screening_functions')
              eri = 0.0_dp
              call calc_contract_eri2(libint_data, cell, rcell, ra, rb, ra, rb, &
                npgfi, npgfj, npgfi, npgfj, l_i, l_j, l_i, l_j, &
                ncoi, ncoj, ncoi, ncoj, ispp%orbnl_zeta(1:npgfi, ishell), &
                jspp%orbnl_zeta(1:npgfj, jshell), &
                ispp%orbnl_zeta(1:npgfi,ishell),  &
                jspp%orbnl_zeta(1:npgfj, jshell), &
                ispp%sphi(1:ncoi,io:io+nso(l_i)-1), &
                jspp%sphi(1:ncoj,jo:jo+nso(l_j)-1), &
                ispp%sphi(1:ncoi,io:io+nso(l_i)-1), &
                jspp%sphi(1:ncoj,jo:jo+nso(l_j)-1), &
                hfx_opts, eri)

              max_val = max(max_val, maxval(dabs(eri)))
              max_val = 2.0_dp*max_val
              max_val = DSQRT(max_val)

              table(1,i) = rb(1)
              if ( max_val == 0.0_dp ) then
                table(2,i) = powell_min_log
              else
                table(2,i) = LOG10((max_val))
              end if

              call de_alloc(eri, name='eri', routine='calc_screening_functions')

            end do
            call run_powell_optimize(table,x,powell_min_log)
            coeffs_set(jshell,ishell,js,is)%x = x
          end do
        end do

      end do
    end do

  end subroutine calc_screening_functions

  ! ***************************************************************************
  !> \brief Driver for the Powell minimizer
  !!
  !! This routine drives the powell minimizer for the data found in a table
  !! storing function values. It constructs an approximate upper bound of
  !! the fitted function.
  !!
  !! This is done in two steps: first, we compute the symmetric weight
  !! to get a good unique initial guess, then we restart for the asymmetric
  !! one. The asymmetric function appears to have several local minima.
  !! Depending on the data to fit the loop over k can make the switch
  !! gradual, but there is seemingly not much need for it.
  !!
  !! \author Xinming Qin
  !!
  !! \par History
  !!      10.2018 Created [Xinming Qin]
  !!      11.2018 Refactored and further documented [Yann Pouillon]
  !!
  !! \param[in] table: data to fit
  !! \param[out] x: fitting coefficients of the form (x(1)*table(1)**2+x(2))
  !! \param[in] fmin: only values of table(2) larger than fmin are considered
  ! ***************************************************************************
  subroutine run_powell_optimize(table, x, fmin)

    use nao2gto_powell

    implicit none

    ! Arguments
    real(dp), intent(in)  :: table(2,0:100)
    real(dp), intent(out) :: x(2)
    real(dp), intent(in)  :: fmin

    ! Local variables
    integer :: i, k
    real(dp) :: f, large_weight, small_weight, weight
    type(opt_state_type) :: opt_state

    ! -------------------------------------------------------------------------

    ! Initial values
    x(1) = 0.0_dp
    x(2) = 0.0_dp

    do k=0,4,4

      small_weight=1.0_dp
      large_weight=small_weight*(10.0_dp**k)

      ! Init optimization run
      opt_state%state = 0
      opt_state%nvar = 2
      opt_state%iprint = 3
      opt_state%unit = 6
      opt_state%maxfun = 100000
      opt_state%rhobeg = 0.1_dp
      opt_state%rhoend = 0.000001_dp

      do

        ! Compute function
        if ( opt_state%state == 2 ) then
          opt_state%f = 0.0_dp
          do i=0,100
            f = x(1)*table(1,i)**2 +  x(2)
            if ( f > table(2,i) ) then
              weight = small_weight
            else
              weight = large_weight
            end if
            if ( table(2,i) > fmin ) then
              opt_state%f = opt_state%f + weight * (f-table(2,i))**2
            end if
          end do
        end if

        if ( opt_state%state == -1 ) exit

        call powell_optimize(opt_state%nvar, x, opt_state)
      end do

      ! Clean-up the mess
      opt_state%state = 8
      call powell_optimize(opt_state%nvar, x, opt_state)

     end do

  end subroutine run_powell_optimize

  ! ***************************************************************************
  !> \brief Build a symmetric upper triangle matrix
  !!
  !! Given a 2d index pair, this function returns a 1d index pair for a
  !! symmetric upper triangle NxN matrix. The compiler should inline this
  !! function, therefore it appears in several modules.
  !!
  !! \author Manuel Guidon
  !!
  !! \par History
  !!      03.2009 Created [Manuel Guidon]
  !!      11.2018 Refactored and further documented [Yann Pouillon]
  !!
  !! \param i: index on the first array dimension
  !! \param j: index on the second array dimension
  !! \param n: matrix size
  !! \return NxN matrix
  ! ***************************************************************************
  pure function get_1d_idx(i, j, n)

    implicit none

    ! Arguments
    integer, intent(in) :: i, j
    integer(int_8), intent(in) :: n

    ! Local variables
    integer(int_8) :: get_1d_idx
    integer(int_8) :: min_ij

    ! -------------------------------------------------------------------------

    min_ij = min(i, j)
    get_1d_idx = min_ij*N + max(i, j) - (min_ij-1)*min_ij/2 - n

  end function get_1d_idx

end module nao2gto_prescreen
