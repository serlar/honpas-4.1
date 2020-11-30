! *** Module: nao2gto_data ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Management of NAO2GTO data
!!
!! This module takes care of managing the initialization, tracking, and
!! destruction, of NAO2GTO data.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.2018 Created [Yann Pouillon]
! *****************************************************************************
module nao2gto_data

  use precision     , only: dp
  use nao2gto_libint, only: Libint_t
  use nao2gto_types

  implicit none

  private

  integer, public :: hfx_call_counter = 0

  real(dp), public :: coeffs_kind_max0 = 0.0_dp
  real(dp), public :: log10_eps_schwarz = 0.0_dp

  type(eri_link_type), pointer, public :: eri_head => null()
  type(eri_link_type), pointer, public :: eri_ptr => null()
  type(eri_link_type), pointer, public :: eri_tail => null()

  real(dp), pointer, public :: eri_prescreen(:,:) => null()

  type(Libint_t)        , public :: hfx_libint
  type(hfx_options_type), public :: hfx_options
  type(hfx_system_type) , public :: hfx_system

  type(pair_list_type), public :: list_ij, list_kl

  integer, pointer, public :: subshell(:) => null()

  type(hfx_screen_coeff_type), dimension(:,:,:,:,:,:), pointer, public :: &
    pair_dist_radii_pgf => null(), sfc_pgf => null()
  type(hfx_screen_coeff_type), dimension(:,:,:,:), pointer, public :: &
    sfc_shell => null()
  type(hfx_screen_coeff_type), dimension(:,:), pointer, public :: &
    sfc_kind => null()

  logical, pointer, public :: um_cut(:,:) => null()

  public :: &
&   nao2gto_libint_end, &
&   nao2gto_libint_init, &
&   nao2gto_system_end, &
&   nao2gto_system_init, &
&   nao2gto_system_update_cell

contains

  ! ***************************************************************************
  !> \brief Terminates the NAO2GTO data structures associated to Libint
  !!
  !! \author Yann Pouillon
  !!
  !! \param[in,out] libint_data: Libint data structure
  ! ***************************************************************************
  subroutine nao2gto_libint_end(libint_data)

    use nao2gto_index,  only: nao2gto_index_free
    use nao2gto_libint, only: Libint_t, terminate_libint
    use nao2gto_utils , only: deallocate_md_ftable

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data

    ! -------------------------------------------------------------------------

    ! Get rid of internal tables
    call deallocate_md_ftable()
    call nao2gto_index_free()
    if ( associated(subshell) ) then
      deallocate(subshell)
      nullify(subshell)
    endif

    ! Get rid of the Libint data structure
    call terminate_libint(libint_data)

  end subroutine nao2gto_libint_end

  ! ***************************************************************************
  !> \brief Initializes the NAO2GTO data structures associated to Libint
  !!
  !! \author Yann Pouillon
  !!
  !! \param[out] libint_data: Libint data structure
  ! ***************************************************************************
  subroutine nao2gto_libint_init(libint_data, l_max, hfx_sys)

    use alloc         , only: re_alloc
    use atmfuncs      , only: lofio, mofio
    use nao2gto_index , only: nao2gto_index_init
    use nao2gto_libint, only: Libint_t, initialize_libint
    use nao2gto_types , only: hfx_system_type
    use nao2gto_utils , only: init_md_ftable

    implicit none

    ! Arguments
    integer           , intent(in)  :: l_max
    type(Libint_t)    , intent(inout) :: libint_data
    type(hfx_system_type), intent(in)  :: hfx_sys

    ! Local variables
    integer :: i, io, ioa, is, l, m

    ! -------------------------------------------------------------------------

    ! Supercell orbital initialisation
    call nao2gto_index_init(hfx_sys%nsc, hfx_sys%nuotot)

    ! Gamma function table and Libint initialisation
    call init_md_ftable(4*l_max)
    call initialize_libint(libint_data, l_max)

    ! Orbital shell initialisation
    call re_alloc(subshell, 1, hfx_sys%norb, &
&     name='subshell', routine='nao2gto_libint_init')
    subshell(:) = 0

    ! 2l+1 orbitals form a shell
    i = 0
    do io=1,hfx_sys%norb
      is = hfx_sys%isa(hfx_sys%iaorb(io))
      ioa = hfx_sys%iphorb(io)
      l = lofio(is,ioa)
      m = mofio(is,ioa)
      if (m .eq. -l ) then
        i = i + 1
      endif
      subshell(io) = i
    enddo

  end subroutine nao2gto_libint_init

  ! ***************************************************************************
  !> \brief Terminates the NAO2GTO data structures associated to SIESTA
  !!
  !! \author Yann Pouillon
  !!
  !! \param[in,out] hfx_sys: information about the system
  ! ***************************************************************************
  subroutine nao2gto_system_end(hfx_sys)

    use nao2gto_types, only: hfx_system_type

    implicit none

    ! Arguments
    type(hfx_system_type), intent(inout) :: hfx_sys

    ! -------------------------------------------------------------------------

    ! Reset cell parameters
    hfx_sys%cell(:,:) = 0.0_dp
    hfx_sys%cell_r(:,:) = 0.0_dp

    ! Make system information unusable
    nullify(hfx_sys%maxnh)
    nullify(hfx_sys%na)
    nullify(hfx_sys%norb)
    nullify(hfx_sys%nspin)
    nullify(hfx_sys%nua)
    nullify(hfx_sys%nuo)
    nullify(hfx_sys%nuotot)
    nullify(hfx_sys%iaorb)
    nullify(hfx_sys%indxua)
    nullify(hfx_sys%iphorb)
    nullify(hfx_sys%isa)
    nullify(hfx_sys%listh)
    nullify(hfx_sys%listhptr)
    nullify(hfx_sys%nsc)
    nullify(hfx_sys%numh)
    nullify(hfx_sys%xa)

  end subroutine nao2gto_system_end

  ! ***************************************************************************
  !> \brief Initializes the NAO2GTO data structures associated to SIESTA
  !!
  !! \author Yann Pouillon
  !!
  !! \param[out] hfx_sys: information about the Hartree-Fock system data
  !! \param[in] maxnh: maximum number of non-zero H matric elements
  !! \param[in] na: number of atoms in the supercell
  !! \param[in] norb: number of orbitals in the supercell
  !! \param[in] nspin: number of spìn degrees of freedom
  !! \param[in] nua: number of atoms in the unit cell
  !! \param[in] nuo: number of orbitals local to node
  !! \param[in] nuotot: number of orbitals in the unit cell
  !! \param[in] iaorb: atomic index of each orbital
  !! \param[in] indxua: ...
  !! \param[in] iphorb: orbital index of each orbital in its atom
  !! \param[in] isa: ...
  !! \param[in] listh: ...
  !! \param[in] listhptr: ...
  !! \param[in] nsc: number of supercells in each direction (diagonal of mscell)
  !! \param[in] numh: ...
  !! \param[in] ucell: unit cell in real space
  !! \param[in] xa: positions of the atoms in the unit cell
  ! ***************************************************************************
  subroutine nao2gto_system_init(hfx_sys, maxnh, na, norb, nspin, nua, &
&                nuo, nuotot, iaorb, indxua, iphorb, isa, &
&                listh, listhptr, nsc, numh, ucell, xa)

    use nao2gto_types, only: hfx_system_type

    implicit none

    ! Arguments
    integer, target, intent(in)  :: maxnh, na, norb, nspin, nua, nuo, nuotot
    integer, target, intent(in)  :: iaorb(norb), indxua(na), iphorb(norb), &
&     isa(na), listh(maxnh), listhptr(nuo), nsc(3), numh(nuo)
    real(dp), target, intent(in) :: ucell(3,3), xa(3,na)
    type(hfx_system_type), intent(out) :: hfx_sys

    ! Local variables
    integer  :: i, io, ioa, is, l, m, ncells
    integer  :: ic, jc

    ! -------------------------------------------------------------------------

    ! Make system information available to Hartree-Fock exchange routines,
    ! avoiding data copy
    hfx_sys%maxnh => maxnh
    hfx_sys%na => na
    hfx_sys%norb => norb
    hfx_sys%nspin => nspin
    hfx_sys%nua => nua
    hfx_sys%nuo => nuo
    hfx_sys%nuotot => nuotot
    hfx_sys%iaorb => iaorb
    hfx_sys%indxua => indxua
    hfx_sys%iphorb => iphorb
    hfx_sys%isa => isa
    hfx_sys%listh => listh
    hfx_sys%listhptr => listhptr
    hfx_sys%nsc => nsc
    hfx_sys%numh => numh
    hfx_sys%xa => xa

    call nao2gto_system_update_cell(hfx_sys, nsc, ucell)

  end subroutine nao2gto_system_init

  ! ***************************************************************************
  !> \brief Updates the supercell used by NAO2GTO data structures
  !!
  !! \author Yann Pouillon
  !!
  !! \param[inout] hfx_sys: information about the Hartree-Fock system data
  !! \param[in] nsc: number of supercells in each direction (diagonal of mscell)
  !! \param[in] ucell: unit cell in real space
  ! ***************************************************************************
  subroutine nao2gto_system_update_cell(hfx_sys, nsc, ucell)

    use nao2gto_types, only: hfx_system_type

    implicit none

    ! Arguments
    integer, intent(in)  :: nsc(3)
    real(dp), intent(in) :: ucell(3,3)
    type(hfx_system_type), intent(inout) :: hfx_sys

    ! Local variables
    integer  :: ic, jc

    ! -------------------------------------------------------------------------

    ! Determine supercell both in real space and reciprocal space
    hfx_sys%cell(:,:) = 0.0_dp
    hfx_sys%cell_r(:,:) = 0.0_dp
    do ic=1,3
      do jc=1,3
        hfx_sys%cell(jc,ic) = ucell(jc,ic) * nsc(ic)
      enddo
    enddo
    call reclat(hfx_sys%cell, hfx_sys%cell_r, 0)

  end subroutine nao2gto_system_update_cell

end module nao2gto_data
