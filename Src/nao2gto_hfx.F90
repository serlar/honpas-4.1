! *** Module: nao2gto_hfx ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Hartree-Fock exchange term of the Hamiltonian
!!
!! In the SIESTA method, the density matrix is stored using the same sparse
!! pattern as the Hamiltonian. This treatment is robust in pure DFT where the
!! effective pure-DFT potential is only determined by the electron density
!! and its gradients.
!!
!! \f$\rho(r) = \sum_{u,v}[Duv*\psi_u(r-R_u)*\psi_v(r-R_v)]\f$,
!! \f$V_{eff}[\rho(r)] /= 0\f$ only when
!! \f$\psi_u(r-R_u)\f$ overlaps with \f$\psi_v(r-R_v)\f$.
!!
!! To calculate \f$rho\f$ and \f$V_{eff}\f$, we only need to store a sparse
!! subset of \f$D_{uv}\f$ that u overlaps with v even though the actual
!! \f$D_{uv}\f$ is unknown before SCF.
!!
!! However, HFX is dependent with the density matrix of real space:
!! \f$\rho(r,r^') = \sum_{u,v} \left[ D_{uv}*\psi_u(r-R_u)*\psi_v(r^'-R_v)
!! \right]\f$
!! \f$\rho(r,r^') \approx exp[-a(r-r^'-R_u-R_v)]\f$, a property of
!! \f$\sqrt(E_g)\f$ in semiconductors.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 11.2016 Edited [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \bug Possible errors from sparse Hamiltonian and density matrix!
! *****************************************************************************
module nao2gto_hfx

  use nao2gto_common

  implicit none

  private

  public :: setup_hfx

contains

  ! ***************************************************************************
  !> \brief Computes the Hartree-Fock exchange part of the Hamiltonian
  !!
  !! Parallel update HFX(u,v) to Hamiltonian,
  !! u: Local PAOs of unitcell per node,
  !! v: Global PAOs of supercell
  !!
  !! HFX(u,v) = (um|vn) * Dm (m,n)
  !! HFX(u,v) 2D structure ---> Hmat(ind) 1D sparse structure
  !!
  !! Hmat(pure DFT)<--- HFX (this subroutine)----> HFX (HSE06)
  !!
  !! \par History
  !!      - 10.2013 Created [Xinming Qin]
  !!      - 11.2016 Edited [Xinming Qin]
  !!      - 01.2018 Replaced subroutine arguments by direct use of SIESTA
  !!        variables [Yann Pouillon]
  !!
  !! \param[in,out] libint_data: initialized Libint data structure
  !! \param[in] hfx_optdata: initialized Hartree-Fock options data structure
  !! \param[in] hfx_sysdata: initialized Hartree-Fock system data structure
  !! \param[in.out] Hmat: initialized hamiltonian where to add the
  !!                      HFX contribution
  ! ***************************************************************************
  subroutine setup_hfx(libint_data, hfx_optdata, hfx_sysdata, Hmat)

    use parallel,        only: IOnode, Node, Nodes
    use parallelsubs,    only: GetNodeOrbs, GlobalToLocalOrb, &
&                              LocalToGlobalOrb, WhichNodeOrb
    use alloc,           only: re_alloc, de_alloc
#ifdef MPI
    use mpi_siesta
#endif

    ! Pre-BSC input parameters, now module variables
    use m_energies,      only: Exc
    use sparse_matrices, only: Dscf, maxnh

    use nao2gto_data  , only: eri_prescreen, hfx_call_counter
    use nao2gto_dm
    use nao2gto_libint, only: Libint_t
    use nao2gto_types , only: hfx_options_type, hfx_system_type

    implicit none

    ! Arguments
    type(Libint_t), intent(inout)      :: libint_data
    type(hfx_options_type), intent(in) :: hfx_optdata
    type(hfx_system_type), intent(in)  :: hfx_sysdata
    real(dp), intent(inout) :: Hmat(hfx_sysdata%maxnh,hfx_sysdata%nspin)

    ! Local variables
    integer  :: io, ispin, j, num_u, num_v, jo, ind
    real(dp) :: E_HFX, time_start, time_end, time_dm_start, time_dm_end
    real(dp) :: spin_factor
    real(dp), dimension(:,:,:), pointer :: H_EXX => null(), DM_tmp => null()
    real(dp), dimension(:,:),   pointer :: P_max => null()
#ifdef MPI
    ! Global buffers for the storage of the sparse matrix
    ! FIXME: The arrays have the save attribute in the original version,
    !        check that pointers are OK
    integer :: BNode, nuog, maxnumh, MPIerror, maxnhg, maxndg, iio
    integer,  dimension(:),     pointer :: numhg => null(), &
&                                          listhptrg => null(), &
&                                          listhg => null()
    integer,  dimension(:),     pointer :: numdg => null(), &
&                                          listdptrg => null(), &
                                           listdg => null()
    real(dp), dimension(:,:),   pointer ::  Dscfg => null()
    real(dp), dimension(:,:,:), pointer ::  Hg_EXX => null()
#endif

    external :: timer

    ! -------------------------------------------------------------------------

    call timer('HFX', 1)

    ! FIXME: counting calls to this routine for debugging
    hfx_call_counter = hfx_call_counter + 1
    if ( Node .eq.  0 ) then
      write(6,'("setup_hfx: call #",i4.4)') hfx_call_counter
    end if

#ifdef MPI
    call re_alloc(numhg, 1, hfx_sysdata%nuotot, name='numhg', &
&     routine='setup_hfx')
    call re_alloc(listhptrg, 1, hfx_sysdata%nuotot, name='listhptrg', &
&     routine='setup_hfx')
    call re_alloc(numdg, 1, hfx_sysdata%nuotot, name='numdg', &
&     routine='setup_hfx')
    call re_alloc(listdptrg, 1, hfx_sysdata%nuotot, name='listdptrg', &
&     routine='setup_hfx')

    ! Globalise numh
    do io=1,hfx_sysdata%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node .eq. BNode ) then
        call GlobalToLocalOrb(io, Node, Nodes, iio)
        numhg(io) = hfx_sysdata%numh(iio)
        numdg(io) = hfx_sysdata%numh(iio)
      endif
      call MPI_Bcast(numhg(io), 1, MPI_integer, BNode, &
&       MPI_Comm_World, MPIerror)
      call MPI_Bcast(numdg(io), 1, MPI_integer, BNode, &
&       MPI_Comm_World, MPIerror)
    enddo

    ! Build global listhptr
    listhptrg(1) = 0
    listdptrg(1) = 0
    do io=2,hfx_sysdata%nuotot
      listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      listdptrg(io) = listdptrg(io-1) + numdg(io-1)
    enddo

    ! Globalise listh
    maxnhg = listhptrg(hfx_sysdata%nuotot) + numhg(hfx_sysdata%nuotot)
    maxndg = listdptrg(hfx_sysdata%nuotot) + numdg(hfx_sysdata%nuotot)
    call re_alloc(listhg, 1, maxnhg, name='listhg', routine='setup_hfx')
    call re_alloc(listdg, 1, maxndg, name='listdg', routine='setup_hfx')

    do io=1,hfx_sysdata%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node .eq. BNode ) then
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        do jo=1,numhg(io)
          listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = &
&           hfx_sysdata%listh(hfx_sysdata%listhptr(iio)+1: &
&             hfx_sysdata%listhptr(iio)+hfx_sysdata%numh(iio))
          listdg(listdptrg(io)+1:listdptrg(io)+numdg(io)) = hfx_sysdata%listh( &
&           hfx_sysdata%listhptr(iio)+1:hfx_sysdata%listhptr(iio)+hfx_sysdata%numh(iio))
        enddo
      endif

      call MPI_Bcast(listhg(listhptrg(io)+1), numhg(io), MPI_integer,&
&       BNode, MPI_Comm_World, MPIerror)
      call MPI_Bcast(listdg(listdptrg(io)+1), numdg(io), MPI_integer,&
&       BNode, MPI_Comm_World, MPIerror)
    enddo
    ! We trans the sparse matrix to full matrix. What's its nuo type?

    call re_alloc(Dscfg, 1, maxndg, 1, hfx_sysdata%nspin, name='Dscfg', &
&     routine='setup_hfx')
    Dscfg(:,:) = 0.0_dp

    ! Globalise Dscf
    do io=1,hfx_sysdata%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node .eq. BNode ) then
        call GlobalToLocalOrb(io, Node, Nodes, iio)
        do ispin=1,hfx_sysdata%nspin
          do jo=1,hfx_sysdata%numh(iio)
            Dscfg(listdptrg(io)+jo,ispin) = &
&             Dscf(hfx_sysdata%listhptr(iio)+jo,ispin)
          enddo
        enddo
      endif
      do ispin=1,hfx_sysdata%nspin
        call MPI_Bcast(Dscfg(listdptrg(io)+1,ispin), numdg(io), &
          MPI_double_precision, BNode, MPI_Comm_World, MPIerror)
      enddo
    enddo
#endif

    ! Build HFX
    call re_alloc(DM_tmp, 1, hfx_sysdata%norb, 1, hfx_sysdata%norb, 1, &
&     hfx_sysdata%nspin, name='DM_tmp', routine='setup_hfx')
    DM_tmp(:,:,:) = 0.0_dp
    call re_alloc(H_EXX, 1, hfx_sysdata%nuotot, 1, hfx_sysdata%norb, 1, &
&     hfx_sysdata%nspin, name='H_EXX', routine='setup_hfx')
    H_EXX(:,:,:) = 0.0_dp
    call re_alloc(P_max, 1, hfx_sysdata%norb, 1, hfx_sysdata%norb, &
&     name='P_max', routine='setup_hfx')
    P_max(:,:) = 0.0_dp

    call cpu_time(time_start)
    call cpu_time(time_dm_start)

#ifdef MPI
    ! Transfer sparse Dm to full matrix
    call sparse_dense(hfx_sysdata%nspin, hfx_sysdata%nuotot, &
&     hfx_sysdata%nuotot, hfx_sysdata%norb, maxndg, numdg, &
&     listdptrg, listdg, Dscfg, DM_tmp)
    call get_pmax_shell(hfx_sysdata%nspin, hfx_sysdata%norb, &
&     hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
&     hfx_sysdata%na, hfx_sysdata%isa, maxndg, numdg, listdptrg, &
&     listdg, DM_tmp, P_max)

    !> Calculate ERIs and store them in RAM
    !!
    !! \bug Was previously done conditionally through build_hfx_potential
    call timer('ERI',1)
	if(node .eq. 0) write(6,'(a, f12.6, a, I6)') "my mark"
    call evaluate_eri(libint_data, hfx_sysdata%nspin, hfx_sysdata%norb, &
&     hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
&     hfx_sysdata%na, hfx_sysdata%isa, hfx_sysdata%cell, hfx_sysdata%cell_r, &
&     hfx_optdata, DM_tmp, P_max, H_EXX)
    call timer('ERI',2)

    ! part H_EXX(norb,norb,nspin) per node : parallel over list_mn and list_uv
    ! u and m are local in this node, to calc H_EXX(u,m,nspin), otherwise it
    ! will be zero !
    call re_alloc(Hg_EXX, 1, hfx_sysdata%nuotot, 1, hfx_sysdata%norb, 1, &
&     hfx_sysdata%nspin, name='Hg_EXX', routine='setup_hfx')
    Hg_EXX(:,:,:) = 0.0_dp

    ! Get all H_EXX
    call MPI_AllReduce(H_EXX(1,1,1), Hg_EXX(1,1,1), &
&     hfx_sysdata%nuotot*hfx_sysdata%norb*hfx_sysdata%nspin, &
&     MPI_double_precision, MPI_Sum, MPI_Comm_World, MPIerror)
    H_EXX(:,:,:) = Hg_EXX(:,:,:)
    call de_alloc(Hg_EXX, name='Hg_EXX', routine='setup_hfx')
#else
    call sparse_dense(hfx_sysdata%nspin, hfx_sysdata%nuotot, &
&     hfx_sysdata%nuotot, hfx_sysdata%norb, hfx_sysdata%maxnh, &
&     hfx_sysdata%numh, hfx_sysdata%listhptr, hfx_sysdata%listh, Dscf, DM_tmp)
    call get_pmax_shell(hfx_sysdata%nspin, hfx_sysdata%norb, &
&     hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
&     hfx_sysdata%na, &
&     hfx_sysdata%isa, hfx_sysdata%maxnh, hfx_sysdata%numh, &
&     hfx_sysdata%listhptr, hfx_sysdata%listh, DM_tmp, P_max)

    !> Calculate ERIs and store them in RAM
    !!
    !! \bug Was previously done conditionally through build_hfx_potential
    call timer('ERI',1)
	write(6,'(a, f12.6, a, I6)') "my mark"
    call evaluate_eri(libint_data, hfx_sysdata%nspin, hfx_sysdata%norb, &
&     hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
&     hfx_sysdata%na, hfx_sysdata%isa, hfx_sysdata%cell, hfx_sysdata%cell_r, &
&     hfx_optdata, DM_tmp, P_max, H_EXX)
    call timer('ERI',2)
#endif

    call cpu_time(time_dm_end)

    write(6,'(a, f12.6, a, I6)') "ERIs time = ", time_dm_end-time_dm_start, &
                                 "(Sec), on Node =", Node

    if ( hfx_sysdata%nspin == 1 ) then
      spin_factor = 1.0_dp
    else
      spin_factor = 2.0_dp
    end if

    ! Here H_EXX and DM is global
    ! Calculate EXX from dense HFX and DM.
    E_HFX = 0.0_dp
    do ispin=1,hfx_sysdata%nspin
      do num_u=1,hfx_sysdata%nuotot
        do num_v=1,hfx_sysdata%norb
          E_HFX = E_HFX - 0.25_dp*spin_factor*H_EXX(num_u,num_v,ispin) &
&               * DM_tmp(num_u,num_v,ispin)*0.25_dp
        enddo
      enddo
    enddo

    Exc = Exc + E_HFX

    if ( Node .eq. 0 ) then
      write(6,'(a,f12.6,a)') 'setup_hfx: HFX energy:',E_HFX*13.60580_dp, "eV"
    endif


#ifdef MPI
    ! Global H_EXX to local Hmat, io: global orbital, iio: local in this
    ! node !
    do ispin=1,hfx_sysdata%nspin
      do iio=1,hfx_sysdata%nuo
        call LocalToGlobalOrb(iio, Node, Nodes, io)
        do j=1,hfx_sysdata%numh(iio)
          ind = hfx_sysdata%listhptr(iio) + j
          jo = hfx_sysdata%listh(ind)
          Hmat(ind,ispin) = Hmat(ind,ispin) &
&                         - H_EXX(io,jo,ispin)*spin_factor*0.5_dp*0.25_dp
        enddo
      enddo
    enddo
#else
    do ispin=1,hfx_sysdata%nspin
      do io=1,hfx_sysdata%nuo
        do j = 1,hfx_sysdata%numh(io)
          ind = hfx_sysdata%listhptr(io) + j
          jo = hfx_sysdata%listh(ind)
          !you need to a combine code  for ns=1, and ns=2
          Hmat(ind,ispin) = Hmat(ind,ispin) &
&                         - H_EXX(io,jo,ispin)*spin_factor*0.5_dp*0.25_dp
        enddo
      enddo
    enddo
#endif

    call cpu_time(time_end)

    if ( Node .eq. 0 ) then
      write(6,'(a, F12.6, a, I6)') "setup_hfx: Build HFX time =",  &
&       time_end-time_start, " (secs) on Node =", Node
    endif

    call de_alloc(H_EXX, name='H_EXX', routine='setup_hfx')
    call de_alloc(DM_tmp, name='DM_tmp', routine='setup_hfx')
    call de_alloc(P_max, name='P_max', routine='setup_hfx')
#ifdef MPI
    call de_alloc(Dscfg, name='Dscfg', routine='setup_hfx')
    call de_alloc(listdg, name='listdg', routine='setup_hfx')
    call de_alloc(listdptrg, name='listdptrg', routine='setup_hfx')
    call de_alloc(listhg, name='listhg', routine='setup_hfx')
    call de_alloc(listhptrg, name='listhptrg', routine='setup_hfx')
    call de_alloc(numdg, name='numdg', routine='setup_hfx')
    call de_alloc(numhg, name='numhg', routine='setup_hfx')
#endif

    call timer('HFX', 2)

  end subroutine setup_hfx

! *****************************************************************************
!> \brief Computes the Hartree-Fock potential corresponding to a given
!!        Hamiltonian
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] nspin: ...
!! \param[in] ncells: ...
!! \param[in] nuotot: ...
!! \param[in] norb: ...
!! \param[in] io: ...
!! \param[in] jo: ...
!! \param[in] ko: ...
!! \param[in] lo: ...
!! \param[in] nao_eri: ...
!! \param[in] DM_tmp: ...
!! \param[in,out] H_EXX: ...
! *****************************************************************************
  subroutine evaluate_eri(libint_data, nspin, norb, iaorb, iphorb, nuotot, na, isa, &
&                cell, rcell, hfx_optdata, DM_tmp, P_max, H_EXX)

    use parallel,       only: IONode, Node, Nodes
    use parallelsubs,   only: GetNodeOrbs, GlobalToLocalOrb, &
                              LocalToGlobalOrb,WhichNodeOrb
    use atm_types,      only: species, species_info, l_max, nco, nso, &
&                               maxn_contract, ncon_max, maxn_orbnl, &
&                               maxnorbs
    use atomlist,       only: indxuo
    use atmfuncs,       only: lofio, mofio
    use listsc_module,  only: listsc
    use nao2gto_index, only: indexsc
    use alloc,          only: re_alloc, de_alloc
    use nao2gto_types
    use nao2gto_data
    use nao2gto_contract,  only : calc_contract_eri
    use nao2gto_libint, only: Libint_t
#ifdef MPI
    use mpi_siesta
!    use nao2gto_parallel
#endif

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data
    integer,  intent(in) :: &
      nspin, na, norb, nuotot,  &
      iaorb(norb), iphorb(norb), isa(na)

    type(hfx_options_type), intent(in) :: hfx_optdata

    real(dp), intent(in) :: cell(3,3), rcell(3,3)
    real(dp), intent(in) :: DM_tmp(norb,norb,nspin)
    real(dp), intent(in) :: P_max(norb,norb)
    real(dp), intent(inout) :: H_EXX(nuotot, norb, nspin)

    ! Local variables
    integer  :: io, jo, ko, lo, is, js, ks, ls, ioa, joa, koa, loa,   &
&               l_i, l_j, l_k, l_l, m_i, m_j, m_k, m_l,               &
&               ncoi, ncoj, ncok, ncol, npgfi, npgfj, npgfk, npgfl,   &
&               nsoi, nsoj, nsok, nsol, num_a, num_b, num_c, num_d,   &
&               i_list_ij, i_list_kl, i_list_ij_local, list_ij_local, &
&               index_ij, index_kl, ncells
    integer :: ishell, jshell, kshell, lshell
    integer(int_8) :: shell_eri_calc, spher_eri_calc, spher_eri_store, &
&                     tot_pgto, neris_tmp
#ifdef MPI
    integer  :: MPIerror, Request, num_loc, Status(MPI_Status_Size)
#endif
    real(dp) :: eps_temp, DM_max, nao_eri
    real(dp) :: rij2, rkl2
    real(dp) :: max_contraction_val, max_val, max_val1, max_val2, &
&               max_val2_set, log10_pmax
    real(dp) :: ri(3), rj(3), rk(3), rl(3)
	integer :: list(2), list_pos(2), send_mes(3)
	integer :: block, count, buf, i, i_list_ijkl, temp_i_list_kl, tag,&
				i_interval, interval, request1, request2, request3, request4, other_msg, &
				task_msg, tar_msg, finish_tag, ask_buf, suc_buf, fall_buf,&
				min_task, end_count, start_pos, end_pos, pos_begin ,pos_end,&
				sou, tar, temp_node, ask_msg
	logical :: flag, flag2, termination
	logical, allocatable :: busy(:)
	real(dp) :: time_all, time_start, time_end

    real(dp), dimension(:,:,:,:), pointer :: eri => null()
    type(species_info), pointer :: ispp => null(), jspp => null(), &
&                                  kspp => null(), lspp => null()
    type(hfx_screen_coeff_type), dimension(:,:), pointer :: &
&     tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2

    ! -------------------------------------------------------------------------

    if ( IONode ) write(*,*) "[nao2gto] dm_trunc?", &
&     hfx_optdata%DM_trunc
    ncells = norb / nuotot
	tag = 0
	buf = 1

! #ifdef MPI
    ! call GetNodeOrbs(list_kl%nelement, Node, Nodes, list_kl_local)
! #endif
	termination = .false.
    neris_tmp = 0_int_8
    shell_eri_calc = 0_int_8
    spher_eri_calc = 0_int_8
    spher_eri_store = 0_int_8
    tot_pgto = 0_int_8
	
	request1 = 111111
	request2 = 222222
	request3 = 333333
	
	other_msg = 1234
	ask_msg = 4321234
	task_msg = 11111
	tar_msg =22222
	finish_tag =33333
	ask_buf = 44444
	suc_buf = 55555
	fall_buf = 66666
	
	min_task = 50
	i_interval = 1
  interval = 5000
  
  time_all = 0
	
	if( IONode ) then
		end_count = nodes-1
		allocate(busy(nodes - 1))
		
		do i =1, nodes-1
			busy(i) = .true.
		enddo

		do while(.true.)
			call MPI_RECV(buf, 1, MPI_integer, MPI_ANY_SOURCE, other_msg, &
					MPI_COMM_WORLD, status, MPIerror)
			sou = status(MPI_SOURCE)
			if(termination) then
				! call MPI_SEND(buf, 1, MPI_integer, sou, &
						! ask_msg, MPI_COMM_WORLD, MPIerror)
				call MPI_SEND(sou, 1, MPI_integer, sou, &
						tar_msg, MPI_COMM_WORLD, MPIerror)
				! print *,"send end meassage",sou
				end_count = end_count -1
				busy(sou) = .false.
			else
				do i =1, nodes-2 
					tar = sou + i 
					if(tar >= nodes) then
						tar = tar - nodes + 1
					endif
					if(busy(tar)) then
						CALL MPI_SEND(tar, 1, MPI_integer, sou, &
							tar_msg, MPI_COMM_WORLD, MPIerror)
						! print *, "tar: ", tar, sou
						call MPI_RECV(buf, 1, MPI_integer, sou, &
							finish_tag, MPI_COMM_WORLD, status, MPIerror)
						if(buf .eq. suc_buf) then
							exit
						endif
					endif
					if(i .eq. (nodes-2)) then
						! call MPI_SEND(sou, 1, MPI_integer, sou, &
							! ask_msg, MPI_COMM_WORLD, MPIerror)
						call MPI_SEND(sou, 1, MPI_integer, sou, &
							tar_msg, MPI_COMM_WORLD, MPIerror)
						! print *,"send end meassage",sou
						busy(sou) = .false.
            end_count = end_count - 1
            call cpu_time(time_start)
						termination = .true.
					endif
				enddo
			endif
		if(end_count .eq. 0) exit
		enddo
		
		deallocate(busy)
    ! call MPI_BARRIER(MPI_COMM_WORLD, MPIerror)
    call cpu_time(time_end)
    time_all = time_end - time_start
    print *,time_all,"terminal time"
	else
		call GetNodeOrbs(list_ij%nelement, Node-1, Nodes-1, end_pos)
		call LocalToGlobalOrb(1, Node-1, Nodes -1, i_list_ij)
		call LocalToGlobalOrb(end_pos, NOde-1, Nodes-1, i_list_kl)
		!print *,"node:", node, " task",i_list_ij, i_list_kl, end_pos
		count = end_pos
		start_pos = 1
		temp_node = Node-1
		
		termination = .false.
		call MPI_IRECV(buf, 1, MPI_integer, MPI_ANY_SOURCE, ask_msg, &
						 MPI_COMM_WORLD, request1, MPIerror)
		do while(.true.)
		i_list_ij_local = start_pos
		!print *,"start wock on node",node, start_pos, end_pos, temp_node
		do while(i_list_ij_local <= end_pos)
			count = end_pos - i_list_ij_local
			call LocalToGlobalOrb(i_list_ij_local, temp_node, Nodes-1, i_list_ij)
			!print *,"xxxxxx", i_list_ij
			! print *, i_list_ij, "xxxxxxxxxx", node
			io = list_ij%element(i_list_ij)%pair(1)
			jo = list_ij%element(i_list_ij)%pair(2)
			ri = list_ij%element(i_list_ij)%r1
			rj = list_ij%element(i_list_ij)%r2
			rij2 = list_ij%element(i_list_ij)%dist2

			is = isa(iaorb(io))
			ispp => species(is)
			ioa = iphorb(io)
			l_i = lofio(is, ioa)
			m_i = mofio(is, ioa)
			npgfi = ispp%orbnl_contract(ispp%orb_index(ioa))
			ncoi = nco(l_i)*npgfi

			js = isa(iaorb(jo))
			jspp => species(js)
			joa = iphorb(jo)
			l_j = lofio(js, joa)
			m_j = mofio(js, joa)
			npgfj = jspp%orbnl_contract(jspp%orb_index(joa))
			ncoj = nco(l_j)*npgfj

			ishell = ispp%orb_index(ioa)
			jshell = jspp%orb_index(joa)
			index_ij = list_ij%element(i_list_ij)%nl_index
			max_val1 = sfc_shell(jshell,ishell,js,is)%x(1)*rij2 + &
			&                sfc_shell(jshell,ishell,js,is)%x(2)

			do i_list_kl = 1, list_kl%nelement
				i_interval = i_interval + 1
				if(MOD(i_interval, interval) .eq. 0) then
					do while(.true.)
						CALL MPI_TEST(request1, flag, status, MPIerror)
						if(flag) THEN
							count = end_pos - i_list_ij_local
							tar = status(MPI_SOURCE)
							if(count > min_task) then
								send_mes(2) = (i_list_ij_local + end_pos)/2 +1
								send_mes(3) = end_pos
								send_mes(1) = send_mes(3) - send_mes(2) + 1
								end_pos = send_mes(2) - 1
								call MPI_SEND(send_mes(1), 3, MPI_integer, tar, &
									task_msg, MPI_COMM_WORLD, MPIerror)
								! print *,"send meassage finished", node,tar, send_mes(1), send_mes(2), send_mes(3)
								!print *,"node",node, i_list_ij_local , end_pos
							else
								send_mes(1) = 0
								call MPI_SEND(send_mes(1), 3, MPI_integer, tar, &
									task_msg, MPI_COMM_WORLD, MPIerror)
								! print *,"send meassage finished", node, status(MPI_SOURCE)
							endif
							call MPI_IRECV(buf, 1, MPI_integer, MPI_ANY_SOURCE, ask_msg, &
									MPI_COMM_WORLD, request1, MPIerror)
						else
							exit
						endif
					enddo
				endif
			
		ko = list_kl%element(i_list_kl)%pair(1)
        lo = list_kl%element(i_list_kl)%pair(2)
        rk = list_kl%element(i_list_kl)%r1
        rl = list_kl%element(i_list_kl)%r2
        rkl2 = list_kl%element(i_list_kl)%dist2

!        write(*,*) "[DEBUG] ko, lo, rk, rl, rkl2 = ", &
!          ko, lo, rk, rl, rkl2

        ks = isa(iaorb(ko))
        kspp => species(ks)
        koa = iphorb(ko)
        l_k = lofio(ks, koa)
        m_k = mofio(ks, koa)
        npgfk = kspp%orbnl_contract(kspp%orb_index(koa))
        ncok =nco(l_k)*npgfk

        ls = isa(iaorb(lo))
        lspp => species(ls)
        loa = iphorb(lo)
        l_l = lofio(ls, loa)
        m_l = mofio(ls, loa)
        npgfl = lspp%orbnl_contract(lspp%orb_index(loa))
        ncol = nco(l_l)*npgfl

        kshell = kspp%orb_index(koa)
        lshell = lspp%orb_index(loa)
        index_kl = list_kl%element(i_list_kl)%nl_index

        if ( index_kl .le. index_ij ) then

          if ( um_cut(io,ko) .and. um_cut(io,lo) .and. um_cut(jo,ko) .and. &
&              um_cut(jo,lo) ) cycle

          max_val2_set = sfc_shell(lshell,kshell,ls,ks)%x(1)*rkl2 + &
&           sfc_shell(lshell,kshell,ls,ks)%x(2)
          max_val = max_val1 + max_val2_set

          eps_temp = eri_prescreen(io, jo)*eri_prescreen(ko, lo)
          eps_temp = dsqrt(eps_temp)

          if ( hfx_optdata%DM_trunc ) then
            DM_max = max( P_max(io, ko), P_max(io, lo), &
                     P_max(jo, ko), P_max(jo, lo) )
            if ( DM_max <= 0.0_dp ) then
              log10_pmax = log_zero
            else
              log10_pmax = log10(DM_max)
            end if
            eps_temp = DM_max*eps_temp
          endif

          if ( eps_temp .gt.hfx_optdata%eps_schwarz) then

            shell_eri_calc = shell_eri_calc + 1

            tmp_R_1 => pair_dist_radii_pgf(:,:,jshell,ishell,js,is)
            tmp_R_2 => pair_dist_radii_pgf(:,:,lshell,kshell,ls,ks)
            tmp_screen_pgf1 => sfc_pgf(:,:,jshell,ishell,js,is)
            tmp_screen_pgf2 => sfc_pgf(:,:,lshell,kshell,ls,ks)

             max_contraction_val = ispp%orbnl_contraction_coeff(ishell) * &
                                   jspp%orbnl_contraction_coeff(jshell) * &
                                   kspp%orbnl_contraction_coeff(kshell) * &
                                   lspp%orbnl_contraction_coeff(lshell) * &
                                   DM_max

            call re_alloc(eri, 1, nso(l_i), 1, nso(l_j), 1, nso(l_k), &
&             1, nso(l_l), name='eri', routine='evaluate_eri')
            eri(:,:,:,:) = 0.0_dp

            call calc_contract_eri( libint_data, cell, rcell, ri, rj, rk, rl, &
&             npgfi, npgfj, npgfk, npgfl, l_i, l_j, l_k, l_l, &
&             ncoi, ncoj, ncok, ncol, &
&             ispp%orbnl_zeta(1:npgfi,ispp%orb_index(ioa)), &
&             jspp%orbnl_zeta(1:npgfj,jspp%orb_index(joa)), &
&             kspp%orbnl_zeta(1:npgfk,kspp%orb_index(koa)), &
&             lspp%orbnl_zeta(1:npgfl,lspp%orb_index(loa)), &
&             ispp%sphi(1:ncoi,ioa:ioa+nso(l_i)-1), &
&             jspp%sphi(1:ncoj,joa:joa+nso(l_j)-1), &
&             kspp%sphi(1:ncok,koa:koa+nso(l_k)-1), &
&             lspp%sphi(1:ncol,loa:loa+nso(l_l)-1), &
&             neris_tmp, max_contraction_val, max_val2_set, log10_pmax, &
&             tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2, &
&             hfx_optdata, eri)

            tot_pgto = tot_pgto + neris_tmp

            do nsoi = 1, nso(l_i)
              do nsoj = 1, nso(l_j)
                do nsok = 1, nso(l_k)
                  do nsol = 1, nso(l_l)
                    spher_eri_calc = spher_eri_calc + 1
                    if ( DM_max*dabs(eri(nsoi,nsoj,nsok,nsol)*2) .gt. &
&                        hfx_optdata%eps_stored ) then
                      spher_eri_store = spher_eri_store + 1
                      num_a = io + nsoi-1
                      num_b = jo + nsoj-1
                      num_c = ko + nsok-1
                      num_d = lo + nsol-1
                      nao_eri = eri(nsoi,nsoj,nsok,nsol)*2

                      call hfx_matrix(nspin, ncells, nuotot, norb, &
&                       num_a, num_b, num_c, num_d, nao_eri, DM_tmp, H_EXX )
                    endif
                  enddo
                enddo
              enddo
           enddo

           call de_alloc(eri, name='eri', routine='evaluate_eri')

          endif ! Schwarz inequality
        endif
		enddo
		i_list_ij_local = i_list_ij_local + 1
		enddo
			
		call MPI_ISEND(ask_buf, 1, MPI_integer, 0, other_msg, MPI_COMM_WORLD, request4 ,MPIerror)
		! print *,"node",node,"ask"
		do while(count .eq. 0)
        
		call MPI_IRECV(tar, 1, MPI_integer, 0, tar_msg, MPI_COMM_WORLD, request2, MPIerror)
			do while(.true.)
				!print *, "cycle on node:", node
				call MPI_TEST(request2, flag2, status, MPIerror)
				
				if(.not.flag2) then
					do while(.true.)
						CALL MPI_TEST(request1, flag, status, MPIerror)
						if(flag) THEN
							send_mes(1) = 0
							call MPI_SEND(send_mes(1), 3, MPI_integer, status(MPI_SOURCE), &
								task_msg, MPI_COMM_WORLD, MPIerror)
							call MPI_IRECV(buf, 1, MPI_integer, MPI_ANY_SOURCE, ask_msg, &
										MPI_COMM_WORLD, request1, MPIerror)
							! print *,"send meassage finished", node, status(MPI_SOURCE)
						else
							 exit
						endif
					enddo
				else
					exit
				endif
			enddo
		! print *,node, "recive tar", tar
		if(tar .eq. node) then
			termination = .true.
			exit
		else if(tar .gt. 0) then
			call MPI_ISEND(ask_buf, 1, MPI_integer, tar, ask_msg, MPI_COMM_WORLD, request3, MPIerror)
			call MPI_WAIT(request3, status, MPIerror)
			call MPI_RECV(send_mes(1), 3, MPI_integer, tar, task_msg, MPI_COMM_WORLD, status, MPIerror)
			count = send_mes(1)
			start_pos = send_mes(2)
			end_pos = send_mes(3)
			temp_node = tar - 1
			if( count > 0 ) then
				call MPI_SEND(suc_buf, 1, MPI_integer, 0, finish_tag, MPI_COMM_WORLD, MPIerror)
			else
				call MPI_SEND(fall_buf, 1, MPI_integer, 0, finish_tag, MPI_COMM_WORLD, MPIerror)
			endif
		endif
		enddo

		if(termination) exit
	enddo 
		call MPI_ISEND(ask_buf, 1, MPI_integer, node, ask_msg, MPI_COMM_WORLD, request3, MPIerror)
		! print *,"node :",node, "quit"
	! call MPI_BARRIER(MPI_COMM_WORLD, MPIerror)
	endif
						
  end subroutine evaluate_eri

! *****************************************************************************
!> \brief Computes the Hartree-Fock potential corresponding to a given
!!        Hamiltonian
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] nspin: ...
!! \param[in] ncells: ...
!! \param[in] nuotot: ...
!! \param[in] norb: ...
!! \param[in] io: ...
!! \param[in] jo: ...
!! \param[in] ko: ...
!! \param[in] lo: ...
!! \param[in] nao_eri: ...
!! \param[in] DM_tmp: ...
!! \param[in,out] H_EXX: ...
! *****************************************************************************
  subroutine hfx_matrix(nspin, ncells, nuotot, norb, io, jo, ko, lo, &
                 nao_eri, DM_tmp, H_EXX)

    use atomlist,      only: indxuo
    use nao2gto_index, only: indexsc
    use nao2gto_data,  only: subshell

    implicit none

    ! Arguments
    integer,  intent(in)    :: nspin, ncells, nuotot, norb, io, jo, ko, lo
    real(dp), intent(in)    :: nao_eri
    real(dp), intent(in)    :: DM_tmp(norb,norb,nspin)
    real(dp), intent(inout) :: H_EXX(nuotot,norb,nspin)

    ! Local variables
    integer  :: ispin, iuo, juo, kuo, luo, llo, jshell, lshell, iushell, &
&               jushell, kushell, lushell, index_ij, index_kl, &
&               io_trans, jo_trans, ko_trans, lo_trans
    real(dp) :: gint

    ! -------------------------------------------------------------------------

    gint = nao_eri
    iuo = io ! u is always u0
    juo = indxuo(jo)
    kuo = indxuo(ko)
    luo = indxuo(lo)
    llo = indexsc(ko, kuo, lo)

    ! num_n have to trans to play with m0 to get num_mn and
    ! campared to num_uv, so there is num_n_1
    iushell = subshell(iuo)
    jushell = subshell(juo)
    kushell = subshell(kuo)
    lushell = subshell(luo)
    jshell  = subshell(jo)
    lshell  = subshell(llo)

    index_ij = ncells*iushell*(iushell-1)/2 + &
&     ((jshell-1)/subshell(nuotot))*iushell + jushell

    index_kl = ncells*kushell*(kushell-1)/2 + &
&     ((lshell-1)/subshell(nuotot))*kushell + lushell

    if ( iushell  .eq. jushell  ) gint = gint*0.5d0
    if ( kushell  .eq. lushell  ) gint = gint*0.5d0
    if ( index_ij .eq. index_kl ) gint = gint*0.5d0

    !! HFX
    do ispin=1,nspin

      !         (u0vR|mR'nR")        =       (u0vR|nR"mR')
      ! = (v0u[-R]|m[R'-R]n[R"-R])   = (v0u[-R]|n[R"-R]m[R'-R])
      ! = (m0n[R"-R']|u[-R']v[R-R']) = (m0n[R"-R']|v[R-R']u[-R'])
      ! = (n0m[R'-R"]|u[-R"]v[R-R"]) = (n0m[R'-R"]|v[R-R"]u[-R"])

      ! 1.VEE(1[0]  2[H] | 3[G]  4[N])  (u0v[R]|m[R']n[R"])
      !   VEE(1[0]  2[H] | 4[N]  3[G])  (u0v[R])|n[R"]m[R'])
      H_EXX(io,ko,ispin) = H_EXX(io,ko,ispin) &
                         + gint*DM_tmp(jo,lo,ispin)
      H_EXX(io,lo,ispin) = H_EXX(io,lo,ispin) &
                         + gint*DM_tmp(jo,ko,ispin)

      ! 2.VEE(2[0] 1[-H]| 3[G-H] 4[N-H])  (v0u[-R]|m[R'-R]n[R"-R])
      !   VEE(2[0] 1[-H]| 4[N-H] 3[G-H])  (v0u[-R]|n[R"-R]m[R'-R])
      io_trans = indexsc(jo, juo, io)
      ko_trans = indexsc(jo, juo, ko)
      lo_trans = indexsc(jo, juo, lo)

      H_EXX(juo,ko_trans,ispin) = H_EXX(juo,ko_trans,ispin) &
&                               + gint*DM_tmp(io_trans,lo_trans,ispin)
      H_EXX(juo,lo_trans,ispin) = H_EXX(juo,lo_trans,ispin) &
&                               + gint*DM_tmp(io_trans,ko_trans,ispin)

      ! 3.VEE(3[0]  4[N-G] | 1[-G] 2[H-G]) (m0n[R"-R']|u[-R']v[R-R'])
      !   VEE(3[0]  4[N-G] |2[H-G] 1[-G] ) (m0n[R"-R']|v[R-R']u[-R'])
      io_trans = indexsc(ko, kuo, io)
      jo_trans = indexsc(ko, kuo, jo)
      lo_trans = indexsc(ko, kuo, lo)

      H_EXX(kuo,io_trans,ispin) = H_EXX(kuo,io_trans,ispin) &
&                               + gint*DM_tmp(lo_trans,jo_trans,ispin)
      H_EXX(kuo,jo_trans,ispin) = H_EXX(kuo,jo_trans,ispin) &
&                               + gint*DM_tmp(lo_trans,io_trans,ispin)

      ! 4.VEE(4[0]  3[G-N] | 1[-N] 2[H-N])  (n0m[R'-R"]|u[-R"]v[R-R"])
      !   VEE(4[0]  3[G-N] | 2[H-N] 1[-N])  (n0m[R'-R"]|v[R-R"]u[-R"])
      io_trans = indexsc(lo, luo, io)
      jo_trans = indexsc(lo, luo, jo)
      ko_trans = indexsc(lo, luo, ko)

      H_EXX(luo,io_trans,ispin) = H_EXX(luo,io_trans,ispin) &
&                               + gint*DM_tmp(ko_trans,jo_trans,ispin)
      H_EXX(luo,jo_trans,ispin) = H_EXX(luo,jo_trans,ispin) &
&                               + gint*DM_tmp(ko_trans,io_trans,ispin)

    enddo

  end subroutine hfx_matrix

end module nao2gto_hfx
