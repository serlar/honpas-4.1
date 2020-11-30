! This file is part of the ONPAS package.
!
!
! Module coded by xmqin, Oct.12, 2013
! Edited by  xmqin, Nov. 27, 2016
!
module nao2gto_hfx_force

  use precision,       only : dp
  use atm_types,       only : l_max, nspecies, &
                              maxn_contract, maxn_orbnl
  use atm_types,       only : species, species_info, nco, nso
  use nao2gto_types
  use alloc,           only :  re_alloc, de_alloc
  use parallel,        only :  Node, Nodes
  use atmfuncs,        only : lofio, mofio
  use listsc_module,   only :  listsc
  use nao2gto_index
  use nao2gto_prescreen
  use atomlist,        only :  indxuo
  use nao2gto_utils
  use nao2gto_contract,    only :  calc_contract_eri
  use nao2gto_libint

  use parallel,       only : Node, Nodes
  use parallelsubs,   only : GetNodeOrbs, GlobalToLocalOrb, &
                             LocalToGlobalOrb,WhichNodeOrb
#ifdef MPI
  use mpi_siesta
#endif
  use atomlist,       only : indxuo
  use nao2gto_contract,  only : calc_contract_deriv_eri

  implicit none

  private

  public :: setup_hfx_force

  real(dp), save :: max_eri_force = 0.0_dp

contains

  subroutine setup_hfx_force(libint_data, hfx_optdata, hfx_sysdata, fal)

    use sparse_matrices, only: Dscf
    use nao2gto_types

    use precision
    use parallel,         only : Node, Nodes
    use parallelsubs,     only : GetNodeOrbs, GlobalToLocalOrb, &
                                 LocalToGlobalOrb, WhichNodeOrb
#ifdef MPI
    use mpi_siesta
#endif
    use alloc,            only : re_alloc, de_alloc
    use nao2gto_dm
  
    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data
    type(hfx_options_type), intent(in) :: hfx_optdata
    type(hfx_system_type), intent(in) :: hfx_sysdata
    real(dp), intent(inout) :: fal(3,hfx_sysdata%nua)

    ! Arguments
    integer, pointer :: &
      maxnd, maxnh, nua, na, norb, nspin, nuo, nuotot,          &
      iaorb(:), indxua(:), iphorb(:), isa(:),           &
      listd(:), listh(:), listdptr(:), listhptr(:),         &
      numd(:), numh(:),  nsc(:)
    real(dp), pointer :: xa(:,:)
    real(dp) :: cell(3,3)
    integer :: i, ia, io, iio, is, ispin, j, n

#ifdef MPI
    ! Global buffers for the storage of sparse matrix
    integer :: BNode, nuog, maxnumh, MPIerror
    integer, save :: maxnhg, maxndg
    integer, dimension(:), allocatable, save      ::   &
      numhg, listhptrg, listhg
    integer, dimension(:), allocatable, save      ::   &
      numdg, listdptrg, listdg
    real(dp), dimension(:,:), allocatable, save   ::  Dscfg
#endif

    integer num_u,num_v,num_m,num_n
    integer iu,ju,jo,ind

    real(dp), dimension(:,:,:), allocatable, save  ::  DM_tmp
    real(dp), dimension(:,:), allocatable, save  :: P_max

    external memory
! ----------------------------------------------------------------------

    ! FIXME: to be removed
    maxnd => hfx_sysdata%maxnh
    maxnh => hfx_sysdata%maxnh
    nua => hfx_sysdata%nua
    na => hfx_sysdata%na
    norb => hfx_sysdata%norb
    nspin => hfx_sysdata%nspin
    nuo => hfx_sysdata%nuo 
    nuotot => hfx_sysdata%nuotot
    iaorb => hfx_sysdata%iaorb
    indxua => hfx_sysdata%indxua
    iphorb => hfx_sysdata%iphorb
    isa => hfx_sysdata%isa
    listd => hfx_sysdata%listh
    listh => hfx_sysdata%listh
    listdptr => hfx_sysdata%listhptr
    listhptr => hfx_sysdata%listhptr
    numd => hfx_sysdata%numh
    numh => hfx_sysdata%numh
    nsc => hfx_sysdata%nsc
    cell(:,:) = hfx_sysdata%cell(:,:)
    xa => hfx_sysdata%xa

#ifdef MPI
     
      allocate(numhg(nuotot))
      call memory('A','I',nuotot,'update_hfx')
      allocate(listhptrg(nuotot))
      call memory('A','I',nuotot,'update_hfx')
      allocate(numdg(nuotot))
      call memory('A','I',nuotot,'update_hfx')
      allocate(listdptrg(nuotot))
      call memory('A','I',nuotot,'update_hfx')
     
!      call cpu_time(time_start) 
! Globalise numh
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
          numdg(io) = numd(iio)
        endif
        call MPI_Bcast(numhg(io),1,MPI_integer,BNode, &
                       MPI_Comm_World,MPIerror)
        call MPI_Bcast(numdg(io),1,MPI_integer,BNode, &
                       MPI_Comm_World,MPIerror)
      enddo

! Build global listhptr
      listhptrg(1) = 0
      listdptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
        listdptrg(io) = listdptrg(io-1) + numdg(io-1)
      enddo

! Globalse listh
      maxnhg = listhptrg(nuotot) + numhg(nuotot)
      maxndg = listdptrg(nuotot) + numdg(nuotot)
      allocate(listhg(maxnhg))
      allocate(listdg(maxndg))
      call memory('A','I',maxnhg,'update_hfx')
      call memory('A','I',maxndg,'update_hfx')

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numhg(io)
            listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = &
                 listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
            listdg(listdptrg(io)+1:listdptrg(io)+numdg(io)) = &
                 listd(listdptr(iio)+1:listdptr(iio)+numd(iio))
          enddo
        endif

        call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer,&
                       BNode,MPI_Comm_World,MPIerror)
        call MPI_Bcast(listdg(listdptrg(io)+1),numdg(io),MPI_integer,&
                       BNode,MPI_Comm_World,MPIerror)
      enddo
! We trans the sparse matrix to full matrix. What's its nuo type?

      allocate(Dscfg(maxndg,nspin))
      call memory('A','D',maxndg*nspin,'update_hfx')
      Dscfg(1:maxndg,1:nspin)=0.0d0

!   Globalise Dscf
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do ispin = 1,nspin
            do jo = 1,numd(iio)
           Dscfg(listdptrg(io)+jo,ispin) = Dscf(listdptr(iio)+jo,ispin)
            enddo
          enddo
        endif
        do ispin = 1,nspin
          call MPI_Bcast(Dscfg(listdptrg(io)+1,ispin),numdg(io), &
                         MPI_double_precision,BNode, &
                         MPI_Comm_World,MPIerror)
        enddo
      enddo
!      call cpu_time(time_end1)
!      if (node.eq.0) then
!      write(6,'(a,f16.6,a)') "global index",time_end1-time_start, "secs"
!      endif

#endif

!--------------------end for global Dscfg and Hmatg--------------------

!      if (allocated(DM_tmp)) then
!         deallocate(DM_tmp)
!      endif

      allocate(DM_tmp(norb,norb,nspin))
      call memory('A','D',norb*norb*nspin,'update_hfx')
      DM_tmp(1:norb,1:norb,1:nspin)=0.0d0

      allocate(P_max(norb,norb))
      call memory('A','D',norb*norb,'update_hfx')
      P_max(1:norb,1:norb)=0.0d0

#ifdef MPI
!---------------tranf sparse Dm to full matrix--------------------------
      call sparse_dense( nspin, nuotot, nuotot, norb, maxndg, numdg, &
                         listdptrg, listdg, Dscfg, DM_tmp )
      call get_pmax_shell( nspin, norb, iaorb, iphorb, nuo, na, &
                           isa, maxndg, numdg, listdptrg, listdg, &
                           DM_tmp, P_max )


!      call cpu_time(time_end2)
!      if (node.eq.0) then
!      write(6,'(a,f16.6,a)')'sparse to dense',time_end2-time_end1,'secs'
!      endif

      call build_hfx_gradient(hfx_optdata, hfx_sysdata, maxnhg, numhg, &
&       listhptrg, listh, DM_tmp, P_max, Fal)
#else

      call sparse_dense( nspin, nuotot, nuotot, norb, maxnd, numd, &
                         listdptr, listd, Dscf, DM_tmp )

      call get_pmax_shell( nspin, norb, iaorb, iphorb, nuo, na, &
                           isa, maxnd, numd, listdptr, listd, &
                           DM_tmp, P_max )
      
      call build_hfx_gradient(hfx_optdata, hfx_sysdata, hfx_sysdata%maxnh, &
&       hfx_sysdata%numh, hfx_sysdata%listhptr, hfx_sysdata%listh, &
&       DM_tmp, P_max, Fal)
#endif

      deallocate(DM_tmp)
      deallocate(P_max)

#ifdef MPI
      deallocate(listhg)
      deallocate(listdg)
      deallocate(listhptrg)
      deallocate(listdptrg)
      deallocate(numhg)
      deallocate(numdg)
      deallocate(Dscfg)
#endif

      end subroutine setup_hfx_force

! ***** HFX_GRADIENT MODULE BEGIN *****

! This file is part of the HONPAS package.
!

! Module coded by xmqin, October 2013
!
! calculate HFX(u,v) 
!
! HFX(u,v) = (um|vn) * Dm (m,n)
!

!------------------------------------------------------------------------------------

      subroutine build_hfx_gradient(hfx_opts, hfx_sys, maxnh, numh, &
&       listhptr, listh, DM_tmp, P_max, Fal)

      use nao2gto_data
      use nao2gto_prescreen

      implicit none

!---------------------------- INPUT  VARIABLES -------------------------------------

      type(hfx_options_type), intent(in) :: hfx_opts
      type(hfx_system_type), intent(in) :: hfx_sys
      integer, intent(in) :: maxnh
      integer, intent(in) :: numh(hfx_sys%nuotot)
      integer, intent(in) :: listhptr(hfx_sys%nuotot)
      integer, intent(in) :: listh(maxnh)

      real(dp),  intent(in)   :: DM_tmp(hfx_sys%norb,hfx_sys%norb,hfx_sys%nspin)
      real(dp),  intent(in)   :: P_max(hfx_sys%norb, hfx_sys%norb)

!--------------------------- INOUT VARIABLES ---------------------------------------
      real(dp), intent(inout) :: Fal(3,hfx_sys%nua)

!------------------------- TEMPOS, INTERNAL VARIABLES ------------------------------
      integer, pointer     :: &
        nua, na, norb, nspin, nuo, nuotot, &
        iaorb(:), indxua(:), iphorb(:), isa(:),  &
        nsc(:)
      real(dp),  pointer   :: xa(:,:)
      real(dp) :: cell(3,3)
      integer ncells, i, ind, j, io, iu, ju, iuo, juo, jo, jjo, ioa, &
              is, l, m, nshells, ispin, ia,js, ishell,jshell,ipgf, jpgf,&
              l_j,m_j,joa

      real(dp) tmax, gint, time_start, time_end

      real(dp)                  ::  scell(3,3), rscell(3,3)
      type(Libderiv_t)           ::  deriv
!      logical, save             ::  frstme = .true.
!      type(hfx_screen_coeff_type), &
!      dimension(:, :, :, :, :, :), pointer, save  :: &
!                                   pair_dist_radii_pgf, &
!                                   sfc_pgf
!      type(hfx_screen_coeff_type), &
!      dimension(:, :, :, :), pointer, save  :: &
!                                   sfc_shell
!      type(hfx_screen_coeff_type), &
!      dimension(:, :), pointer, save  :: &
!                                   sfc_kind



      external memory, timer

    ! FIXME: to be removed
    nua => hfx_sys%nua
    na => hfx_sys%na
    norb => hfx_sys%norb
    nspin => hfx_sys%nspin
    nuo => hfx_sys%nuo 
    nuotot => hfx_sys%nuotot
    iaorb => hfx_sys%iaorb
    indxua => hfx_sys%indxua
    iphorb => hfx_sys%iphorb
    isa => hfx_sys%isa
    nsc => hfx_sys%nsc
    xa => hfx_sys%xa
    cell(:,:) = hfx_sys%cell(:,:)

!
! Find supercell and ncells

      ncells=nsc(1)*nsc(2)*nsc(3)

      do  i = 1,3
        do j = 1,3
          scell(j,i) = cell(j,i) * nsc(i)
        enddo
      enddo

      call reclat(scell, rscell, 0)

! If this routine is called first time, read and print hfx input information
!       call read_hfx_info(hfx_opts)
!       call print_hfx_info(hfx_opts)

! supercell orbital initialisation
        call nao2gto_index_init(nsc, nuotot)

! gamma function table and Libint initialisation 

!        call init_md_ftable(4*l_max)
!        call initialize_libint(lib, l_max)
        call initialize_libderiv(deriv,l_max)        

!  orbital shell initialisation
        nullify(subshell)
        call re_alloc(subshell,1,norb)
        subshell(1:norb) = 0
! 2l+1 orbitals form a shell
        i = 0
        do io = 1,norb
           is = isa(iaorb(io))
           ioa = iphorb(io)
           l = lofio(is,ioa)
           m = mofio(is,ioa)
           if(m.eq.-l) then
              i = i + 1
              subshell(io) = i
           else
              subshell(io) = i
           endif
        enddo
        nshells = subshell(norb)
!        if(node.eq.0) &
!           write(6,'(A, I8, A, I8)') &
!           'Supercell orbitals = ', norb, ' shells = ', nshells
!      endif

        
!      if(frstme .or. (.not.samexa)) then
        nullify(eri_prescreen)
        call re_alloc( eri_prescreen, 1, norb, 1, norb,  &
               name='eri_prescreen', routine='build_hfx_gradient' )  
        eri_prescreen(1:norb,1:norb)=0.0d0

!!  To calculate ERI prescreen matrix

        call calc_prescreen_eri( hfx_libint, hfx_options, hfx_system, &
&         maxnh, numh, listhptr, listh, max_eri_force)


        nullify(um_cut)
        allocate(um_cut(norb,norb))
        um_cut(1:norb,1:norb)=.true.

        do io = 1,norb
           is = isa(iaorb(io))
           ioa = iphorb(io)
           l = lofio(is,ioa)
           m = mofio(is,ioa)
           if(m.ne.-l) cycle

           iu=indxuo(io)
           do j=1,numh(iu)
              ind=listhptr(iu) +j
              ju=listh(ind)
              jo=listsc(io,iu,ju)
              js=isa(iaorb(jo))
              joa=iphorb(jo)
              l_j=lofio(js,joa)
              m_j=mofio(js,joa)
              if(m_j.ne.-l_j) cycle
              um_cut(io,jo)=.false.
           ! write(8,*) num_v,num_n
           enddo
       enddo

        nullify(pair_dist_radii_pgf)
        nullify(sfc_pgf)
        nullify(sfc_shell)
        nullify(sfc_kind)

        allocate(pair_dist_radii_pgf(maxn_contract, maxn_contract, &
                       maxn_orbnl, maxn_orbnl, nspecies, nspecies))
        allocate(sfc_pgf(maxn_contract, maxn_contract, &
                       maxn_orbnl, maxn_orbnl, nspecies, nspecies))
        allocate(sfc_shell(maxn_orbnl, maxn_orbnl, &
                                                 nspecies, nspecies))
        allocate(sfc_kind(nspecies, nspecies))

      DO is = 1,nspecies
        DO js = 1,nspecies
           sfc_kind(js,is)%x(:) = 0.0_dp
          DO ishell=1,maxn_orbnl
            DO jshell=1,maxn_orbnl
               sfc_shell(jshell,ishell,js,is)%x(:) = 0.0_dp
              DO ipgf=1,maxn_contract
                DO jpgf=1,maxn_contract
                   pair_dist_radii_pgf(jpgf,ipgf,jshell,&
                                       ishell,js,is)%x(:) = 0.0_dp
                   sfc_pgf(jpgf,ipgf,jshell, &
                                            ishell,js,is)%x(:) = 0.0_dp
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

      call calc_pair_dist_radii(hfx_options, pair_dist_radii_pgf)

      call calc_screening_functions(hfx_libint, hfx_options, scell, rscell, &
                                    sfc_pgf, sfc_shell,&
                                    pair_dist_radii_pgf)

      coeffs_kind_max0=MAXVAL(sfc_shell(:,:,:,:)%x(2))
!      write(6,*) "max0 = ", coeffs_kind_max0

!! build left and right orbital pair_list
       if (associated(list_ij%element)) then
         deallocate(list_ij%element)
         list_ij%element => null()
       endif
       if (associated(list_kl%element)) then
         deallocate(list_kl%element)
         list_kl%element => null()
       endif  
         
       allocate(list_ij%element(ncells*nuotot*(nuotot+1)/2))
       allocate(list_kl%element(ncells**2*nuotot*(nuotot+1)/2))

       
        call build_pair_list(hfx_sys%nua, hfx_sys%na, hfx_opts, hfx_sys, &
                              eri_prescreen, max_eri_force, list_ij)

        call build_pair_list(hfx_sys%na, hfx_sys%na, hfx_opts, hfx_sys, &
                              eri_prescreen, max_eri_force, list_kl)

        log10_eps_schwarz = LOG10(hfx_options%eps_schwarz)



        if(node.eq.0) then
          write(6,*) ncells*nuotot*(nuotot+1)/2, ncells**2*nuotot*(nuotot+1)/2
          write(6,'(a,2x,I12)')  'list_ij:',list_ij%nelement
          write(6,'(a,2x,I12)')  'list_kl:',list_kl%nelement
          write(6,'(a,f12.9)')   'max_eri_force', max_eri_force
        endif

!      endif
      
      call timer('HFX_gradient',1)    
!      call cpu_time(time_start)
      call  evaluate_gradient( deriv, nspin, norb, iaorb, iphorb, nuotot, nua, na, isa,    &
                               scell, rscell, eri_prescreen, hfx_options, DM_tmp, P_max, &
                               Fal )


!      call cpu_time(time_end)
!      if(node.eq.0) then
!      write(6,'(a, f12.6, a)') "dERIs time = ", node, time_end-time_start, " s"
!      endif
      call timer('HFX_gradient',2)

      end subroutine build_hfx_gradient

! ***** HFX_GRADIENT MODULE END *****

! ***** GRADIENT MODULE BEGIN *****

! Writen by xmqin, October 2013

      subroutine evaluate_gradient( deriv, nspin, norb, iaorb, iphorb, nuotot,  &
                              nua, na, isa, cell, rcell, eri_prescreen_arg, &
                              hfx_opts, DM_tmp, P_max,  &
!                              radii_pgf, screen_coeffs_pgf, screen_coeffs_shell,  &
                              Fal )

      use nao2gto_common
      use nao2gto_data

      implicit none

! ----------------------------INPUT--------------------------------
      !type(Libint_t)        :: lib
      type(Libderiv_t)      :: deriv
      integer,  intent(in) :: &
        nspin, nua, na, norb, nuotot,  &
        iaorb(norb), iphorb(norb), isa(na)

      type(hfx_options_type) :: hfx_opts

      real(dp), intent(in) :: eri_prescreen_arg(norb,norb), DM_tmp(norb,norb,nspin)
      real(dp), intent(in) :: P_max(norb,norb)

      real(dp), intent(in) :: cell(3,3), rcell(3,3)
      real(dp), intent(inout) :: Fal(3, nua)


      TYPE(hfx_screen_coeff_type), &
        DIMENSION(:,:,:,:,:,:), POINTER        :: screen_coeffs_pgf, radii_pgf
      TYPE(hfx_screen_coeff_type), &
        DIMENSION(:,:,:,:), POINTER            :: screen_coeffs_shell


!------------------------------------------------------------------
!------------------------------------------------------------------
      type(species_info), pointer :: ispp, jspp, kspp, lspp
      real(dp)  ri(3), rj(3), rk(3), rl(3), rij2,  rkl2

      real(dp), dimension(:,:,:,:,:), allocatable :: eri_deriv
      real(dp) eps_temp, DM_max,symm_factor

      real(dp) max_contraction_val, max_val, max_val1, &
               max_val2, pmax_entry

      TYPE(hfx_screen_coeff_type), &
        DIMENSION(:,:), POINTER                :: tmp_R_1, tmp_R_2, &
                                     tmp_screen_pgf1, tmp_screen_pgf2
      REAL(dp)       :: max_val2_set, log10_pmax

      
      real(dp) :: nao_eri_deriv(12)
      integer io, jo, ko, lo, is, js, ks, ls, ioa, joa, koa, loa,   &
              l_i, l_j, l_k, l_l, m_i, m_j, m_k, m_l,               &
              ncoi, ncoj, ncok, ncol, npgfi, npgfj, npgfk, npgfl,   &
              nsoi, nsoj, nsok, nsol, num_a, num_b, num_c, num_d,   &
              i_list_ij, i_list_kl, i_list_kl_local, list_kl_local, &
              index_ij, index_kl, ia, ja, ka, la , ishell, jshell,  &
              kshell, lshell

              
      integer ispin, ncells

      real(dp) time_start, time_end

#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer num_loc
#endif

!      external  memory, timer

      ncells = norb / nuotot

#ifdef MPI
      call GetNodeOrbs(list_kl%nelement, Node, Nodes, list_kl_local)
#endif

      do i_list_ij = 1, list_ij%nelement
         io = list_ij%element(i_list_ij)%pair(1)
         jo = list_ij%element(i_list_ij)%pair(2)
         ri = list_ij%element(i_list_ij)%r1 
         rj = list_ij%element(i_list_ij)%r2
         rij2 = list_ij%element(i_list_ij)%dist2

         ia = iaorb(io)
         ja = iaorb(jo)

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
                    sfc_shell(jshell,ishell,js,is)%x(2)

#ifdef MPI
         do i_list_kl_local = 1, list_kl_local

            call LocalToGlobalOrb(i_list_kl_local, Node, Nodes, i_list_kl)
#else
         do i_list_kl = 1, list_kl%nelement
#endif
            ko = list_kl%element(i_list_kl)%pair(1)
            lo = list_kl%element(i_list_kl)%pair(2)
            rk = list_kl%element(i_list_kl)%r1
            rl = list_kl%element(i_list_kl)%r2
            rkl2 =  list_kl%element(i_list_kl)%dist2
            ka = iaorb(ko) 
            la = iaorb(lo)

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

            if( index_kl .le. index_ij ) then

            if(um_cut(io,ko).and.um_cut(io,lo).and.um_cut(jo,ko) &
              .and.um_cut(jo,lo) )  cycle

            if(ia.eq.ja.and.ia.eq.ka.and.ka.eq.la) cycle ! Four centers along to the same atom

            max_val2_set = (sfc_shell(lshell,kshell,ls,ks)%x(1)*rkl2 + &
                        sfc_shell(lshell,kshell,ls,ks)%x(2) )

            max_val = max_val1 + max_val2_set

            eps_temp = eri_prescreen_arg(io, jo)*eri_prescreen_arg(ko, lo)
            eps_temp = dsqrt(eps_temp)

            if(hfx_opts%DM_trunc) then 
               DM_max = 2.0d0*max( P_max(io, ko)*P_max(jo, lo), &
                        P_max(io, lo)*P_max(jo, ko) )

               IF(DM_max <= 0.0_dp) THEN
                  log10_pmax = log_zero
               ELSE
                  log10_pmax = LOG10(DM_max)
               END IF

               eps_temp = DM_max*eps_temp
            endif

                
               if( eps_temp .gt.hfx_opts%eps_schwarz) then  ! hfx_parameter%eps_schwarz) then

                  tmp_R_1 => pair_dist_radii_pgf(:,:,jshell,ishell,js,is)
                  tmp_R_2 => pair_dist_radii_pgf(:,:,lshell,kshell,ls,ks)
                  tmp_screen_pgf1 => sfc_pgf(:,:,jshell,ishell,js,is)
                  tmp_screen_pgf2 => sfc_pgf(:,:,lshell,kshell,ls,ks)

                  max_contraction_val = ispp%orbnl_contraction_coeff(ishell) * &
                                        jspp%orbnl_contraction_coeff(jshell)* &
                                        kspp%orbnl_contraction_coeff(kshell)* &
                                        lspp%orbnl_contraction_coeff(lshell)* &
                                        DM_max


                   allocate( eri_deriv(nso(l_i),nso(l_j),nso(l_k),nso(l_l),12) )
                   eri_deriv = 0.0d0

! FIXME: uncomment this call after upgrading calc_contract_deriv_eri
!                   call calc_contract_deriv_eri( deriv, cell, rcell, ri, rj, rk, rl,   &
!                                            npgfi, npgfj, npgfk, npgfl,  &
!                                            l_i, l_j, l_k, l_l,          &
!                                            ncoi, ncoj, ncok, ncol,      &
!                           ispp%orbnl_zeta(1:npgfi,ispp%orb_index(ioa)), &
!                           jspp%orbnl_zeta(1:npgfj,jspp%orb_index(joa)), &
!                           kspp%orbnl_zeta(1:npgfk,kspp%orb_index(koa)), &
!                           lspp%orbnl_zeta(1:npgfl,lspp%orb_index(loa)), &
!                           ispp%sphi(1:ncoi, ioa:ioa+nso(l_i)-1),        &
!                           jspp%sphi(1:ncoj, joa:joa+nso(l_j)-1),        &
!                           kspp%sphi(1:ncok, koa:koa+nso(l_k)-1),        &
!                           lspp%sphi(1:ncol, loa:loa+nso(l_l)-1),        &
!                           hfx_opts, eri_deriv,max_contraction_val, &
!                           max_val2_set,log10_eps_schwarz, log10_pmax,   &
!                           tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2 )
!                   call cpu_time(time_end)
   
                do nsoi = 1, nso(l_i)
                   do nsoj = 1, nso(l_j)
                      do nsok = 1, nso(l_k)
                         do nsol = 1, nso(l_l)                   

                          if(DM_max*maxval(dabs(eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2)) &
                              .gt.hfx_opts%eps_stored ) then 
                               num_a = io + nsoi-1
                               num_b = jo + nsoj-1
                               num_c = ko + nsok-1
                               num_d = lo + nsol-1
                               nao_eri_deriv(1:12) =  &
                                  eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2 

                               call  hfx_gradient_matrix( nspin, ncells, nuotot, norb, &
                                               iaorb, nua, num_a, num_b, num_c, num_d, &
                                               nao_eri_deriv, DM_tmp, Fal )
                             endif

                         enddo
                      enddo
                   enddo
                enddo

                deallocate(eri_deriv)

                endif ! Schwarz inequality
            endif
         enddo !mn
       enddo   !uv
          
      end subroutine evaluate_gradient

      subroutine hfx_gradient_matrix( nspin, ncells, nuotot, norb, &
                                      iaorb, nua, io, jo, ko, lo,  &
                                      nao_eri_deriv, DM_tmp, Fal )

      use nao2gto_data

      implicit none

!---------------------------- INPUT  VARIABLES -------------------------------------
      integer, intent(in)     :: nspin, ncells, nuotot, norb, iaorb(norb), &
                                 nua, io, jo, ko, lo
      real(dp), intent(in)    :: nao_eri_deriv(12)
      real(dp), intent(in)    :: DM_tmp(norb,norb,nspin)
!--------------------------- INOUT VARIABLES ---------------------------------------
      real(dp), intent(inout) :: Fal(1:3,nua)

!------------------------- TEMPOS, INTERNAL VARIABLES ------------------------------
      real(dp) gint_deriv(1:12)
      integer iuo, juo, kuo, luo, llo,  &
              ishell, jshell, kshell, lshell, iushell,  &
              jushell, kushell, lushell, index_ij, index_kl, &
              io_trans, jo_trans, ko_trans, lo_trans, &
              ia, ja, ka, la

      integer ispin
 
            gint_deriv = nao_eri_deriv
            iuo = io ! u is always u0
            juo = indxuo(jo)
            kuo = indxuo(ko)
            luo = indxuo(lo)
            llo = indexsc(ko, kuo, lo)
            ! lo have to trans to play with m0 to get kl and
            ! campared to ij, so there is llo

            ia=iaorb(iuo)
            ja=iaorb(juo)
            ka=iaorb(kuo)
            la=iaorb(luo)
                                                 
            iushell = subshell(iuo)
            jushell = subshell(juo)
            kushell = subshell(kuo)
            lushell = subshell(luo)
            jshell  = subshell(jo)
            lshell  = subshell(llo)

            index_ij = ncells*iushell*(iushell-1)/2 + &
                       ((jshell-1)/subshell(nuotot))*iushell + jushell 

            index_kl = ncells*kushell*(kushell-1)/2 + &
                       ((lshell-1)/subshell(nuotot))*kushell + lushell

            if( iushell .eq. jushell )   gint_deriv(1:12) = gint_deriv(1:12)*0.5d0
            if( kushell .eq. lushell )   gint_deriv(1:12) = gint_deriv(1:12)*0.5d0
            if( index_ij .eq. index_kl ) gint_deriv(1:12) = gint_deriv(1:12)*0.5d0

!! HFX
        do ispin=1,nspin
!
!         (u0vR|mR'nR")        =       (u0vR|nR"mR')
! = (v0u[-R]|m[R'-R]n[R"-R])   = (v0u[-R]|n[R"-R]m[R'-R])
! = (m0n[R"-R']|u[-R']v[R-R']) = (m0n[R"-R']|v[R-R']u[-R'])
! = (n0m[R'-R"]|u[-R"]v[R-R"]) = (n0m[R'-R"]|v[R-R"]u[-R"])
!
       ! 1.VEE(1[0]  2[H] | 3[G]  4[N])  (u0v[R]|m[R']n[R"])
       !   VEE(1[0]  2[H] | 4[N]  3[G])  (u0v[R])|n[R"]m[R'])
      Fal(1:3,ia) = Fal(1:3,ia) &
       + 0.25d0*0.25d0*gint_deriv(1:3)*DM_tmp(jo,lo,ispin)*DM_tmp(io,ko,ispin) &
       + 0.25d0*0.25d0*gint_deriv(1:3)*DM_tmp(jo,ko,ispin)*DM_tmp(io,lo,ispin)

      Fal(1:3,ja) = Fal(1:3,ja) &
       + 0.25d0*0.25d0*gint_deriv(4:6)*DM_tmp(jo,lo,ispin)*DM_tmp(io,ko,ispin) &
       + 0.25d0*0.25d0*gint_deriv(4:6)*DM_tmp(jo,ko,ispin)*DM_tmp(io,lo,ispin)

      Fal(1:3,ka) =  Fal(1:3,ka) &
       + 0.25d0*0.25d0*gint_deriv(7:9)*DM_tmp(jo,lo,ispin)*DM_tmp(io,ko,ispin) &
       + 0.25d0*0.25d0*gint_deriv(7:9)*DM_tmp(jo,ko,ispin)*DM_tmp(io,lo,ispin)

      Fal(1:3,la) =  Fal(1:3,la) &
       + 0.25d0*0.25d0*gint_deriv(10:12)*DM_tmp(jo,lo,ispin)*DM_tmp(io,ko,ispin) &
       + 0.25d0*0.25d0*gint_deriv(10:12)*DM_tmp(jo,ko,ispin)*DM_tmp(io,lo,ispin)

                      
      ! 2.VEE(2[0] 1[-H]| 3[G-H] 4[N-H])  (v0u[-R]|m[R'-R]n[R"-R])
      !   VEE(2[0] 1[-H]| 4[N-H] 3[G-H])  (v0u[-R]|n[R"-R]m[R'-R])

         io_trans = indexsc( jo, juo, io )
         ko_trans = indexsc( jo, juo, ko )
         lo_trans = indexsc( jo, juo, lo )

      Fal(1:3,ja) = Fal(1:3,ja) &
       + 0.25d0*0.25d0*gint_deriv(4:6)*DM_tmp(io_trans,lo_trans,ispin) &
         *DM_tmp(juo,ko_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(4:6)*DM_tmp(io_trans,ko_trans,ispin) &
         *DM_tmp(juo,lo_trans,ispin)

      Fal(1:3,ia) = Fal(1:3,ia) &
       + 0.25d0*0.25d0*gint_deriv(1:3)*DM_tmp(io_trans,lo_trans,ispin) &
         *DM_tmp(juo,ko_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(1:3)*DM_tmp(io_trans,ko_trans,ispin) &
         *DM_tmp(juo,lo_trans,ispin)

      Fal(1:3,ka) = Fal(1:3,ka) &
       + 0.25d0*0.25d0*gint_deriv(7:9)*DM_tmp(io_trans,lo_trans,ispin) &
         *DM_tmp(juo,ko_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(7:9)*DM_tmp(io_trans,ko_trans,ispin) &
         *DM_tmp(juo,lo_trans,ispin)

      Fal(1:3,la) = Fal(1:3,la) &
       + 0.25d0*0.25d0*gint_deriv(10:12)*DM_tmp(io_trans,lo_trans,ispin) &
         *DM_tmp(juo,ko_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(10:12)*DM_tmp(io_trans,ko_trans,ispin) &
         *DM_tmp(juo,lo_trans,ispin)

      ! 3.VEE(3[0]  4[N-G] | 1[-G] 2[H-G]) (m0n[R"-R']|u[-R']v[R-R'])
      !   VEE(3[0]  4[N-G] |2[H-G] 1[-G] ) (m0n[R"-R']|v[R-R']u[-R'])

         io_trans = indexsc( ko, kuo, io )
         jo_trans = indexsc( ko, kuo, jo )
         lo_trans = indexsc( ko, kuo, lo )

      Fal(1:3,ka) = Fal(1:3,ka) &
       + 0.25d0*0.25d0*gint_deriv(7:9)*DM_tmp(lo_trans,jo_trans,ispin) &
         *DM_tmp(kuo,io_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(7:9)*DM_tmp(lo_trans,io_trans,ispin) &
         *DM_tmp(kuo,jo_trans,ispin)

      Fal(1:3,ia) = Fal(1:3,ia) &
       + 0.25d0*0.25d0*gint_deriv(1:3)*DM_tmp(lo_trans,jo_trans,ispin) &
         *DM_tmp(kuo,io_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(1:3)*DM_tmp(lo_trans,io_trans,ispin) &
         *DM_tmp(kuo,jo_trans,ispin)

      Fal(1:3,ja) = Fal(1:3,ja) &
       + 0.25d0*0.25d0*gint_deriv(4:6)*DM_tmp(lo_trans,jo_trans,ispin) &
         *DM_tmp(kuo,io_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(4:6)*DM_tmp(lo_trans,io_trans,ispin) &
         *DM_tmp(kuo,jo_trans,ispin)

      Fal(1:3,la) = Fal(1:3,la) &
       + 0.25d0*0.25d0*gint_deriv(10:12)*DM_tmp(lo_trans,jo_trans,ispin) &
         *DM_tmp(kuo,io_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(10:12)*DM_tmp(lo_trans,io_trans,ispin) &
         *DM_tmp(kuo,jo_trans,ispin)
      
! 4.VEE(4[0]  3[G-N] | 1[-N] 2[H-N])  (n0m[R'-R"]|u[-R"]v[R-R"])
!   VEE(4[0]  3[G-N] | 2[H-N] 1[-N])  (n0m[R'-R"]|v[R-R"]u[-R"])

         io_trans = indexsc( lo, luo, io )
         jo_trans = indexsc( lo, luo, jo )
         ko_trans = indexsc( lo, luo, ko )

      Fal(1:3,la) = Fal(1:3,la) &
       + 0.25d0*0.25d0*gint_deriv(10:12)*DM_tmp(ko_trans,jo_trans,ispin) &
         *DM_tmp(luo,io_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(10:12)*DM_tmp(ko_trans,io_trans,ispin) &
         *DM_tmp(luo,jo_trans,ispin)

      Fal(1:3,ka) = Fal(1:3,ka) &
       + 0.25d0*0.25d0*gint_deriv(7:9)*DM_tmp(ko_trans,jo_trans,ispin) &
         *DM_tmp(luo,io_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(7:9)*DM_tmp(ko_trans,io_trans,ispin) &
         *DM_tmp(luo,jo_trans,ispin)

      Fal(1:3,ja) = Fal(1:3,ja) &
       + 0.25d0*0.25d0*gint_deriv(4:6)*DM_tmp(ko_trans,jo_trans,ispin) &
         *DM_tmp(luo,io_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(4:6)*DM_tmp(ko_trans,io_trans,ispin) &
         *DM_tmp(luo,jo_trans,ispin)

      Fal(1:3,ia) = Fal(1:3,ia) &
       + 0.25d0*0.25d0*gint_deriv(1:3)*DM_tmp(ko_trans,jo_trans,ispin) &
         *DM_tmp(luo,io_trans,ispin) &
       + 0.25d0*0.25d0*gint_deriv(1:3)*DM_tmp(ko_trans,io_trans,ispin) &
         *DM_tmp(luo,jo_trans,ispin)


        enddo

      end subroutine hfx_gradient_matrix

! ***** GRADIENT MODULE END *****

end module nao2gto_hfx_force
