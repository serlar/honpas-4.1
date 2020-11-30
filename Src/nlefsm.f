! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_nlefsm

      use precision,     only : dp

      implicit none

      public :: nlefsm, nlefsm_SO_off

      private

      
      CONTAINS

      subroutine nlefsm( scell, nua, na, isa, xa, indxua,
     .                   maxnh, maxnd, lasto, lastkb, iphorb, 
     .                   iphKB, numd, listdptr, listd, numh, 
     .                   listhptr, listh, nspin, Dscf, Enl, 
     .                   fa, stress, H , matrix_elements_only)
C *********************************************************************
C Calculates non-local (NL) pseudopotential contribution to total 
C energy, atomic forces, stress and hamiltonian matrix elements.
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, June 1997.
C **************************** INPUT **********************************
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer isa(na)          : Species index of each atom
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C integer indxua(na)       : Index of equivalent atom in unit cell
C integer maxnh            : First dimension of H and listh
C integer maxnd            : Maximum number of elements of the
C                            density matrix
C integer lasto(0:na)      : Position of last orbital of each atom
C integer lastkb(0:na)     : Position of last KB projector of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom,
C                            where no=lasto(na)
C integer iphKB(nokb)      : Index of each KB projector in its atom,
C                            where nokb=lastkb(na)
C integer numd(nuo)        : Number of nonzero elements of each row of the
C                            density matrix
C integer listdptr(nuo)    : Pointer to the start of each row (-1) of the
C                            density matrix
C integer listd(maxnd)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer numh(nuo)        : Number of nonzero elements of each row of the
C                            hamiltonian matrix
C integer listhptr(nuo)    : Pointer to the start of each row (-1) of the
C                            hamiltonian matrix
C integer listh(maxnh)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer nspin            : Number of spinor components (really min(nspin,2))
C                            If computing only matrix elements, it
C                            can be set to 1.
C logical matrix_elements_only:
C integer Dscf(maxnd,nspin): Density matrix. Not touched if computing
C                            only matrix elements.
C ******************* INPUT and OUTPUT *********************************
C real*8 fa(3,na)          : NL forces (added to input fa)
C real*8 stress(3,3)       : NL stress (added to input stress)
C real*8 H(maxnh,nspin)    : NL Hamiltonian (added to input H)
C **************************** OUTPUT *********************************
C real*8 Enl               : NL energy
C *********************************************************************
C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,  only : GlobalToLocalOrb
      use atm_types,     only : nspecies
      use atomlist,      only : in_kb_orb_u_range
      use atmfuncs,      only : rcut, epskb, orb_gindex, kbproj_gindex
      use atmfuncs,      only : nofis, nkbfis
      use chemical,      only : is_floating
      use neighbour,     only : iana=>jan, r2ki=>r2ij, xki=>xij
      use neighbour,     only : mneighb, reset_neighbour_arrays
      use alloc,         only : re_alloc, de_alloc
      use m_new_matel,   only : new_matel

      integer, intent(in) ::
     .   maxnh, na, maxnd, nspin, nua

      integer, intent(in)  ::
     .  indxua(na), iphKB(*), iphorb(*), isa(na),  
     .  lasto(0:na), lastkb(0:na), listd(maxnd), listh(maxnh),
     .  numd(*), numh(*), listdptr(*), listhptr(*)

      real(dp), intent(in) :: scell(3,3), Dscf(maxnd,nspin),
     .                        xa(3,na)
      real(dp), intent(inout) :: fa(3,nua), stress(3,3)
      real(dp), intent(inout) :: H(maxnh,nspin)
      real(dp), intent(out)   :: Enl
      logical, intent(in)     :: matrix_elements_only

      real(dp) ::   volcel
      external ::   timer, volcel

C Internal variables ................................................
C maxno  = maximum number of basis orbitals overlapping a KB projector
! This should be estimated beforehand, to avoid checks,
! or a "guard" region implemented for a single check at the end
      
      integer, save ::  maxno = 500
  
      integer
     .  ia, ikb, ina, ind, ino,
     .  io, iio, ioa, is, ispin, ix, ig, kg,
     .  j, jno, jo, jx, ka, ko, koa, ks, kua,
     .  nkb, nna, nno, no, nuo, nuotot, maxkba
      integer :: natoms_k_over, max_nno_used

      integer, dimension(:), pointer :: iano, iono

      real(dp)
     .  Cijk, epsk, fik, rki, rmax, rmaxkb, rmaxo, 
     .  Sik, Sjk, volume

      real(dp), dimension(:), pointer :: Di, Vi
      real(dp), dimension(:,:), pointer :: Ski, xno
      real(dp), dimension(:,:,:), pointer :: grSki

      logical ::   within
      logical, dimension(:), pointer ::  listed, listedall
      !
      real(dp), allocatable :: rorbmax(:), rkbmax(:)
C ......................

C Start time counter
      call timer( 'nlefsm', 1 )

C Find unit cell volume
      volume = volcel( scell ) * nua / na

C Find maximum range and maximum number of KB projectors
      maxkba = 0

      allocate(rorbmax(nspecies),rkbmax(nspecies))
      do is = 1, nspecies

         ! Species orbital range
         rorbmax(is) = 0.0_dp
         do io = 1, nofis(is)
            rorbmax(is) = max(rorbmax(is), rcut(is,io))
         enddo
         
         ! Species KB range
         io = nkbfis(is)
         rkbmax(is) = 0.0_dp
         do ikb = 1, io
            rkbmax(is) = max(rkbmax(is), rcut(is,-ikb))
         enddo
         maxkba = max(maxkba,io)

      enddo
      rmaxo = maxval(rorbmax(1:nspecies))
      rmaxkb = maxval(rkbmax(1:nspecies))
      ! Calculate max extend
      rmax = rmaxo + rmaxkb

C Initialize arrays Di and Vi only once

      no = lasto(na)
      nuotot = lasto(nua)
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)

C Allocate local memory

      nullify( Vi )
      call re_alloc( Vi, 1, no, 'Vi', 'nlefsm' )
      Vi(1:no) = 0.0_dp         ! OK. Later on, any non-zero elements
                                ! will be zero-ed out explicitly
      nullify( listed )
      call re_alloc( listed, 1, no, 'listed', 'nlefsm' )
      listed(1:no) = .false.
      nullify( listedall )
      call re_alloc( listedall, 1, no, 'listedall', 'nlefsm' )
      listedall(1:no) = .false.

      if (.not. matrix_elements_only) then
         nullify( Di )
         call re_alloc( Di, 1, no, 'Di', 'nlefsm' )
         Di(1:no) = 0.0_dp
      endif

      Enl = 0.0d0

!     Make list of all orbitals needed for this node
      
      do io = 1,nuo
        ! we need this process's orbitals...
        call LocalToGlobalOrb(io,Node,Nodes,iio)
        listedall(iio) = .true.
        ! ... and those with which they interact
        do j = 1,numh(io)
          jo = listh(listhptr(io)+j)
          listedall(jo) = .true.
        enddo
      enddo

C Allocate local arrays that depend on saved parameters
      nullify( iano )
      call re_alloc( iano, 1, maxno, 'iano', 'nlefsm' )
      nullify( iono )
      call re_alloc( iono, 1, maxno, 'iono', 'nlefsm' )
      nullify( xno )
      call re_alloc( xno, 1, 3, 1, maxno, 'xno',  'nlefsm' )
      nullify( Ski )
      call re_alloc( Ski, 1, maxkba, 1, maxno, 'Ski', 'nlefsm' )
      nullify( grSki )
      call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno, 'grSki',
     &               'nlefsm' )

C     Initialize neighb subroutine
      call mneighb( scell, rmax, na, xa, 0, 0, nna )

! Loop on atoms with KB projectors
! All processes will be doing this loop over atoms.
! This is one reason for non-scalability
!
!        And what happens if there is no supercell?
!        How do we count out-of-unit-cell interactions?
!        ... they are automatically accounted for, in
!        the same way as the Hmu_nu terms themselves.
!
      natoms_k_over = 0
      max_nno_used = 0
      do ka = 1,na
!        Only the atoms within the proper
!        distance of a unit cell orbital (in our process) should
!        be considered, not the whole supercell.
!        This array was initialized in hsparse
         
        if (.not. in_kb_orb_u_range(ka)) CYCLE
        
        ks = isa(ka)
        ! Cycle also if ghost-orbital species...
        if (is_floating(ks)) CYCLE
        
        kua = indxua(ka)  ! Used only if forces and energies are comp.

C       Find neighbour atoms
        call mneighb( scell, rmax, na, xa, ka, 0, nna )

        nno = 0
        do ina = 1,nna
          rki = sqrt(r2ki(ina))
          ia = iana(ina)
          is = isa(ia)
          !     Early exit if too far
          !     This duplicates the test in hsparse...
          if (rki - rkbmax(ks) - rorbmax(is) > 0.d0) CYCLE

          ! Loop over orbitals close enough to overlap
          do io = lasto(ia-1)+1,lasto(ia)

C           Only calculate if needed locally in our MPI process
             if (.not. listedall(io)) CYCLE
             
             ioa = iphorb(io)
             ! rki_minus_rc_orb= rki - rcut(is,ioa)

             ! Find if orbital is within range
             ! This can be done with rkbmax(ks):
             within = (rki-rkbmax(ks)) < rcut(is,ioa)
             if (.not. within) CYCLE
              
!             Find overlap between neighbour orbitals and KB projectors

             if (nno.eq.maxno) call increase_maxno()
              
              nno = nno + 1  ! Update number of overlaps to keep
              iono(nno) = io 
              iano(nno) = ia
              do ix = 1,3
                 xno(ix,nno) = xki(ix,ina)
              enddo

! For each overlap family we keep the individual
! KB-orb matrix elements
! This will store some zeros sometimes, as some
! of the KBs might not actually overlap
! We could re-check the distances...
! Not worth it, as then we would have different
! numbers of matrix elements for different orbitals, and
! the bookeeping would get messy
              
              
              ikb = 0
              ig = orb_gindex(is,ioa)
              do ko = lastkb(ka-1)+1,lastkb(ka)
                 koa = iphKB(ko)
                 ! if ( rki_minus_rc_orb > rcut(ks,koa) CYCLE
                 ikb = ikb + 1
                 ! epsk_sqrt = sqrt(epskb(ks,koa))
                 kg = kbproj_gindex(ks,koa)
                 call new_MATEL( 'S', kg, ig, xki(1:3,ina),
     &                Ski(ikb,nno), grSki(1:3,ikb,nno) )
              enddo

           enddo ! loop over orbitals

        enddo ! loop over neighbor atoms

!     Now we check which of the overlaps of our atom's KB's involve
!     two orbitals: one in the unit cell, and handled by our process,
!     and the other unrestricted
        
        max_nno_used = max(max_nno_used, nno)
        do ino = 1,nno    ! loop over overlaps
          ia = iano(ino)
          if (ia > nua) CYCLE  ! We want the 1st orb to be in the unit cell

          io = iono(ino)
          ! Note that if ia is in the unit cell, io is <= nuo,
          ! so that this call makes sense
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio == 0) CYCLE

          !  Scatter filter of desired matrix elements
          do j = 1,numh(iio)
             ind = listhptr(iio)+j
             jo = listh(ind)
             listed(jo) = .true.
             if (.not. matrix_elements_only) then
                do ispin = 1, nspin ! Both spins add up... (nspin=spin%spinor here)
                   Di(jo) = Di(jo) + Dscf(ind,ispin)
                enddo
             endif
          enddo

! Find matrix elements with other neighbour orbitals
! Note that several overlaps might contribute to the
! same matrix element, hence the additions above (Dscf) and below (H)
          
          do jno = 1,nno
             jo = iono(jno)
             ! Check whether there is H_io_jo...
             if (.not. listed(jo)) CYCLE 

! Loop on KB projectors again. Note that ikb and ko run
! in step. ko is only needed for the Epskb factor.
! maybe we can store it with the value of the projector.

             ikb = 0
             do ko = lastkb(ka-1)+1,lastkb(ka)
                ikb = ikb + 1
                koa = iphKB(ko)
                epsk = epskb(ks,koa)
                Sik = Ski(ikb,ino)
                Sjk = Ski(ikb,jno)
                Vi(jo) = Vi(jo) + epsk * Sik * Sjk
                ! We should distinguish "energy-only" and
                ! "forces-and-stress"
                if (.not. matrix_elements_only) then
                   Cijk = Di(jo) * epsk
                   Enl = Enl + Cijk * Sik * Sjk
                   do ix = 1,3
                      fik = 2.d0 * Cijk * Sjk * grSki(ix,ikb,ino)
                      fa(ix,ia)  = fa(ix,ia)  - fik
                      fa(ix,kua) = fa(ix,kua) + fik
                      do jx = 1,3
                         stress(jx,ix) = stress(jx,ix) +
     &                        xno(jx,ino) * fik / volume
                      enddo
                   enddo
                endif

             enddo
          
          enddo ! loop over second orbitals

C         Pick up contributions to H and restore Di and Vi
          do j = 1,numh(iio)
             ind = listhptr(iio)+j
             jo = listh(ind)
             do ispin = 1,nspin  ! Only diagonal parts
                H(ind,ispin) = H(ind,ispin) + Vi(jo)
             enddo
             Vi(jo) = 0.0d0     ! See initial zero-out at top
             listed(jo) = .false.
             if (.not. matrix_elements_only) Di(jo) = 0.0d0
          enddo

       enddo  ! loop over 1st orbitals
       natoms_k_over = natoms_k_over + 1
      enddo   ! loop over atoms holding KB projectors

      if (Node == 0) then
         ! For future diagnostics
         ! Currently only the root process outputs info
         write(6,"(a,2i8)")
     $     "No. of atoms with KB's overlaping orbs in proc 0." //
     $     " Max # of overlaps:", natoms_k_over, max_nno_used
      endif
      
C     Deallocate local memory
!      call new_MATEL( 'S', 0, 0, 0, 0, xki, Ski, grSki )
      call reset_neighbour_arrays( )
      call de_alloc( grSki, 'grSki', 'nlefsm' )
      call de_alloc( Ski, 'Ski', 'nlefsm' )
      call de_alloc( xno, 'xno', 'nlefsm' )
      call de_alloc( iono, 'iono', 'nlefsm' )
      call de_alloc( iano, 'iano', 'nlefsm' )
      call de_alloc( listedall, 'listedall', 'nlefsm' )
      call de_alloc( listed, 'listed', 'nlefsm' )
      call de_alloc( Vi, 'Vi', 'nlefsm' )
      if (.not. matrix_elements_only) then
         call de_alloc( Di, 'Di', 'nlefsm' )
      endif
      
      deallocate(rkbmax,rorbmax)

      call timer( 'nlefsm', 2 )
      
      CONTAINS
      
      subroutine increase_maxno()
      
! if too small then increase array sizes
      maxno = maxno + 10
      call re_alloc( iano, 1, maxno, 'iano', 'nlefsm',
     &        .true. )
      call re_alloc( iono, 1, maxno, 'iono', 'nlefsm',
     &        .true. )
      call re_alloc( xno, 1, 3, 1, maxno, 'xno',  'nlefsm',
     &        .true. )
      call re_alloc( Ski, 1, maxkba, 1, maxno, 'Ski',
     &        'nlefsm', .true. )
      call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno,
     &        'grSki', 'nlefsm', .true. )

      end subroutine increase_maxno
      
      end subroutine nlefsm

C nlefsm_SO_off calculates the KB elements to the total Hamiltonian 
C (including the SO contribution)
C when Off-Site Spin Orbit is included in the calculation 
      subroutine nlefsm_SO_off( scell, nua, na, isa, xa, indxua,
     .                      maxnh, maxnd, lasto, lastkb, iphorb, 
     .                      iphKB, numd, listdptr, listd, numh, 
     .                      listhptr, listh, nspin, Enl, Enl_offsiteSO, 
     .                      fa, stress, H0 , H0_off,
     &                      matrix_elements_only,
     &                      ill_defined_sr_so_split)


C *********************************************************************
C Calculates non-local (NL) pseudopotential contribution to total 
C energy, atomic forces, stress and hamiltonian matrix elements.
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, June 1997.
C Modified by R. Cuadrado and J. I. Cerda for the 
C off-site Spin-Orbit, January 2017.
C **************************** INPUT **********************************
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer isa(na)          : Species index of each atom
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C integer indxua(na)       : Index of equivalent atom in unit cell
C integer maxnh            : First dimension of H and listh
C integer maxnd            : Maximum number of elements of the
C                            density matrix
C integer lasto(0:na)      : Position of last orbital of each atom
C integer lastkb(0:na)     : Position of last KB projector of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom,
C                            where no=lasto(na)
C integer iphKB(nokb)      : Index of each KB projector in its atom,
C                            where nokb=lastkb(na)
C integer numd(nuo)        : Number of nonzero elements of each row of the
C                            density matrix
C integer listdptr(nuo)    : Pointer to the start of each row (-1) of the
C                            density matrix
C integer listd(maxnd)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer numh(nuo)        : Number of nonzero elements of each row of the
C                            hamiltonian matrix
C integer listhptr(nuo)    : Pointer to the start of each row (-1) of the
C                            hamiltonian matrix
C integer listh(maxnh)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer nspin            : Number of spin_grid components (really min(nspin,4))
C                            If computing only matrix elements, it
C                            can be set to 1.
C logical matrix_elements_only:
C integer Dscf(maxnd,nspin): Density matrix. Not touched if computing
C                            only matrix elements.
C ******************* INPUT and OUTPUT *********************************
C real*8 fa(3,na)          : NL forces (added to input fa)
C real*8 stress(3,3)       : NL stress (added to input stress)
C real*8 H0(maxnh,nspin)    : NL Hamiltonian (added to input H)
C complex*16 H0_off(maxnh,nspin) : NL off-site Hamiltonian (added to input H)
C **************************** OUTPUT *********************************
C real*8 Enl               : NL energy
C *********************************************************************
C
C  Modules
C
      use parallel,        only : Node, Nodes
      use parallelsubs,    only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,    only : GlobalToLocalOrb
      use atmfuncs,        only : rcut, orb_gindex, kbproj_gindex
      use atmfuncs,        only : epskb
      use neighbour,       only : iana=>jan, r2ki=>r2ij, xki=>xij
      use neighbour,       only : mneighb, reset_neighbour_arrays
      use alloc,           only : re_alloc, de_alloc
      use m_new_matel,     only : new_matel
      use atm_types,       only : species_info, species
      use sparse_matrices, only : Dscf, xijo
      use siesta_options,  only : split_sr_so
      use chemical,        only : is_floating
      use atomlist,        only : in_kb_orb_u_range
      
      use fdf
      
      integer, intent(in) ::
     .   maxnh, na, maxnd, nspin, nua

      integer, intent(in)  ::
     .  indxua(na), iphKB(*), iphorb(*), isa(na),  
     .  lasto(0:na), lastkb(0:na), listd(maxnd), listh(maxnh),
     .  numd(*), numh(*), listdptr(*), listhptr(*)

      real(dp), intent(in)    :: scell(3,3),  xa(3,na)
      real(dp), intent(inout) :: fa(3,nua), stress(3,3)
      real(dp), intent(inout) :: H0(maxnh)
      complex(dp), intent(inout) :: H0_off(maxnh,4)
                                           
      real(dp), intent(out)   :: Enl, Enl_offsiteSO 
      logical, intent(in)     :: matrix_elements_only
      logical, intent(out)    :: ill_defined_sr_so_split

      real(dp) ::   volcel
      external ::   timer, volcel

C Internal variables ................................................
C maxno  = maximum number of basis orbitals overlapping a KB projector

      integer, save ::  maxno = 2000
  
      integer
     .  ia, ikb, ina, ind, ino, indt,
     .  io, iio, ioa, is, ispin, ix, ig, kg,
     .  j, jno, jo, jx, ka, ko, koa, ks, kua,
     .  nkb, nna, nno, no, nuo, nuotot, maxkba

      integer :: l, koa1, koa2, i 

      integer, dimension(:), pointer :: iano, iono

      real(dp)
     .  Cijk, fik, rki, rmax, rmaxkb, rmaxo, 
     .  volume, CVj, epsk(2), Vit  

      real(dp)           :: r1(3), r2(3)

      real(dp), dimension(:), pointer :: Di, Vi ! This is Vion
      real(dp), dimension(:,:), pointer :: Ski, xno
      real(dp), dimension(:,:,:), pointer :: grSki

      complex(dp)          :: V_sot(2,2), F_so(3,2,2)
      complex(dp), pointer :: V_so(:,:,:), Ds(:,:,:)


      logical ::   within
      logical, dimension(:), pointer ::  listed, listedall

      complex(dp) :: E_offsiteSO(4)
      
      type(species_info), pointer        :: spp

      integer :: nd, ndn
      real(dp) :: Vit_saved

C ------------------------------------------------------------

C Start time counter
      call timer( 'nlefsm', 1 )

C Find unit cell volume
      volume = volcel( scell ) * nua / na

C Find maximum range
      rmaxo = 0.0d0
      rmaxkb = 0.0d0
      do ia = 1,na
        is = isa(ia)
        do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphKB(ikb)
          rmaxkb = max( rmaxkb, rcut(is,ioa) )
        enddo
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rmaxo = max( rmaxo, rcut(is,ioa) )
        enddo
      enddo
      rmax = rmaxo + rmaxkb

C Initialize arrays Di and Vi only once
      no = lasto(na)
      nuotot = lasto(nua)

      call GetNodeOrbs(nuotot,Node,Nodes,nuo)

C Allocate local memory

      nullify( Vi ) ! This is Vion
      call re_alloc( Vi, 1, no, 'Vi', 'nlefsm_SO_off' )
      Vi(1:no) = 0.0_dp

      nullify( listed )
      call re_alloc( listed, 1, no, 'listed', 'nlefsm_SO_off' )
      listed(1:no) = .false.
      nullify( listedall )
      call re_alloc( listedall, 1, no, 'listedall', 'nlefsm_SO_off' )
      listedall(1:no) = .false.

      allocate( V_so(2,2,no) )

      if (.not. matrix_elements_only) then
       nullify( Di ) 
       call re_alloc( Di, 1, no, 'Di', 'nlefsm_SO_off' )
       Di(1:no) = 0.0_dp
       allocate( Ds(2,2,no) )
      endif

C Make list of all orbitals needed for this node

      do io = 1,nuo
        call LocalToGlobalOrb(io,Node,Nodes,iio)
        listedall(iio) = .true.
        do j = 1,numh(io)
          jo = listh(listhptr(io)+j)
          listedall(jo) = .true.
        enddo
      enddo

C Find maximum number of KB projectors of one atom = maxkba
      maxkba = 0
      do ka = 1,na
        nkb = lastkb(ka) - lastkb(ka-1)
        maxkba = max(maxkba,nkb)
      enddo

C Allocate local arrays that depend on saved parameters
      nullify( iano )
      call re_alloc( iano, 1, maxno, 'iano', 'nlefsm_SO_off' )
      nullify( iono )
      call re_alloc( iono, 1, maxno, 'iono', 'nlefsm_SO_off' )
      nullify( xno )
      call re_alloc( xno, 1, 3, 1, maxno, 'xno',  'nlefsm_SO_off' )
      nullify( Ski )
      call re_alloc( Ski, 1, maxkba, 1, maxno, 'Ski','nlefsm_SO_off')
      nullify( grSki )
      call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno, 'grSki',
     &               'nlefsm_SO_off' )

C     Initialize neighb subroutine
      call mneighb( scell, rmax, na, xa, 0, 0, nna )

      ill_defined_sr_so_split = .false.
      nd= 0
      ndn= 0
      Enl = 0.0d0
      E_offsiteSO(1:4)=dcmplx(0.0d0,0.0d0)
      Enl_offsiteSO = 0.0d0
C     Loop on atoms with KB projectors      


      do ka = 1,na              ! Supercell atoms

       ! Use info from hsparse to remove "far" atoms, as in nlefsm
       if (.not. in_kb_orb_u_range(ka)) CYCLE

       ks = isa(ka)             ! Specie index of atom ka
       if (is_floating(ks)) CYCLE

       kua = indxua(ka)         ! Equivalent atom in the UC
       nkb = lastkb(ka) - lastkb(ka-1) ! number of KB projs of atom ka

C      Find neighbour atoms
       call mneighb( scell, rmax, na, xa, ka, 0, nna )

C      Find neighbour orbitals
       Ski(:,:) = 0.0_dp
       nno = 0
       iano(:)=0
       iono(:)=0
       do ina = 1,nna  ! Neighbour atoms
        ia = iana(ina) ! Atom index of ina (the neighbour to ka)
        is = isa(ia)   ! Specie index of atom ia
        rki = sqrt(r2ki(ina)) ! Square distance

        ! Use this to further filter
        !!     if (rki - rkbmax(ks) - rorbmax(is) > 0.d0) CYCLE
        do io = lasto(ia-1)+1,lasto(ia) ! Orbitals of atom ia

C        Only calculate if needed locally

        if (listedall(io)) then
          ioa = iphorb(io)  ! Orbital index of orbital io in 
C                             neighbour atom ina

C         Find if orbital is within range
          within = .false.
          do ko = lastkb(ka-1)+1,lastkb(ka)
           koa = iphKB(ko)
           if ( rki .lt. rcut(is,ioa)+rcut(ks,koa) ) 
     &            within = .true.
          enddo

C         Find overlap between neighbour orbitals and KB projectors
          if (within) then
C          Check maxno - if too small then increase array sizes
           if (nno.eq.maxno) then
            maxno = maxno + 100
            call re_alloc( iano, 1, maxno, 'iano', 'nlefsm_SO_off',
     &                    .true. )
            call re_alloc( iono, 1, maxno, 'iono', 'nlefsm_SO_off',
     &                    .true. )
            call re_alloc( xno, 1, 3, 1, maxno,'xno','nlefsm_SO_off',
     &                    .true. )
            call re_alloc( Ski, 1, maxkba, 1, maxno, 'Ski',
     &                    'nlefsm_SO_off', .true. )
            call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno,
     &                    'grSki', 'nlefsm_SO_off', .true. )
           endif
           nno = nno + 1  ! Number of neighbour orbitals
           iono(nno) = io ! io orbital of atom ina (neighbour to ka)
           iano(nno) = ia ! atom index of the neighbour ina
           xno(1:3,nno) = xki(1:3,ina)
           ikb = 0
           do ko = lastkb(ka-1)+1,lastkb(ka) !Generic positions of kb's
            ikb = ikb + 1
            ioa = iphorb(io)
            koa = iphKB(ko)            ! koa = -ikb 

            !if ( koa.ne.-ikb ) then
            ! write(6,*) 'koa ERROR: koa,ikb=',koa,ikb
            ! stop
            !endif
            kg = kbproj_gindex(ks,koa)
            ig = orb_gindex(is,ioa)
            call new_MATEL( 'S', kg, ig, xki(1:3,ina),
     &                     Ski(ikb,nno), grSki(1:3,ikb,nno) )
           enddo
          endif  ! Within

         endif    
        enddo    ! neighbour AO
       enddo     ! neighbour atoms

C----- Loop on neighbour orbitals
       do ino = 1,nno
        io = iono(ino)
        ia = iano(ino)

        call GlobalToLocalOrb(io,Node,Nodes,iio)

        if (iio.gt.0) then
C        Valid orbital
         if (ia .le. nua) then
          if (.not. matrix_elements_only) then
           !Scatter density matrix row of orbital io
           Ds(1:2,1:2,1:no) = dcmplx(0.0d0,0.0d0) 
           do j = 1,numh(iio)
            ind = listhptr(iio)+j  ! jptr
            jo = listh(ind)       ! j
            Di(jo) = 0.0_dp
            do ispin = 1,min(2,nspin) ! Only diagonal parts 
             Di(jo) = Di(jo) + Dscf(ind,ispin)
            enddo
            ! Should we *add to* Ds, as above for Di?
            Ds(1,1,jo) = dcmplx(Dscf(ind,1), Dscf(ind,5)) ! D(ju,iu)
            Ds(2,2,jo) = dcmplx(Dscf(ind,2), Dscf(ind,6))  ! D(jd,id)
            Ds(1,2,jo) = dcmplx(Dscf(ind,3), Dscf(ind,4))  ! D(ju,id)
            Ds(2,1,jo) = dcmplx(Dscf(ind,7),-Dscf(ind,8))  ! D(jd,iu)
           enddo
          endif

C-------- Scatter filter of desired matrix elements
          do j = 1,numh(iio)
           jo = listh(listhptr(iio)+j)
           listed(jo) = .true.
          enddo

C  Loading V_ion/V_so
          Vi(1:no) = 0.0_dp
          V_so(1:2,1:2,1:no) = dcmplx(0.0d0,0.0d0)

C-------- Find matrix elements with other neighbour orbitals
          do jno = 1,nno
           jo = iono(jno)

           if (listed(jo)) then

C---------- Loop on KB projectors
            ko  = lastkb(ka-1)
            KB_loop: do
             koa = -iphKB(ko+1)
             spp => species(ks)
             l = spp%pj_l(koa)

c----------- Compute Vion
             if ( l.eq.0 ) then
              epsk(1) = epskb(ks,koa)
              Vit_saved = epsk(1) * Ski(koa,ino) * Ski(koa,jno)
              if (split_sr_so) then
                 Vit = Vit_saved
              else
                 ! Move 'ionic' part to V_so diagonal
                 V_so(1,1,jo)= V_so(1,1,jo) + Vit_saved
                 V_so(2,2,jo)= V_so(2,2,jo) + Vit_saved
                 Vit = 0.0_dp
              endif
              Vi(jo) = Vi(jo) + Vit
              if (.not. matrix_elements_only) then
               if (split_sr_so) then
                  Enl = Enl +  Di(jo) * Vit
               else
                  ! Move energy contribution to "Eso"
                  E_offsiteSO(1) = E_offsiteSO(1) + Vit_saved*Ds(1,1,jo) 
                  E_offsiteSO(2) = E_offsiteSO(2) + Vit_saved*Ds(2,2,jo) 
               endif
               CVj  = epsk(1) * Ski(koa,jno)
               Cijk = 2.0_dp * Di(jo) * CVj
               do ix = 1,3
                fik = Cijk * grSki(ix,koa,ino)
                fa(ix,ia)  = fa(ix,ia)  - fik
                fa(ix,kua) = fa(ix,kua) + fik
                do jx = 1,3
                stress(jx,ix) = stress(jx,ix) +
     &                           xno(jx,ino) * fik / volume
                enddo
               enddo
              endif
              ko = ko + 1 

c----------- Compute Vion from j+/-1/2 and V_so
             else  ! l /= 0
              ! First proj of l-1/2 block
              koa1 = -iphKB(ko+1)
              ! Last proj of l+1/2 block
              koa2 = -iphKB(ko+2*(2*l+1))
              epsk(1) = epskb(ks,koa1)
              epsk(2) = epskb(ks,koa2)
              if (epsk(1)*epsk(2) < 0) then
                 ill_defined_sr_so_split = .true.
              endif

              call calc_Vj_offsiteSO( l, epsk, Ski(koa1:koa2,ino), 
     &                       Ski(koa1:koa2,jno), grSki(:,koa1:koa2,ino),
     &                       grSki(:,koa1:koa2,jno), Vit, V_sot, F_so,
     &                       split_sr_so)
              Vi(jo) = Vi(jo) + Vit
              V_so(1:2,1:2,jo)= V_so(1:2,1:2,jo) + V_sot(1:2,1:2)
         
c------------ Forces & SO contribution to E_NL
              if (.not. matrix_elements_only) then     
               Enl = Enl +  Di(jo) * Vit

               E_offsiteSO(1) = E_offsiteSO(1) + V_sot(1,1)*Ds(1,1,jo) ! V(iu,ju)*D(ju,iu)
               E_offsiteSO(2) = E_offsiteSO(2) + V_sot(2,2)*Ds(2,2,jo) ! V(id,jd)*D(jd,id)
               E_offsiteSO(3) = E_offsiteSO(3) + V_sot(1,2)*Ds(2,1,jo) ! V(iu,jd)*D(jd,iu)
               E_offsiteSO(4) = E_offsiteSO(4) + V_sot(2,1)*Ds(1,2,jo) ! V(id,ju)*D(ju,id)

               ! These forces include both the 'ion' and 'SO' parts
               do ix = 1,3
                fik = 2.0_dp*dreal(Ds(1,1,jo)*F_so(ix,1,1) +
     &                             Ds(2,2,jo)*F_so(ix,2,2) +
     &                             Ds(2,1,jo)*F_so(ix,1,2) +
     &                             Ds(1,2,jo)*F_so(ix,2,1) )
                fa(ix,ia)  = fa(ix,ia)  - fik
                fa(ix,kua) = fa(ix,kua) + fik
                do jx = 1,3
                 stress(jx,ix) = stress(jx,ix) +
     &                           xno(jx,ino) * fik / volume
                enddo
               enddo
              endif
              ko = ko+2*(2*l+1)    ! Point to next group of l-/+ 1/2 blocks
             endif
             if ( ko.ge.lastkb(ka) ) exit KB_loop
            enddo KB_loop

           endif  ! listed
          enddo ! jno orbitals

C-------- Pick up contributions to H and restore Di and Vi
          do j = 1,numh(iio)
           ind = listhptr(iio)+j
           jo = listh(ind)
           H0(ind) = H0(ind) + Vi(jo)
           H0_off(ind,1) = H0_off(ind,1) + V_so(1,1,jo)
           H0_off(ind,2) = H0_off(ind,2) + V_so(2,2,jo)
           H0_off(ind,3) = H0_off(ind,3) + V_so(1,2,jo)
           H0_off(ind,4) = H0_off(ind,4) + V_so(2,1,jo)


C     Careful with this Vi()
           Vi(jo) = 0.0d0
           listed(jo) = .false.
          enddo
         endif ! Atom in UC?

        endif   ! iio .gt. 0 
       enddo     ! ino AOs
      enddo       ! atoms with KB projectors loop

      if (.not. matrix_elements_only) then
       Enl_offsiteSO = sum( dreal(E_offsiteSO(1:4)) )
      endif

C     Deallocate local memory

      call reset_neighbour_arrays( )
      call de_alloc( grSki, 'grSki', 'nlefsm_SO_off' )
      call de_alloc( Ski, 'Ski', 'nlefsm_SO_off' )
      call de_alloc( xno, 'xno', 'nlefsm_SO_off' )
      call de_alloc( iono, 'iono', 'nlefsm_SO_off' )
      call de_alloc( iano, 'iano', 'nlefsm_SO_off' )

      call de_alloc( listedall, 'listedall', 'nlefsm_SO_off' )
      call de_alloc( listed, 'listed', 'nlefsm_SO_off' )
      call de_alloc( Vi, 'Vi', 'nlefsm_SO_off' )
      deallocate( V_so )

      if (.not. matrix_elements_only) then
         call de_alloc( Di, 'Di', 'nlefsm_SO_off' )
         deallocate( Ds )
      endif

      call timer( 'nlefsm', 2 )

      end subroutine nlefsm_SO_off

c-----------------------------------------------------------------------
c
!> Evaluates:
!!   <i|V_NL|j>, where V_NL= Sum_{j,mj} |V,j,mj><V,j,mj|
c
c-----------------------------------------------------------------------
      subroutine calc_Vj_offsiteSO( l, epskb, Ski, Skj, grSki, grSkj,
     &                       V_ion, V_so, F_so, split_sr_so )

      implicit none

      ! Constructs the NL operator from the l+/- 1/2 information, which
      ! is passed in two blocks. Hence the re-dimensioning of the dummy
      ! arguments:
      !      Ski(-l:l,2) means that there are '2' blocks of 2l+1 values
      !      and so on.
      ! The two epskb values correspond to the l-1/2 and l+1/2 blocks, 
      ! respectively (all the projectors in a block have the same value)

      ! On output, V_so and V_sr (V_ion) are separated, unless split_sr_so is .false.
      ! The force contributions are not separated

      integer     , intent(in)  :: l
      real(dp)    , intent(in)  :: epskb(2)
      real(dp)    , intent(in)  :: Ski(-l:l,2), Skj(-l:l,2)
      real(dp)    , intent(in)  :: grSki(3,-l:l,2), grSkj(3,-l:l,2)
      real(dp)    , intent(out) :: V_ion
      complex(dp) , intent(out) :: F_so(3,2,2)
      complex(dp) , intent(out) :: V_so(2,2)
      logical     , intent(in)  :: split_sr_so

      integer    :: J, ij, imj, m, is, facpm
      real(dp)   :: aj, amj, al, a2l1, fac, facm,
     &              epskpm, V_iont, cp, cm

      real(dp)   :: cg(2*(2*l+1),2)
      complex(dp):: u(-l:l,-l:l)
      complex(dp):: SVi(2), SVj(2), grSVi(3,2)

      external   :: die
c-----------------------------------------------------------------------

c---- set constants and factors
      al   = dble(l)
      a2l1 = dble( 2*l+1 )

c---- load Clebsch-Gordan coefficients; cg(J,+-)
      !
      ! NOTE that this code will not work for l=0, since in this case
      ! there is a single j subspace (j=0.5)
      if ( l == 0 ) then
         call die("Code in calc_Vj_offsiteSO does not work for l=0")
      endif
      
      J = 0
      cg(:,:) = 0.0_dp

      ! ij ranges over the two j subspaces
      ! The loop could have used: 'do ik = -1,1,2' for easier reading

      do ij = 1, 2
       aj = al + (2*ij-3)*0.5d0        ! j(ij=1)=l-1/2; j(ij=2)=l+1/2

       ! This way of computing a sign is very fragile. Better below (and integer)
       !facpm= (-1.0d0)**(aj-al-0.5d0) ! +/- sign: j=l-1/2: (-1)**(-1)=-1 ;  j = l+1/2: (-1)**0 = 1
       facpm = 2*ij - 3

       do imj = 1, nint(2*aj)+1 ! Degeneracy for j: 2j+1. Safe nint as aj is always half integer
        amj = -aj + (imj-1)     ! mj value: -j, -j+1, ..., j-1, j
        J = J+1                 ! Combined (j,mj) index  (great notation!)

        cp = sqrt( (al+0.5d0+amj)/a2l1 )
        cm = sqrt( (al+0.5d0-amj)/a2l1 )
        if ( ij.eq. 1 ) then
         cg(J,1) =  cm*facpm ! <j-|up>
         cg(J,2) =  cp       ! <j-|down>
        else
         cg(J,1) =  cp*facpm ! <j+|up>
         cg(J,2) =  cm       ! <j+|down>
        endif
       enddo
      enddo

c---- Ski(M)= <l,M|i> ; Si(m)= <l,m|i> = u(m,-M)*Ski(-M) + u(m,M)*Ski(M)
      fac = 1.0d0/sqrt(2.0d0)
      u(:,:) = cmplx(0.0d0,0.0d0)
      u(0,0)= cmplx(1.0d0,0.0d0)
      do m =  1, l
       facm = fac*(-1.0d0)**m
       u(-m,+M) = cmplx(1.0d0,0.0d0)*fac
       u(-m,-M) = cmplx(0.0d0,1.0d0)*fac  ! J. Cerda
       u(+m,+M) = cmplx(1.0d0,0.0d0)*facm
       u(+m,-M) =-cmplx(0.0d0,1.0d0)*facm ! J. Cerda
      enddo

c---- Load V_so
      V_so= cmplx(0.0d0,0.0d0)
      F_so= cmplx(0.0d0,0.0d0)
      J = 0
      do ij = 1, 2
       aj = al + (2*ij-3)*0.5d0        ! j value
       do imj = 1, nint(2*aj)+1        ! Degeneracy for j
        amj = -aj + imj-1              ! mj value
        J = J+1                        ! Combined (j,mj) index

        SVi(1:2)= cmplx(0.0d0,0.0d0)
        SVj(1:2)= cmplx(0.0d0,0.0d0)
        grSVi(1:3,1:2)= cmplx(0.0d0,0.0d0)
        do is = 1, 2  ! spin loop

         ! select correct m. Safe nint as amj is always a half-integer
         if ( is.eq.1 ) then
          m = nint(amj-0.5d0)    ! up   => m=mj-1/2
         else
          m = nint(amj+0.5d0)    ! down => m=mj+1/2
         endif
     
         if ( abs(m).le.l ) then
          SVi(is)= Ski(+M,ij)*u(+m,M)
          SVj(is)= Skj(+M,ij)*u(+m,M)
          grSVi(1:3,is)= grSki(1:3,+M,ij)*u(+m,M)
          if ( m.ne.0 ) then
           SVi(is)= SVi(is) + Ski(-M,ij)*u(+m,-M)
           SVj(is)= SVj(is) + Skj(-M,ij)*u(+m,-M)
           grSVi(1:3,is)= grSVi(1:3,is) +
     &                    grSki(1:3,-M,ij)*u(+m,-M)
          endif
          SVi(is) = SVi(is) * cg(J,is)
          SVj(is) = SVj(is) * cg(J,is)
 
          grSVi(1:3,is) = grSVi(1:3,is) * cg(J,is)
         endif
        enddo ! is

        ! Note that these involve just one epskb at a time.
        
c       up-up = <i,+|V,J><V,J|j,+>
        V_so(1,1)  = V_so(1,1)  + SVi(1) * epskb(ij) * conjg(SVj(1))
        F_so(:,1,1)= F_so(:,1,1)+ grSVi(:,1) * epskb(ij) * conjg(SVj(1))

c       down-down = <i,-|V,J><V,J|j,->
        V_so(2,2)  = V_so(2,2)  + SVi(2) * epskb(ij) * conjg(SVj(2))
        F_so(:,2,2)= F_so(:,2,2)+ grSVi(:,2) * epskb(ij) * conjg(SVj(2))

c       up-down = <i,+|V,J><V,J|j,->
        V_so(1,2)  = V_so(1,2)  + SVi(1) * epskb(ij) * conjg(SVj(2))
        F_so(:,1,2)= F_so(:,1,2)+ grSVi(:,1) * epskb(ij) * conjg(SVj(2))

c       down-up= <i,-|V,J><V,J|j,+>
        V_so(2,1)  = V_so(2,1)  + SVi(2) * epskb(ij) * conjg(SVj(1))
        F_so(:,2,1)= F_so(:,2,1)+ grSVi(:,2) * epskb(ij) * conjg(SVj(1))

       enddo ! mj
      enddo ! ij

cc--- debugging
      if ( cdabs(V_so(1,2)+conjg(V_so(2,1))).gt.1.0d-4 ) then
       write(6,'(a,2f12.6)') 'V_so(1,2)=',V_so(1,2)
       write(6,'(a,2f12.6)') 'V_so(2,1)=',V_so(2,1)
       call die('calc_Vj_LS: ERROR')
      endif

      if (split_sr_so) then

      ! Attempt to get SR and SO parts of the total lj non-local contribution
      ! (Note that there is an extra contribution to SR coming from l=0, computed
      ! directly in the caller routine)

      ! This is some kind of average, following Hemstreet, of the l+/- 1/2 
      ! components. (In contrast to Hamann's more involved procedure)

      !  v_sr = ( (l+1) v_j+ + l v_j- ) / (2l+1)
      !  and the v_j+ and v_j- are "square roots" of the KB projectors      

      !!  epskpm = sqrt( epskb(1)*epskb(2) ) ! This sqrt should be guarded with an abs()
      !!  epskpm = sign(epskpm,epskb(1))   ! The value of epskpm with the sign of epskb(1) WHY?

      ! This is the only acceptable solution from symmetry arguments
      ! (and within the fragility of the approach)
      epskpm  = sqrt(abs(epskb(1)*epskb(2)))
      if (epskb(1)*epskb(2) > 0) then
         ! same sign
         epskpm  = epskpm * sign(1.0_dp,epskb(1)) 
      else
         ! different sign: this possibility was not taken into account in previous versions,
         ! and probably implies that the Hemstreet ansatz is ill-defined
         ! We now flag this as a warning (logic in caller) instead of accepting this heuristic blindly.
         epskpm  = - epskpm
      endif
      
      V_ion = 0.0d0
      do M = -l, l
       V_iont = ( l**2     * Ski(M,1)*epskb(1)*Skj(M,1)
     &          + (l+1)**2 * Ski(M,2)*epskb(2)*Skj(M,2)
     &          + l*(l+1)  * Ski(M,1)*epskpm  *Skj(M,2)
     &          + l*(l+1)  * Ski(M,2)*epskpm  *Skj(M,1) )/(a2l1**2)
       V_ion = V_ion + V_iont
      enddo
   
      !---- substract out V_ion from V_so
      V_so(1,1) = V_so(1,1) - cmplx(1.0d0,0.0d0)*V_ion
      V_so(2,2) = V_so(2,2) - cmplx(1.0d0,0.0d0)*V_ion

      else
         !---- Keep all the NL contribution in V_so
         V_ion = 0.0_dp
      endif

      end subroutine calc_Vj_offsiteSO

      end module m_nlefsm
