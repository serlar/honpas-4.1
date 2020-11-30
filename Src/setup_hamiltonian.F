! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_setup_hamiltonian
      private
      public :: setup_hamiltonian
      CONTAINS

      subroutine setup_hamiltonian( iscf )

      USE siesta_options
      use sparse_matrices, only: H_kin_1D, H_vkb_1D
      use sparse_matrices, only: H_ldau_2D
      use sparse_matrices, only: H_so_on_2D, H_so_off_2D
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H, S, Hold
      use sparse_matrices, only: Dscf, Escf, xijo
      use class_dSpData1D,  only: val
      use class_dSpData2D,  only: val
      use class_zSpData2D,  only: val

      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     .                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     .                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use ldau_specs,   only: switch_ldau     ! This variable determines whether
                                              !   the subroutine to compute the
                                              !   Hubbard terms should be called
                                              !   or not
      use m_ldau,       only: hubbard_term    ! Subroutine that compute the
                                              !   Hubbard terms
      use m_dhscf,      only: dhscf
      use m_stress
      use m_energies
      use parallel, only: Node
      use m_steps, only: istp
      use m_ntm

      use m_spin,          only: spin 

      use m_dipol
      use alloc, only: re_alloc, de_alloc
      use m_hsx, only: write_hsx
      use sys, only: die, bye
      use m_partial_charges, only: want_partial_charges
      use files, only : filesOut_t    ! derived type for output file names
      use m_rhog, only: rhog_in, rhog
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
      use nao2gto_data, only: hfx_libint, hfx_options, hfx_system
      use nao2gto_hfx,  only: setup_hfx
      use nao2gto_hfx_force,  only: setup_hfx_force

      implicit none
      integer, intent(in) :: iscf
      real(dp)            :: stressl(3,3)
      real(dp), pointer   :: fal(:,:)   ! Local-node part of atomic F
#ifdef MPI
      real(dp)            :: buffer1
#endif
      integer             :: io, is, ispin
      integer             :: ifa     ! Calc. forces?      0=>no, 1=>yes
      integer             :: istr    ! Calc. stress?      0=>no, 1=>yes
      integer             :: ihmat   ! Calc. hamiltonian? 0=>no, 1=>yes
      real(dp)            :: g2max
      type(filesOut_t)    :: filesOut  ! blank output file names
      logical             :: use_rhog_in

      real(dp), pointer   :: H_vkb(:), H_kin(:), H_ldau(:,:)
      real(dp), pointer   :: H_so_on(:,:)
      complex(dp), pointer:: H_so_off(:,:)

      complex(dp):: Dc
      integer :: ind, i, j

!------------------------------------------------------------------------- BEGIN

      call timer('setup_H',1)

      ! Nullify pointers
      nullify(fal)

!$OMP parallel default(shared), private(ispin,io)

!     Save present H matrix
!$OMP do collapse(2)
      do ispin = 1, spin%H
         do io = 1,maxnh
            Hold(io,ispin) = H(io,ispin)
         end do
      end do
!$OMP end do

!$OMP single
      H_kin => val(H_kin_1D)
      H_vkb => val(H_vkb_1D)
      
      if ( spin%SO_onsite ) then
        ! Sadly some compilers (g95), does
        ! not allow bounds for pointer assignments :(
        H_so_on => val(H_so_on_2D)
        
      else if ( spin%SO_offsite ) then
        H_so_off => val(H_so_off_2D)
        
      end if
!$OMP end single ! keep wait

      ! Initialize diagonal Hamiltonian
      do ispin = 1, spin%spinor
!$OMP do
        do io = 1,maxnh
          H(io,ispin) = H_kin(io) + H_vkb(io)
        end do
!$OMP end do nowait
      end do

      if ( spin%SO_onsite ) then
!$OMP do collapse(2)
        do ispin = 3 , spin%H
          do io = 1,maxnh
            H(io,ispin) = H_so_on(io,ispin-2)
          end do
        end do
!$OMP end do nowait
         
      else

!$OMP do collapse(2)
        do ispin = 3 , spin%H
          do io = 1,maxnh
            H(io,ispin) = 0._dp
          end do
        end do
!$OMP end do nowait
         
      end if
         
! ..................

! Non-SCF part of total energy .......................................
! Note that these will be "impure" for a mixed Dscf

! If mixing the charge, Dscf is the previous step's DM_out. Since
! the "scf" components of the energy are computed with the (mixed)
! charge, this introduces an inconsistency. In this case the energies
! coming out of this routine need to be corrected.
! 
!$OMP single
      Ekin = 0.0_dp
      Enl = 0.0_dp
      Eso = 0.0_dp
!$OMP end single ! keep wait
      
!$OMP do collapse(2), reduction(+:Ekin,Enl)
      do ispin = 1, spin%spinor
        do io = 1,maxnh
          Ekin = Ekin + H_kin(io) * Dscf(io,ispin)
          Enl  = Enl  + H_vkb(io) * Dscf(io,ispin)
        end do
      end do
!$OMP end do nowait

!
!  Dc IS NOT the dense matrix, it is just a complex number 
! (per each io index) used as an artifact to multiply the 
! elements of the H_SO times the corresponding elements of 
! DM in a such way that the result gives Re{Tr[H_SO*DM]}.
!

      if ( spin%SO_offsite ) then
        do io = 1, maxnh

!-------- Eso(u,u)
          Dc = cmplx(Dscf(io,1),Dscf(io,5), dp)
          Eso = Eso + real( H_so_off(io,1)*Dc, dp)
!-------- Eso(d,d)
          Dc = cmplx(Dscf(io,2),Dscf(io,6), dp)
          Eso = Eso + real( H_so_off(io,2)*Dc, dp)
!-------- Eso(u,d)
          Dc = cmplx(Dscf(io,3),Dscf(io,4), dp)
          Eso = Eso + real( H_so_off(io,4)*Dc, dp)
!-------- Eso(d,u)
          Dc = cmplx(Dscf(io,7),-Dscf(io,8), dp)
          Eso = Eso + real( H_so_off(io,3)*Dc, dp)

        end do

      else if ( spin%SO_onsite ) then
        
!$OMP do reduction(+:Eso)
         do io = 1, maxnh
           Eso = Eso + H_so_on(io,1)*Dscf(io,7) +
     &         H_so_on(io,2)*Dscf(io,8)+ H_so_on(io,5)*Dscf(io,3) +
     &         H_so_on(io,6)*Dscf(io,4)- H_so_on(io,3)*Dscf(io,5) -
     &         H_so_on(io,4)*Dscf(io,6)
         end do
!$OMP end do nowait

       end if

!$OMP end parallel

#ifdef MPI
      ! Global reduction of Ekin, Enl
      call globalize_sum(Ekin,buffer1)
      Ekin = buffer1
      call globalize_sum(Enl,buffer1)
      Enl = buffer1
      if ( spin%SO ) then
         ! Global reduction of Eso 
         call globalize_sum(Eso,buffer1)
         Eso = buffer1
      end if
#endif

!     Non-SCF part of total energy
      call update_E0()

! Hubbard term for LDA+U: energy, forces, stress and matrix elements ....
      if( switch_ldau ) then
        if ( spin%NCol ) then
         call die('LDA+U cannot be used with non-collinear spin.')
        end if
        if ( spin%SO ) then
         call die('LDA+U cannot be used with spin-orbit coupling.')
        end if

        call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'setup_hamiltonian' )

        H_ldau => val(H_ldau_2D)
        call hubbard_term(scell, na_u, na_s, isa, xa, indxua,
     .                    maxnh, maxnh, lasto, iphorb, no_u, no_l, 
     .                    numh, listhptr, listh, numh, listhptr, listh,
     .                    spin%spinor, Dscf, Eldau, DEldau, H_ldau, 
     .                    fal, stressl, H, iscf,
     .                    matrix_elements_only=.true.)

#ifdef MPI
        ! Global reduction of energy terms
        call globalize_sum(Eldau,buffer1)
        Eldau = buffer1
        ! DEldau should not be globalized
        ! as it is based on globalized occupations
#endif
        Eldau = Eldau + DEldau

        call de_alloc( fal, 'fal', 'setup_hamiltonian' ) 
      endif
! ..................


! Add SCF contribution to energy and matrix elements ..................
      g2max = g2cut

      call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'setup_hamiltonian' )

      ifa  = 0
      istr = 0
      ihmat = 1
      if ((hirshpop .or. voropop)
     $     .and. partial_charges_at_every_scf_step) then
         want_partial_charges = .true.
      endif
      use_rhog_in =  (mix_charge .and. iscf > 1)
  
      call dhscf( spin%Grid, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua, 
     .            ntm, ifa, istr, ihmat, filesOut,
     .            maxnh, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, stress, fal, stressl,
     .            use_rhog_in)

      if ( spin%SO_offsite ) then

! H(:, [5, 6]) are not updated in dhscf, see vmat for details.
        
!------- H(u,u)
         H(:,1) = H(:,1) + real(H_so_off(:,1), dp)
         H(:,5) =         dimag(H_so_off(:,1))
!------- H(d,d)
         H(:,2) = H(:,2) + real(H_so_off(:,2), dp)
         H(:,6) =         dimag(H_so_off(:,2))
!------- H(u,d)
         H(:,3) = H(:,3) + real(H_so_off(:,3), dp)
         H(:,4) = H(:,4) +dimag(H_so_off(:,3))
!------- H(d,u)
         H(:,7) = H(:,7) + real(H_so_off(:,4), dp)
         H(:,8) = H(:,8) -dimag(H_so_off(:,4))

      endif

!-----------------------------------------------------------------------
! add Hartree-Fock exchange potential to Hamiltonian (nuo, norb)
! xmqin, January 2014
!----------------------------------------------------------------------
      call timer( 'DHFX', 1 )
!      if(node.eq.0) write(6,*) 'Hybrid DFT calculations'
      if ( hfx_wanted ) then
        write(*,*) "CALL SETUP_HFX"
        call setup_hfx(hfx_libint, hfx_options, hfx_system, H)
      endif
      call timer( 'DHFX', 2 )
!-----------------------------------------------------------------------
! add Hartree-Fock exchange gradients to force (1:3, nua)
! xmqin, January 2014
!----------------------------------------------------------------------
!      call timer( 'DFHFX', 1 )
!      if ( hfx_wanted ) then
!        call setup_hfx_force(hfx_libint, hfx_options, hfx_system, fal)
!      endif
!      call timer( 'DFHFX', 2 )
!----------------------------------------------------------------------

      ! This statement will apply to iscf = 1, for example, when
      ! we do not use rhog_in. Rhog here is always the charge used to
      ! build H, that is, rhog_in.
      if (mix_charge) rhog_in = rhog

      want_partial_charges = .false.
      call de_alloc( fal, 'fal', 'setup_hamiltonian' ) 
!  It is wasteful to write over and over H and S, as there are
!  no different files.
! Save Hamiltonian and overlap matrices ............................
! Only in HSX format now.  Use Util/HSX/hsx2hs to generate an HS file

      if (savehs .or. write_coop) then
         call write_hsx( no_u, no_s, spin%H, indxuo,
     &        maxnh, numh, listhptr, listh, H, S, qtot,
     &        temp, xijo)
      endif

      call timer('setup_H',2)
#ifdef SIESTA__PEXSI
      if (node==0) call memory_snapshot("after setup_H")
#endif

      if ( h_setup_only ) then
         call timer( 'all', 2 )  ! New call to close the tree
         call timer( 'all', 3 )
         call bye("H-Setup-Only requested")
         STOP
      endif

!------------------------------------------------------------------------- END
      END subroutine setup_hamiltonian
      END module m_setup_hamiltonian
