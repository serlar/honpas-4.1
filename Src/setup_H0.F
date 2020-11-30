! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_setup_H0
      private
      public :: setup_H0
      CONTAINS

      subroutine setup_H0(G2max)
      
C     Computes non-self-consistent part of the Hamiltonian
C     and initializes data structures on the grid.
      
      USE siesta_options,  only: g2cut, split_sr_so
      USE siesta_options,  only: hfx_wanted
      use sparse_matrices, only: H_kin_1D, H_vkb_1D
      use sparse_matrices, only: H_so_on_2D, H_so_off_2D
      use sparse_matrices, only: Dscf

      use m_nlefsm,        only: nlefsm_SO_off
      use m_spin,          only: spin

      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, indxuo, datm, 
     &                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     &                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_nlefsm,     only: nlefsm
      use m_kinefsm,    only: kinefsm
      use m_naefs,      only: naefs
      use m_dnaefs,     only: dnaefs
      use m_dhscf,      only: dhscf_init
      use m_energies,   only: Eions, Ena, DEna, Emm, Emeta, Eso
      use m_ntm
      use m_spin,       only: spin
      use spinorbit,    only: spinorb
      use alloc, only: re_alloc, de_alloc
      use class_dSpData1D, only: val
      use class_dSpData2D, only: val
      use class_zSpData2D, only: val
      use atm_types, only: l_max
      use nao2gto_data, only: hfx_libint, hfx_options, hfx_system,
     &                        nao2gto_libint_init, nao2gto_system_init,
     &                        nao2gto_system_update_cell
      use nao2gto_prescreen, only: init_prescreen_eri

#ifdef MPI
      use m_mpi_utils, only: globalize_sum
      use parallel,    only: Node
#endif

      implicit none
      real(dp), intent(inout) :: g2max
      
      real(dp) :: dummy_stress(3,3), dummy_fa(1,1), dummy_dm(1,1)
      real(dp) :: dummy_E
      integer  :: ia, is
      integer  :: ic, jc, ncells
      real(dp) :: max_eri

      real(dp)    :: dummy_Eso
      integer     :: ispin, i, j 
      complex(dp) :: Dc
#ifdef MPI
      real(dp) :: buffer1
#endif

      real(dp), pointer :: H_val(:), H_so_on(:,:)
      complex(dp), pointer :: H_so_off(:,:)

      logical :: ill_defined_sr_so_split
      
#ifdef DEBUG
      call write_debug( '    PRE setup_H0' )
#endif

!----------------------------------------------------------------------BEGIN
      call timer('Setup_H0',1)

C     Self-energy of isolated ions
      Eions = 0.0_dp
      do ia = 1,na_u
        is = isa(ia)
        Eions = Eions + uion(is)
      enddo

!     In these routines, add a flag to tell them NOT to compute
!     forces and stresses in this first pass, only energies.

!     Neutral-atom: energy

      call naefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &           isa, Ena, dummy_fa, dummy_stress,
     &           forces_and_stress=.false.)
      call dnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, DEna, dummy_fa, dummy_stress,
     &            forces_and_stress=.false.) 
      Ena = Ena + DEna
 
C     Metadynamics energy
      if (lMetaForce) then
        call meta(xa,na_u,ucell,Emeta,dummy_fa,dummy_stress,
     $           .false.,.false.)
      endif

C     Add on force field contribution to energy
      call twobody( na_u,xa,isa,ucell,Emm,
     &              ifa=0,fa=dummy_fa,istr=0,stress=dummy_stress)

!
!     Now we compute matrix elements of the Kinetic and Non-local
!     parts of H

!     Kinetic: matrix elements only
      H_val => val(H_kin_1D)
!$OMP parallel workshare default(shared)
      H_val(:) = 0.0_dp
!$OMP end parallel workshare
      call kinefsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     &             maxnh, maxnh, lasto, iphorb, isa, 
     &             numh, listhptr, listh, numh, listhptr, listh, 
     &             1,
     &             dummy_dm, dummy_E, dummy_fa, dummy_stress,
     &             H_val,
     &             matrix_elements_only=.true.)

!     Non-local-pseudop:  matrix elements only
      H_val => val(H_vkb_1D)
!$OMP parallel workshare default(shared)
      H_val(:) = 0.0_dp
!$OMP end parallel workshare

      Eso = 0.0d0

      if ( .not.spin%SO_offsite ) then
       call nlefsm(scell, na_u, na_s, isa, xa, indxua, 
     &             maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     &             numh, listhptr, listh, numh, listhptr, listh, 
     &             1,
     &             dummy_dm, dummy_E, dummy_fa, dummy_stress,
     &             H_val,
     &             matrix_elements_only=.true.) 
      else
        H_so_off => val(H_so_off_2D)
        H_so_off = dcmplx(0._dp, 0._dp)
        call nlefsm_SO_off(scell, na_u, na_s, isa, xa, indxua,
     &                 maxnh, maxnh, lasto, lastkb, iphorb, iphKB,
     &                 numh, listhptr, listh, numh, listhptr, listh,
     &                 spin%Grid,
     &                 dummy_E, dummy_Eso, dummy_fa,
     &                 dummy_stress, H_val, H_so_off,
     &                 matrix_elements_only=.true.,
     &                 ill_defined_sr_so_split=ill_defined_sr_so_split)

        if (ill_defined_sr_so_split) then
           if (split_sr_so) then
              call message("WARNING",
     $        "The Enl-Eso split in energies might not be accurate")
           else
              call message("INFO",
     $        "The Enl-Eso split in energies " //
     $        "would not have been accurate")
           endif
        endif

!
!  Dc IS NOT the dense matrix, it is just a complex number 
! (per each io index) used as an artifact to multiply the 
! elements of the H_SO times the corresponding elements of 
! DM in a such way that the result gives Re{Tr[H_SO*DM]}.
!

        do i = 1, maxnh

!-------- Eso(u,u)
          Dc = cmplx(Dscf(i,1),Dscf(i,5), dp)
          Eso = Eso + real( H_so_off(i,1)*Dc, dp)
!-------- Eso(d,d)
          Dc = cmplx(Dscf(i,2),Dscf(i,6),dp)
          Eso = Eso + real( H_so_off(i,2)*Dc, dp)
!-------- Eso(u,d)
          Dc = cmplx(Dscf(i,3),Dscf(i,4), dp)
          Eso = Eso + real( H_so_off(i,4)*Dc, dp)
!-------- Eso(d,u)
          Dc = cmplx(Dscf(i,7),-Dscf(i,8), dp)
          Eso = Eso + real( H_so_off(i,3)*Dc, dp)

        enddo

#ifdef MPI
! Global reduction of Eso 
      call globalize_sum(Eso,buffer1)
      Eso = buffer1
#endif
      endif
! ..................


! If in the future the spin-orbit routine is able to compute
! forces and stresses, then "last" will be needed. If we are not
! computing forces and stresses, calling it in the first iteration
! should be enough
!
      if ( spin%SO_onsite ) then
         H_so_on => val(H_so_on_2D)
!$OMP parallel workshare default(shared)
         H_so_on(:,:) = 0._dp
!$OMP end parallel workshare
         call spinorb(no_u,no_l,iaorb,iphorb,isa,indxuo,
     &        maxnh,numh,listhptr,listh,Dscf,H_so_on,Eso)
      end if

C     This will take care of possible changes to the mesh and atomic-related
C     mesh structures for geometry changes
      g2max = g2cut
      call dhscf_init( spin%Grid, no_s, iaorb, iphorb,
     &                 no_l, no_u, na_u, na_s,
     &                 isa, xa, indxua, ucell,
     &                 mscell, G2max, ntm,
     &                 maxnh, numh, listhptr, listh, datm,
     &                 dummy_fa, dummy_stress)

!     ----- Begin HFX part -----
      if ( hfx_wanted ) then
!       Prepare atom-dependent data for Hartree-Fock exchange
!       FIXME: should not be called if atoms and cells are identical
!              to previous step
        call nao2gto_system_update_cell(hfx_system, nsc, ucell)
        call init_prescreen_eri(hfx_libint, hfx_options, hfx_system)
      endif ! hfx_wanted
!     ----- End HFX part -----

      call timer('Setup_H0',2)

#ifdef DEBUG
      call write_debug( '    POS setup_H0' )
#endif

!---------------------------------------------------------------------- END
      END subroutine setup_H0
      END module m_setup_H0
