! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_compute_energies
  implicit none
  public :: compute_energies
CONTAINS
  subroutine  compute_energies(iscf)

    ! This routine computes the Harris energy from E_KS(DM_in):
    !
    !    E_Harris = E_KS(DM_in) + Tr[H_in*(DM_out-DM_in)]
    !
    ! and possibly E_KS(DM_out) as
    !
    !    E_KS(DM_out) = Tr[H_in*DM_out] + E_HXC(DM_out)
    !
    ! Note that E_KS(DM_in), as computed in setup_hamiltonian, is not
    ! variational if mixing the DM, as the kinetic energy term is
    ! using a DM (DM_in) which is not derived from wave-functions. It
    ! is also not consistent if we are using the charge density as
    ! primary variable for mixing. In this case the Enl and Ekin parts
    ! are computed with the previous step's DM, whereas the "scf" part
    ! comes from the (mixed) rho.
    !
    ! E_KS(DM_out) computed in setup_hamiltonian with DM_out is
    ! variational. Rather than calling again setup_hamiltonian, it is
    ! enough to compute the "scf" part of the energy by calling dhscf
    ! with DM_out.

    ! The routine checks the mixing type and proceeds accordingly

      use precision,       only: dp
      use fdf,             only: fdf_get
      use siesta_options,  only: g2cut, Temp
      use siesta_options,  only: mix_charge, mixH
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H
      use sparse_matrices, only: Dscf, Dold
      use m_dhscf,         only: dhscf
      use m_energies
      use atomlist,        only: no_u, iaorb, iphkb, qtot, indxuo, datm,   &
                                 lastkb, no_s, rmaxv, indxua, iphorb, lasto, &
                                 rmaxo, no_l
      use m_ntm,           only: ntm

      use m_spin,          only: spin
      use parallel,        only: IONode 

      use m_dipol,         only: dipol
      use siesta_geom,     only: na_u, na_s, xa, isa
      use m_rhog,          only: rhog
#ifdef MPI
      use m_mpi_utils,     only: globalize_sum
#endif

      integer, intent(in)   :: iscf

      integer               :: ispin, io
#ifdef MPI
      real(dp) :: buffer1
#endif

      logical  :: mixDM


      mixDM = (.not. (mixH .or. mix_charge))

!     Compute the band-structure energy

      call compute_EBS()

      ! These energies were calculated in the latest call to
      ! setup_hamiltonian, using as ingredient D_in 

      ! Ecorrec comes from O(N)...
      ! DUext is the energy of the charge from DM_in in a field.
      ! Emad, Emm, Emeta are extra terms that are added for
      ! consistency of the total energy.

      ! E0 = Ena + Ekin + Enl + Eso - Eions

      call update_DEna()
      call update_Etot()


! Harris energy
! It does not make sense if mixing the charge density, as there is no DM_in
! If mixing the Hamiltonian its usefulness is questionable also.

      if (mix_charge) then   ! possibly add mixH here
         EHarrs = 0.0_dp
      else
         call compute_DEharr()
         Eharrs = Etot + DEharr
      endif

! Possible correction to Etot if mixing the DM. This is purely
! cosmetic, to show uniformly converged values during the scf
! cycle. The final energy will be completely correct if the DM is not
! mixed after scf convergence.

! The correction is mandatory if mixing the charge density. In this
! case the call to dhscf is needed to generate rho(G)_out

! If mixing H the KS energy is already variational, as it is
! computed with DM_out

      if (mixDM) then
         if (fdf_get("SCF.Want.Variational.EKS",.false.)) then
            call compute_correct_EKS()
         endif
      else if (mixH) then
         ! not needed
      else if (mix_charge) then
         call compute_correct_EKS()
      endif

      call update_FreeE(Temp)

 CONTAINS

   subroutine compute_EBS()

     real(dp)    :: Ebs_SO(4)
     complex(dp) :: Ebs_Daux(2,2), Ebs_Haux(2,2)

      integer :: i, j, ind

      Ebs = 0.0_dp

! Modifed for Off-Site Spin-orbit coupling by R. Cuadrado, Feb. 2018
!
!*****************************************************************************
!  Note about Ebs and E_Harris calculation for the Spin-Orbit:
!*****************************************************************************
!
!  E_bs and E_Harris are calculated by means of the following 
! complex matrix multiplication: 
!
!                 E_bs = Re { Tr[ H * DM ] } 
!             E_Harris = Re { Tr[ H * (DM-DM_old) ] }  
!
!  In the following: DM/H(1,1) --> up/up     <--> uu
!                    DM/H(2,2) --> down/down <--> dd
!                    DM/H(1,2) --> up/down   <--> ud
!                    DM/H(2,1) --> down/up   <--> du
!
!  Using DM/H components, E_bs would be sum(E_bs(1:4)), where
!
!   E_bs(1)=Re{sum_ij(H_ij(1,1)*D_ji(1,1))}=Re{sum_ij(H_ij^uu*(DM_ij^uu)^*)}
!   E_bs(2)=Re{sum_ij(H_ij(2,2)*D_ji(2,2))}=Re{sum_ij(H_ij^dd*(DM_ij^dd)^*)}
!   E_bs(3)=Re{sum_ij(H_ij(1,2)*D_ji(2,1))}=Re{sum_ij(H_ij^ud*(DM_ij^ud)^*)}
!   E_bs(4)=Re{sum_ij(H_ij(2,1)*D_ji(1,2))}=Re{sum_ij(H_ij^du*(DM_ij^du)^*)}
!
!         since, due to overall hermiticity,  DM_ij^ab = (DM_ji^ba)^*
!
!   The trace operation is then an extended dot product over the "ij"
!   sparse index, which can also be conveniently done in parallel, as
!   each processor handles the same indexes in H and the DM. Only a
!   global reduction is needed at the end.
      
!  Same comments are valid for the E_Harris calculation.
!
!*****************************************************************************
!
      if ( spin%SO ) then

        Ebs_SO = 0.0_dp  
        Ebs_Daux = dcmplx(0.0_dp, 0.0_dp)
        Ebs_Haux = dcmplx(0.0_dp, 0.0_dp)

        do io = 1, maxnh

          Ebs_Haux(1,1) = dcmplx(H(io,1),H(io,5))  
          Ebs_Haux(2,2) = dcmplx(H(io,2),H(io,6))  
          Ebs_Haux(1,2) = dcmplx(H(io,3),-H(io,4)) 
          Ebs_Haux(2,1) = dcmplx(H(io,7),H(io,8)) 

          Ebs_Daux(1,1) = dcmplx(Dscf(io,1),Dscf(io,5)) 
          Ebs_Daux(2,2) = dcmplx(Dscf(io,2),Dscf(io,6)) 
          Ebs_Daux(1,2) = dcmplx(Dscf(io,3),-Dscf(io,4))
          Ebs_Daux(2,1) = dcmplx(Dscf(io,7),Dscf(io,8)) 


          Ebs_SO(1) = Ebs_SO(1) + real( Ebs_Haux(1,1)*dconjg(Ebs_Daux(1,1)) )
          Ebs_SO(2) = Ebs_SO(2) + real( Ebs_Haux(2,2)*dconjg(Ebs_Daux(2,2)) )
          Ebs_SO(3) = Ebs_SO(3) + real( Ebs_Haux(1,2)*dconjg(Ebs_Daux(1,2)) )
          Ebs_SO(4) = Ebs_SO(4) + real( Ebs_Haux(2,1)*dconjg(Ebs_Daux(2,1)) )

        enddo

        Ebs = sum ( Ebs_SO )

      else if ( spin%NCol ) then
        do io = 1,maxnh
          Ebs    = Ebs    + H(io,1) * ( Dscf(io,1)  )   &
                          + H(io,2) * ( Dscf(io,2)  )   &
                 + 2.0_dp * H(io,3) * ( Dscf(io,3)  )   &
                 + 2.0_dp * H(io,4) * ( Dscf(io,4)  )
        enddo
      else if ( spin%Col )  then
        do io = 1,maxnh
          Ebs    = Ebs    + H(io,1) * Dscf(io,1)  &
                          + H(io,2) * Dscf(io,2)
        enddo
      else if ( spin%none ) then
        do io = 1,maxnh
          Ebs    = Ebs + H(io,1) *  Dscf(io,1)
        enddo
      endif

#ifdef MPI
!     Global reduction 
      call globalize_sum(Ebs,buffer1)
      Ebs = buffer1
#endif
    end subroutine compute_EBS

    subroutine compute_DEharr()

      real(dp)    :: DEharr_SO(4)
      complex(dp) :: DEharr_Daux(2,2), DEharr_Haux(2,2), DEharr_Daux_old(2,2)

      DEharr = 0.0_dp

      if ( spin%SO ) then

        DEharr_SO = 0.0_dp 
        DEharr_Daux     = dcmplx(0.0_dp, 0.0_dp)
        DEharr_Daux_old = dcmplx(0.0_dp, 0.0_dp)
        DEharr_Haux     = dcmplx(0.0_dp, 0.0_dp)

        do io = 1, maxnh

          DEharr_Haux(1,1) = dcmplx( H(io,1),H(io,5) )
          DEharr_Haux(2,2) = dcmplx( H(io,2),H(io,6) )
          DEharr_Haux(1,2) = dcmplx( H(io,3),-H(io,4) )
          DEharr_Haux(2,1) = dcmplx( H(io,7),H(io,8) )

          DEharr_Daux(1,1) = dcmplx( Dscf(io,1),Dscf(io,5) )
          DEharr_Daux(2,2) = dcmplx( Dscf(io,2),Dscf(io,6) )
          DEharr_Daux(1,2) = dcmplx( Dscf(io,3),-Dscf(io,4) )
          DEharr_Daux(2,1) = dcmplx( Dscf(io,7),Dscf(io,8) )

          DEharr_Daux_old(1,1) = dcmplx( Dold(io,1),Dold(io,5) )
          DEharr_Daux_old(2,2) = dcmplx( Dold(io,2),Dold(io,6) ) 
          DEharr_Daux_old(1,2) = dcmplx( Dold(io,3),-Dold(io,4) ) 
          DEharr_Daux_old(2,1) = dcmplx( Dold(io,7),Dold(io,8) )


          DEharr_SO(1) = DEharr_SO(1) &
             + real( DEharr_Haux(1,1)*dconjg(DEharr_Daux(1,1)-DEharr_Daux_old(1,1)) ) 

          DEharr_SO(2) = DEharr_SO(2) &
             + real( DEharr_Haux(2,2)*dconjg(DEharr_Daux(2,2)-DEharr_Daux_old(2,2)) )

          DEharr_SO(3) = DEharr_SO(3) &
             + real( DEharr_Haux(1,2)*dconjg(DEharr_Daux(1,2)-DEharr_Daux_old(1,2)) )

          DEharr_SO(4) = DEharr_SO(4) &
             + real( DEharr_Haux(2,1)*dconjg(DEharr_Daux(2,1)-DEharr_Daux_old(2,1)) )

         enddo

         DEharr = sum ( DEharr_SO )

      else if ( spin%NCol ) then
        do io = 1,maxnh
          DEharr = DEharr + H(io,1) * ( Dscf(io,1) - Dold(io,1) )  &
                          + H(io,2) * ( Dscf(io,2) - Dold(io,2) )  &
                 + 2.0_dp * H(io,3) * ( Dscf(io,3) - Dold(io,3) )  &
                 + 2.0_dp * H(io,4) * ( Dscf(io,4) - Dold(io,4) )
        enddo
      elseif ( spin%Col )  then
        do io = 1,maxnh
          DEharr = DEharr + H(io,1) * ( Dscf(io,1) - Dold(io,1) )  &
                          + H(io,2) * ( Dscf(io,2) - Dold(io,2) )
        enddo
      elseif ( spin%none ) then
        do io = 1,maxnh
          DEharr = DEharr + H(io,1) * ( Dscf(io,1) - Dold(io,1) )
        enddo
      endif

#ifdef MPI
      !     Global reduction of DEharr
      call globalize_sum( DEharr, buffer1 )
      DEharr = buffer1
#endif
    end subroutine compute_DEharr
      
    subroutine compute_correct_EKS()

      use files, only : filesOut_t    ! derived type for output file names
      use class_dSpData1D, only : val
      use class_dSpData2D, only : val
      use class_zSpData2D, only : val
      use sparse_matrices, only: H_kin_1D, H_vkb_1D
      use sparse_matrices, only: H_so_on_2D, H_so_off_2D


      type(filesOut_t)  :: filesOut  ! blank output file names
      real(dp), pointer :: H_vkb(:), H_kin(:), H_so_on(:,:)
      complex(dp), pointer :: H_so_off(:,:)


      complex(dp) :: Hc, Dc
      real(dp)    :: dummy_stress(3,3), dummy_fa(1,1)
      real(dp)    :: dummy_E, g2max, dummy_H(1,1)
      integer     :: ihmat, ifa, istr
      
      ! Compute E_KS(DM_out)

      g2max = g2cut
      ifa  = 0
      istr = 0
      ihmat = 0

      ! Pass DM_out to compute E_HXC(out)

      ! Remove unwanted arguments...

      call dhscf( spin%Grid, no_s, iaorb, iphorb, no_l,      &
                  no_u, na_u, na_s, isa, xa, indxua,                        &
                  ntm, ifa, istr, ihmat, filesOut,                          &
                  maxnh, numh, listhptr, listh, Dscf, Datm,                 &
                  maxnh, dummy_H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, &
                  Exc, Dxc, dipol, dummy_stress, dummy_fa, dummy_stress)
      ! (when mix_charge is true, rhog will contain the output rho(G))

      call update_DEna()

!     Compute Tr[H_0*DM_out] = Ekin + Enl + Eso with DM_out

      H_vkb => val(H_vkb_1D)
      H_kin => val(H_kin_1D)
      
      Ekin = 0.0_dp
      Enl  = 0.0_dp
      do ispin = 1, spin%spinor
         do io = 1,maxnh
            Ekin = Ekin + H_kin(io) * Dscf(io,ispin)
            Enl  = Enl  + H_vkb(io) * Dscf(io,ispin)
         enddo
      enddo

#ifdef MPI
      ! Global reduction of Ekin, Enl
      call globalize_sum( Ekin, buffer1 )
      Ekin = buffer1
      call globalize_sum( Enl, buffer1 )
      Enl = buffer1
#endif

      Eso = 0._dp

      if ( spin%SO_offsite ) then
         H_so_off => val(H_so_off_2D)

         ! The computation of the trace is different here, as H_so_off has
         ! a different structure from H and the DM.
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
         ! Sadly some compilers (g95), does
         ! not allow bounds for pointer assignments :(
         H_so_on => val(H_so_on_2D)
         do io = 1,maxnh
            Eso = Eso + H_so_on(io,1)*Dscf(io,7) + H_so_on(io,2)*Dscf(io,8) &
                 + H_so_on(io,5)*Dscf(io,3) + H_so_on(io,6)*Dscf(io,4) &
                 - H_so_on(io,3)*Dscf(io,5) - H_so_on(io,4)*Dscf(io,6)
          end do
          
      end if

#ifdef MPI
      if (spin%SO) then
         ! Global reduction of Eso
         call globalize_sum( Eso, buffer1 )
         Eso = buffer1
      endif
#endif
      
      ! E0 = Ena + Ekin + Enl + Eso - Eions

      ! Clarify: Ecorrec (from O(N))
      call update_Etot()
      
    end subroutine compute_correct_EKS

  end subroutine compute_energies

end module m_compute_energies
