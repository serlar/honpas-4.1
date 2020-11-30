! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

module m_local_DOS

  implicit none
  private

  public :: init_local_DOS
  public :: local_DOS

contains

  subroutine init_local_DOS( ucell )

    use precision, only: dp
    use siesta_options
    use fdf,         only: fdf_isblock
    use kpoint_ldos_m, only: setup_kpoint_ldos
    use parallel,    only: IOnode

    real(dp), intent(in) :: ucell(3,3)

    !-------------------------------------------------------------------------BEGIN
    !     Compute the projected density of states
    do_ldos = fdf_isblock('LocalDensityOfStates')
    if ( .not. do_ldos ) return

    if ( isolve /= SOLVE_DIAGON ) then
      if (.not.((isolve == SOLVE_MINIM).and. minim_calc_eigenvalues)) then
        if (IONode) then
          write(*,*) 'siesta: ERROR: LDOS implemented only with diagon'
        end if
        do_ldos = .false.
      end if
    end if

    if ( .not. do_ldos ) return

    call setup_kpoint_ldos( ucell )
      
  end subroutine init_local_DOS

  subroutine local_DOS( )

    use units, only: eV
    use alloc, only: re_alloc
    use m_energies
    use sparse_matrices
    use siesta_options
    use siesta_geom
    use atomlist,       only: indxuo, indxua           
    use atomlist,       only: qtot, qtots, no_u, no_l
    use atomlist,       only: iphorb                   
    use atomlist,       only: datm, no_s, iaorb        
    use fdf
    use sys,            only: die                      
    use files,          only: slabel     ! system label
    use files,          only: filesOut_t ! derived type for output file names
    use kpoint_scf_m, only: kpoint_scf
    use kpoint_ldos_m, only: kpoint_ldos, gamma_ldos
    use parallel,       only: IOnode                   
    use files,          only : label_length            
    use m_ntm
    use m_forces,       only: fa
    use m_energies,     only: Ef, Efs
    use m_eo
    use m_spin,         only: nspin
    use m_spin,         only: spinor_dim 
    use m_diagon,       only: diagon
    use m_dhscf,        only: dhscf

    integer :: dummy_iscf = 1

    real(dp):: e1  ! Lower bound of energy range
    real(dp):: e2  ! Upper bound of energy range

    real(dp)  :: dummy_str(3,3), dummy_strl(3,3)  ! for dhscf call
    real(dp)  :: dummy_dipol(3)

    real(dp)  :: factor, g2max, dummy_Entrop

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    type(filesOut_t)           :: filesOut  ! blank output file names

    if ( .not. do_ldos ) return

#ifdef DEBUG
    call write_debug( '  PRE local_DOS' )
#endif

    ! Find local density of states
    if ( fdf_block('LocalDensityOfStates',bfdf) ) then

      ! Find the desired energy range
      if (.not. fdf_bline(bfdf,pline)) &
          call die('local_DOS: ERROR in LocalDensityOfStates block')

      if ( IONode ) write(*,'(/a)') 'siesta: LDOS info'

      if ( fdf_bmatch(pline, 'nvvn') ) then
        ! EF e1 e2 unit
        if ( .not. leqi(fdf_bnames(pline,1), 'Ef') ) then
          call die('local_DOS: ERROR in LocalDensityOfStates block, first name *must* be EF or not set')
        end if
        if ( IONode ) &
            write(*,'(a)') 'siesta: Shifting energies with respect to Fermi-level'
        
        factor = fdf_convfac( fdf_bnames(pline,2), 'Ry' )
        e1 = fdf_bvalues(pline,1)*factor + Ef 
        e2 = fdf_bvalues(pline,2)*factor + Ef
        
      else if ( fdf_bmatch(pline, 'vvn') ) then
        
        factor = fdf_convfac( fdf_bnames(pline,1), 'Ry' )
        e1 = fdf_bvalues(pline,1)*factor
        e2 = fdf_bvalues(pline,2)*factor
        
      else
        call die('local_DOS: ERROR in LocalDensityOfStates block!')
      end if

      ! Close block
      call fdf_bclose(bfdf)

      if ( IONode ) then
        write(*,'(a,tr1,f8.3," -- ",f8.3)') 'siesta: E1 -- E2 [eV]:', e1/eV, e2/eV
      end if

      ! If the k points have been set specifically for the LDOS then use this set
      if ( kpoint_ldos%N > kpoint_scf%N ) then
        call re_alloc(eo,1,no_u,1,spinor_dim,1,kpoint_ldos%N,name="eo", &
            routine="local_dos")
        call re_alloc(qo,1,no_u,1,spinor_dim,1,kpoint_ldos%N,name="qo", &
            routine="local_dos")
      end if

      ! Find the density matrix for states between e1 and e2
      call diagon(no_s, spinor_dim, no_l, maxnh, maxnh, no_u, &
          numh, listhptr, listh, numh, listhptr, listh, &
          H, S, qtot, fixspin, qtots, temp, e1, e2, &
          xijo, indxuo, gamma_ldos, &
          kpoint_ldos%N, kpoint_ldos%k, kpoint_ldos%w, &
          eo, qo, Dscf, Escf, ef, efs, dummy_Entrop, no_u, &
          occtol, dummy_iscf, neigwanted)

      ! Find the LDOS in the real space mesh
      filesOut%rho = trim(slabel) // '.LDOS'
      g2max = g2cut
      call dhscf( nspin, no_s, iaorb, iphorb, no_l, &
          no_u, na_u, na_s, isa, xa_last, indxua, &
          ntm, 0, 0, 0, filesOut, &
          maxnh, numh, listhptr, listh, Dscf, Datm, maxnh, H, &
          Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc, &
          dummy_dipol, dummy_str, fa, dummy_strl )
      
      ! next to last argument is dummy here,
      ! as no forces are calculated
      ! todo: make all these optional

    else

      call die('LDOS: something went terribly wrong')

    end if

#ifdef DEBUG
    call write_debug( '  POS local_DOS' )
#endif
  end subroutine local_DOS

end module m_local_DOS
