! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_projected_DOS

  use precision

  implicit none

  private

  public :: init_projected_DOS, projected_DOS

contains

  subroutine init_projected_DOS( ucell )
    
    use precision, only: dp
    use siesta_options
    use fdf,         only: fdf_isblock
    ! This is to get the reference kpoints in case PDOS.kgrid* has not
    ! been specified
    use kpoint_pdos_m, only: setup_kpoint_pdos
    use parallel, only: IOnode, Nodes
    
    real(dp), intent(in) :: ucell(3,3)
    
    ! Compute the projected density of states
    do_pdos = fdf_isblock('ProjectedDensityOfStates')
    if ( .not. do_pdos ) return
    
    if (isolve.ne.SOLVE_DIAGON) then
      if (.not.((isolve.eq.SOLVE_MINIM).and. minim_calc_eigenvalues)) then
        if (IONode) then
          write(*,*) 'siesta: ERROR: PDOS implemented only with diagon'
        end if
        do_pdos = .false.
      end if
    end if
    
    if ( .not. do_pdos ) return
    
    ! Determine whether the projected density of states is to be computed
    ! on a different grid to the SCF calculation
    call setup_kpoint_pdos( ucell )
      
  end subroutine init_projected_DOS
  
  subroutine projected_DOS()

    use sparse_matrices
    USE siesta_options
    use alloc,       only : re_alloc
    use atomlist,    only : indxuo, no_s, no_u, no_l
    use fdf
    use sys,         only : die
    use kpoint_scf_m, only: kpoint_scf
    use Kpoint_pdos_m, only: kpoint_pdos, gamma_pdos
    use parallel,    only: IOnode
    use m_energies, only: Ef
    use m_eo
    use m_spin,      only: h_spin_dim, spinor_dim
    use units,       only: eV
      
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    
    real(dp) :: factor
    logical  :: dummy ! Logical to hold return value from call to fdf_block
    integer  :: nhist ! Number of histogram intervals in projected DOS
    real(dp) :: e1    ! Lower bound of energy range
    real(dp) :: e2    ! Upper bound of energy range
    real(dp) :: sigma ! Energy width used to convolute partial DOS

    ! Compute the projected density of states

    if ( .not. do_PDOS ) return

#ifdef DEBUG
    call write_debug( '  PRE projected_DOS' )
#endif

    ! Call fdf_block to get iu - presence has already been tested in init_projected_PDOS
    if ( fdf_block('ProjectedDensityOfStates',bfdf) ) then

      ! Find the desired energy range
      if ( .not. fdf_bline(bfdf,pline) ) &
          call die('projected_DOS: ERROR in ProjectedDensityOfStates block')

      if ( IONode ) write(*,'(/a)') 'siesta: PDOS info'

      if ( fdf_bmatch(pline, 'nvvvin') ) then

        if ( .not. leqi(fdf_bnames(pline,1), 'Ef') ) then
          call die('projected_DOS: ERROR in ProjectedDensityOfStates block, first name *must* be EF or not set')
        end if
        if ( IONode ) &
            write(*,'(a)') 'siesta: Shifting energies with respect to Fermi-level'

        factor = fdf_convfac( fdf_bnames(pline,2), 'Ry' )
        e1 = fdf_breals(pline,1) * factor + Ef
        e2 = fdf_breals(pline,2) * factor + Ef

      else if ( fdf_bmatch(pline, 'vvvin') ) then
        
        factor = fdf_convfac( fdf_bnames(pline,1), 'Ry' )
        e1 = fdf_breals(pline,1) * factor
        e2 = fdf_breals(pline,2) * factor
        
      else
        call die("projected_DOS: ERROR Wrong format in PDOS block, not enough reals/integer/names")
      end if

      ! Get sigma and n-hist
      sigma = fdf_breals(pline,3) * factor
      nhist = fdf_bintegers(pline,1)

      ! Close block
      call fdf_bclose(bfdf)

      if ( IOnode ) then
        write(*,'(a,f8.3," -- ",2(f8.3,tr1),i0)') 'siesta: E1 -- E2, sigma [eV], nhist: ', &
            e1/eV, e2/eV, sigma/eV, nhist
      end if

      ! If the k points have been set specifically for the PDOS then use this set
      if ( kpoint_pdos%N > kpoint_scf%N ) then
        call re_alloc(eo,1,no_u,1,spinor_dim,1,kpoint_pdos%N,name="eo", routine="projected_dos")
      end if
      
      call pdos( no_s, h_spin_dim, spinor_dim, no_l, &
          maxnh, no_u, numh, listhptr, listh, H, S, &
          e1, e2, sigma, nhist, xijo, indxuo, gamma_pdos, &
          kpoint_pdos%N, kpoint_pdos%k, kpoint_pdos%w, eo, no_u)

    else

      call die('PDOS: something went terribly wrong')

    end if
    
#ifdef DEBUG
    call write_debug( '  POS projected_DOS' )
#endif

  end subroutine projected_DOS

end module m_projected_DOS
