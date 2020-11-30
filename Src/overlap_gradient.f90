! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module overlap_gradient_m
    
  implicit none

  private
  public :: overlap_gradient
  
contains
  
  subroutine overlap_gradient(na_u, lasto, isa, nnz, xijo, grS_2D)
    
    use precision, only : dp
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use geom_helper, only : iaorb, ucorb
    use atmfuncs,      only : orb_gindex
    use m_new_matel,   only : new_matel
    
    ! *********************************************************************
    ! Calculate the gradient of the overlap matrix
    ! Energies in Ry. Lengths in Bohr.
    ! Writen by Nick Papior, July '18
    ! **************************** INPUT **********************************
    ! integer na_u             : Number of atoms in unit cell
    ! integer lasto(0:na_u)    : Last orbital index of each atom
    ! integer isa(na_u)        : Species index of each atom
    ! integer nnz              : Number of non-zero elements in the sparse matrix
    ! real*8  xij(3,nnz)       : orbital distance vectors
    ! ********************** INPUT and OUTPUT *****************************
    ! real*8  grS_2D(:,3)      : Gradient of overlap matrix
    ! *********************************************************************
    
    integer, intent(in) :: na_u, nnz
    integer, intent(in) :: isa(na_u), lasto(0:na_u)
    real(dp), intent(in) :: xijo(3,nnz)
    type(dSpData2D), intent(inout) :: grS_2D
    
    ! Internal variables ......................................................
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: s
    real(dp), pointer :: grS(:,:)

    integer :: no_l, no_u
    integer, pointer :: ncol(:), l_ptr(:), col(:)
    integer :: lio, ind
    integer :: io, ia, ioa, is, ig
    integer :: jo, ja, joa, js, jg
    
    real(dp) :: Sij, grSij(3)

    ! Start timer
    call timer( 'overlap_grad', 1 )

    dit => dist(grS_2D)
    s => spar(grS_2D)
    grS => val(grS_2D)
    call attach(s, nrows=no_l, nrows_g=no_u, n_col=ncol, list_ptr=l_ptr, list_col=col)

    do lio = 1, no_l
      io = index_local_to_global(dit, lio)
      ! Orbital on atom
      ia = iaorb(io,lasto) ! atom-index
      ioa = io - lasto(ia-1)
      is = isa(ia) ! species number
      ! orbital g-index
      ig = orb_gindex(is,ioa)

      ! Loop connection graph
      do ind = l_ptr(lio) + 1, l_ptr(lio) + ncol(lio)
        ! Orbital on atom
        jo = ucorb(col(ind),no_u)
        ja = iaorb(jo,lasto) ! atom-index
        joa = jo - lasto(ja-1)
        js = isa(ja) ! species number
        ! orbital g-index
        jg = orb_gindex(js,joa)

        call new_MATEL('S', ig, jg, xijo(:,ind), Sij, grSij )
        grS(ind, 1) = grSij(1)
        grS(ind, 2) = grSij(2)
        grS(ind, 3) = grSij(3)
        
      end do
      
    end do
    
    ! Finish timer
    call timer( 'overlap_grad', 2 )
    
  end subroutine overlap_gradient

end module overlap_gradient_m




