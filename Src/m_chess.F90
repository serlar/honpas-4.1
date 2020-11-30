module m_chess

  use precision,      only : dp
#ifdef MPI
  use mpi_siesta,     only : mpi_integer, mpi_double_precision, mpi_status_size
  use parallelsubs,   only : GetNodeOrbs, GlobalToLocalOrb, WhichNodeOrb
  use parallelsubs,   only : set_BlockSizeDefault
#endif
  use parallel, only: Node, Nodes, BlockSize, IOnode

#ifdef SIESTA__CHESS
  ! The following modules used to be in sparsematrix_base...
  use dynamic_memory
  use dictionaries
  use yaml_output
  use yaml_strings
  use wrapper_MPI
  use wrapper_linalg
  use f_utils
  use numerics
  use time_profiling

  use sparsematrix_base
  use sparsematrix_highlevel, only: matrix_fermi_operator_expansion
  use foe_base, only: foe_data
#endif

  private

#ifdef SIESTA__CHESS

  !Public routines
  public :: CheSS_wrapper
  public :: CheSS_init
  public :: CheSS_finalize
  public :: get_CheSS_parameter
  public :: set_CheSS_parameter

  ! Global variables, private to this module
  integer :: nvctr, nvctr_kernel, nvctr_mult
  type(sparse_matrix),dimension(2) :: smat
  type(matrices) :: mat_h, mat_s, mat_k, mat_ek
  type(matrices),dimension(1) :: mat_ovrlpminusonehalf
  integer,dimension(:),allocatable :: numh_global, conversion_lookup
  integer,dimension(:),allocatable :: numh_global_kernel, conversion_lookup_kernel
  integer,dimension(:),allocatable :: numh_global_mult, conversion_lookup_mult
  type(foe_data) :: foe_obj, ice_obj
  real(mp) :: chess_buffer_kernel, chess_buffer_mult, chess_betax
  real(mp) :: chess_fscale, chess_fscale_lowerbound, chess_fscale_upperbound
  real(mp) :: chess_evlow_h, chess_evhigh_h, chess_evlow_s, chess_evhigh_s
  type(dictionary),pointer :: dict_timing_info
  external :: gather_timings

  contains



    subroutine CheSS_init(node, nodes, nhmax, nhmax_kernel, nhmax_mult, h_dim, &
               nbasis, nbasis_kernel, nbasis_mult, BlockSize, nspin, qs, &
               listh, listh_kernel, listh_mult, numh, numh_kernel, numh_mult)
      use sparsematrix_base
      use sparsematrix_highlevel, only: sparse_matrix_init_from_data_ccs, matrices_init
      use sparsematrix_init, only: write_sparsematrix_info, init_matrix_taskgroups_wrapper
      use foe_common, only: init_foe
      implicit none

      ! Calling arguments
      integer,intent(in) :: node, nodes, nhmax, nhmax_kernel, nhmax_mult, h_dim
      integer,intent(in) :: nbasis, nbasis_kernel, nbasis_mult, BlockSize, nspin
      real(dp),dimension(1:nspin),intent(in) :: qs
      integer,dimension(nhmax),intent(in) :: listh
      integer,dimension(nbasis),intent(in) :: numh
      integer,dimension(nhmax_kernel),intent(in) :: listh_kernel
      integer,dimension(nbasis_kernel),intent(in) :: numh_kernel
      integer,dimension(nhmax_mult),intent(in) :: listh_mult
      integer,dimension(nbasis_mult),intent(in) :: numh_mult

      ! Local variables
      integer,dimension(:),allocatable :: row_ind_global, col_ptr_global
      integer,dimension(:),allocatable :: row_ind_global_kernel, col_ptr_global_kernel
      integer,dimension(:),allocatable :: row_ind_global_mult, col_ptr_global_mult
      integer :: i

      ! Initialize flib
      call f_lib_initialize()

      if (IOnode) then
          call yaml_new_document()
      end if

      ! Initialize the sparsematrix error handling and timing.
      call sparsematrix_init_errors()
      call sparsematrix_initialize_timing_categories()

      call f_timing_reset(filename='time.yaml',master=(node==0),verbose_mode=.false.)

      ! Gather together the matrices and descriptors from all MPI tasks
      nvctr = nhmax
      call mpiallred(nvctr, 1, mpi_sum, comm=mpi_comm_world)
      row_ind_global = f_malloc0(nvctr,id='row_ind_global')
      numh_global = f_malloc0(h_dim,id='numh_global')
      col_ptr_global = f_malloc0(h_dim,id='col_ptr_global')
      conversion_lookup = f_malloc(nhmax,id='conversion_lookup')

      nvctr_kernel = nhmax_kernel
      call mpiallred(nvctr_kernel, 1, mpi_sum, comm=mpi_comm_world)
      row_ind_global_kernel = f_malloc0(nvctr_kernel,id='row_ind_global_kernel')
      numh_global_kernel = f_malloc0(h_dim,id='numh_global_kernel')
      col_ptr_global_kernel = f_malloc0(h_dim,id='col_ptr_global_kernel')
      conversion_lookup_kernel = f_malloc(nhmax_kernel,id='conversion_lookup_kernel')

      nvctr_mult = nhmax_mult
      call mpiallred(nvctr_mult, 1, mpi_sum, comm=mpi_comm_world)
      row_ind_global_mult = f_malloc0(nvctr_mult,id='row_ind_global_mult')
      numh_global_mult = f_malloc0(h_dim,id='numh_global_mult')
      col_ptr_global_mult = f_malloc0(h_dim,id='col_ptr_global_mult')
      conversion_lookup_mult = f_malloc(nhmax_mult,id='conversion_lookup_mult')



      call gather_numh(node, nodes, h_dim, nbasis, BlockSize, numh, numh_global)
      call gather_numh(node, nodes, h_dim, nbasis_kernel, BlockSize, numh_kernel, numh_global_kernel)
      call gather_numh(node, nodes, h_dim, nbasis_mult, BlockSize, numh_mult, numh_global_mult)

      !!write(*,*) 'numh_global',numh_global
      !!write(*,*) 'numh_global_mult',numh_global_mult

      !!write(*,*) 'listh',listh
      !!write(*,*) 'listh_mult',listh_mult
      !!write(*,*) 'nhmax, nhmax_mult', nhmax, nhmax_mult

      call gather_row_ind(node, nodes, h_dim, nhmax, nvctr, &
           BlockSize, numh_global, listh, row_ind_global)
      call gather_row_ind(node, nodes, h_dim, nhmax_kernel, nvctr_kernel, &
           BlockSize, numh_global_kernel, listh_kernel, row_ind_global_kernel)
      call gather_row_ind(node, nodes, h_dim, nhmax_mult, nvctr_mult, &
           BlockSize, numh_global_mult, listh_mult, row_ind_global_mult)

      !!write(*,*) 'row_ind_global',row_ind_global
      !!write(*,*) 'row_ind_global_mult',row_ind_global_mult

      call get_gol_ptr_global(h_dim, numh_global, col_ptr_global)
      call get_gol_ptr_global(h_dim, numh_global_kernel, col_ptr_global_kernel)
      call get_gol_ptr_global(h_dim, numh_global_mult, col_ptr_global_mult)

      call sparse_matrix_init_from_data_ccs(node, nodes, mpi_comm_world, &
           h_dim, nvctr, row_ind_global, col_ptr_global, smat(1), &
           init_matmul=.false.)

      call sparse_matrix_init_from_data_ccs(node, nodes, mpi_comm_world, &
           h_dim, nvctr_kernel, row_ind_global_kernel, col_ptr_global_kernel, smat(2), &
           init_matmul=.true., nvctr_mult=nvctr_mult, row_ind_mult=row_ind_global_mult, col_ptr_mult=col_ptr_global_mult)

      call init_matrix_taskgroups_wrapper(node, nodes, mpi_comm_world, .false., 2, smat)

      if (IOnode) then
          call yaml_mapping_open('Sparse matrices')
          call write_sparsematrix_info(smat(1), 'Overlap and hamiltonian matrix')
          call write_sparsematrix_info(smat(2), 'Density kernel matrix')
          call yaml_mapping_close()
      end if

      ! Initialize the object holding some parameters for FOE.
      call init_foe(node, nodes, nspin, qs, foe_obj, &
           fscale=chess_fscale, &
           fscale_lowerbound=chess_fscale_lowerbound, &
           fscale_upperbound=chess_fscale_upperbound, &
           evlow=chess_evlow_h, &
           evhigh=chess_evhigh_h, &
           betax=chess_betax)
      ! Initialize the same object for the calculation of the inverse. Charge does not really make sense here...
      call init_foe(node, nodes, nspin, qs, ice_obj, &
           evlow=chess_evlow_s, &
           evhigh=chess_evhigh_s, &
           betax=chess_betax)

      call f_free(row_ind_global)
      call f_free(col_ptr_global)
      call f_free(row_ind_global_kernel)
      call f_free(col_ptr_global_kernel)
      call f_free(row_ind_global_mult)
      call f_free(col_ptr_global_mult)

      ! Allocate the matrices
      call matrices_init(smat(1), mat_h, matsize=SPARSE_TASKGROUP)
      call matrices_init(smat(1), mat_s, matsize=SPARSE_TASKGROUP)
      call matrices_init(smat(2), mat_k, matsize=SPARSE_TASKGROUP)
      call matrices_init(smat(2), mat_ek, matsize=SPARSE_TASKGROUP)
      call matrices_init(smat(2), mat_ovrlpminusonehalf(1), matsize=SPARSE_TASKGROUP)

      call f_timing_checkpoint(ctr_name='INIT',mpi_comm=mpiworld(),nproc=mpisize(), &
           gather_routine=gather_timings)

    end subroutine CheSS_init


    subroutine CheSS_finalize()
      use sparsematrix_base
      use foe_base, only: foe_data_deallocate
      implicit none

      call f_timing_checkpoint(ctr_name='CALC',mpi_comm=mpiworld(),nproc=mpisize(), &
           gather_routine=gather_timings)

      call f_free(numh_global)
      call f_free(conversion_lookup)
      call f_free(numh_global_kernel)
      call f_free(conversion_lookup_kernel)
      call f_free(numh_global_mult)
      call f_free(conversion_lookup_mult)
      call deallocate_sparse_matrix(smat(1))
      call deallocate_sparse_matrix(smat(2))
      call foe_data_deallocate(ice_obj)
      call foe_data_deallocate(foe_obj)

      call deallocate_matrices(mat_h)
      call deallocate_matrices(mat_s)
      call deallocate_matrices(mat_k)
      call deallocate_matrices(mat_ek)
      call deallocate_matrices(mat_ovrlpminusonehalf(1))

      ! Gather the timings
      call build_dict_info(node, nodes, dict_timing_info)
      call f_timing_stop(mpi_comm=mpi_comm_world, nproc=nodes, &
           gather_routine=gather_timings, dict_info=dict_timing_info)
      call dict_free(dict_timing_info)

       if (IOnode) then
           call yaml_release_document()
           call f_lib_finalize()
       end if

    end subroutine CheSS_finalize


    subroutine CheSS_wrapper(CalcE,PreviousCallDiagon,iscf,istp,nbasis,nspin,&
               h_dim,nhmax,numh,listhptr,listh,&
               qs,h_sparse,s_sparse,&
               d_sparse,ed_sparse,Ef)
      use sparsematrix_base
      use foe_base, only: foe_data, foe_data_deallocate, foe_data_get_real
      use foe_common, only: init_foe
      use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_ccs, &
                                        sparse_matrix_init_from_file_ccs, &
                                        matrices_get_values, matrices_set_values, &
                                        sparse_matrix_init_from_data_ccs, &
                                        ccs_data_from_sparse_matrix, ccs_matrix_write, &
                                        matrix_matrix_multiplication, matrix_fermi_operator_expansion, &
                                        trace_AB, trace_A
      use sparsematrix_init, only: write_sparsematrix_info
      use sparsematrix, only: write_matrix_compressed, transform_sparse_matrix, get_minmax_eigenvalues
      use futile

      implicit none
    
      !**** INPUT ***********************************!
    
      logical, intent(in) :: CalcE              ! Calculate the energy-density matrix from the existing coeffs.?
      logical, intent(in) :: PreviousCallDiagon ! Previous SCF iteration solved by diagonalization?
    
      integer, intent(in) :: iscf               ! SCF iteration num.
      integer, intent(in) :: istp               ! MD iteration num.
      integer, intent(in) :: nbasis             ! dimension of numh and listhptr
      integer, intent(in) :: nspin              ! num. of spins
      integer, intent(in) :: h_dim              ! num. of AOs (global)
      integer, intent(in) :: nhmax              ! first dimension of listh and sparse matrices
      integer, intent(in) :: numh(1:nbasis)     ! num. of nonzero elements of each row of sparse matrices
      integer, intent(in) :: listhptr(1:nbasis) ! pointer to start of row in listh
      integer, intent(in) :: listh(1:nhmax)     ! list of nonzero elements of each row of sparse matrices
      real(dp), intent(in) :: qs(1:nspin)       ! num. of electrons per spin
      real(dp), intent(in) :: h_sparse(1:nhmax,1:nspin) ! hamiltonian matrix (sparse)
      real(dp), intent(in) :: s_sparse(1:nhmax)         ! overlap matrix (sparse)
    
      !**** OUTPUT **********************************!
      real(dp), intent(out) :: d_sparse(1:nhmax,1:nspin) ! density matrix (sparse)
      real(dp), intent(out) :: ed_sparse(1:nhmax,1:nspin) ! energy-density matrix (sparse)
      real(dp), intent(out) :: Ef ! Fermi energy

      !**** LOCAL ***********************************!
      integer :: ibasis, ind, i, ii, ist, jproc, jjproc, ist_global, jorb, iblock
      integer,dimension(:),allocatable :: lookup, row_ind_global, col_ptr_global
      real(dp),dimension(:),allocatable :: charge, s_work, h_work, k_work, k_small
      real(dp) :: energy, eval_min, eval_max



      if (nspin/=1) then
          call f_err_throw('The Chess coupling is only tested for nspin=1')
      end if

      s_work = f_malloc(nhmax,id='s_work')
      h_work = f_malloc(nhmax,id='h_work')
      k_work = f_malloc(nhmax,id='k_work')
      k_small = f_malloc(nvctr,id='k_small')
      call order_matrices(nbasis, nhmax, numh, listh, s_sparse, h_sparse, conversion_lookup, s_work, h_work)

      call get_matrices(node, nodes, h_dim, BlockSize, nhmax, nvctr, numh_global, &
           s_work, h_work, mat_s%matrix_compr, mat_h%matrix_compr)


      !!call get_minmax_eigenvalues(node, smat(1), mat_s, eval_min, eval_max)
      !!call get_minmax_eigenvalues(node, smat(1), mat_h, eval_min, eval_max)

      ! Only calculate S^-1/2 in the very first iteration, as later it does not change any more
      call matrix_fermi_operator_expansion(node, nodes, mpi_comm_world, &
           foe_obj, ice_obj, smat(1), smat(1), smat(2), &
           mat_s, mat_h, mat_ovrlpminusonehalf, mat_k, energy, &
           calculate_minusonehalf=(iscf==1), foe_verbosity=1, symmetrize_kernel=.true., &
           calculate_energy_density_kernel=.true., energy_kernel=mat_ek)

      Ef = foe_data_get_real(foe_obj,"ef",1)

      ! Transform the density kernel to the small SIESTA sparsity format
      call transform_sparse_matrix(node, smat(1), smat(2), SPARSE_FULL, 'large_to_small', &
           lmat_in=mat_k%matrix_compr, smat_out=k_small)
      call set_matrices(node, nodes, h_dim, BlockSize, nhmax, nvctr, numh_global, k_small, k_work)
      call unorder_matrices(nhmax, conversion_lookup, k_work, d_sparse)

      call transform_sparse_matrix(node, smat(1), smat(2), SPARSE_FULL, 'large_to_small', &
           lmat_in=mat_ek%matrix_compr, smat_out=k_small)
      call set_matrices(node, nodes, h_dim, BlockSize, nhmax, nvctr, numh_global, k_small, k_work)
      call unorder_matrices(nhmax, conversion_lookup, k_work, ed_sparse)


      call f_free(s_work)
      call f_free(h_work)
      call f_free(k_work)
      call f_free(k_small)


    end subroutine CheSS_wrapper


    subroutine gather_numh(node, nodes, h_dim, nbasis, BlockSize, numh, numh_global)
      use sparsematrix_base
      implicit none

      ! Calling arguments
      integer,intent(in) :: node, nodes, h_dim, nbasis, BlockSize
      integer,dimension(nbasis),intent(in) :: numh
      integer,dimension(h_dim),intent(out) :: numh_global

      ! Local variables
      integer :: ist, jproc, ist_global, jorb, jjproc

      ist = 0
      jproc=0
      ist_global = 0
      do jorb=1,h_dim,BlockSize
          jjproc = mod(jproc,nodes)
          if (node==jjproc) then
              !numh_global(ist_global+1:min(ist_global+BlockSize,nbasis)) = numh(ist+1:min(ist+BlockSize,h_dim))
              numh_global(ist_global+1:min(ist_global+BlockSize,h_dim)) = numh(ist+1:min(ist+BlockSize,nbasis))
              ist = ist + BlockSize
          end if
          ist_global = ist_global + BlockSize
          jproc = jproc + 1
      end do
      call mpiallred(numh_global, mpi_sum, comm=mpi_comm_world)

    end subroutine gather_numh

    subroutine gather_row_ind(node, nodes, h_dim, nhmax, nvctr, BlockSize, numh_global, listh, row_ind_global)
      use sparsematrix_base
      implicit none

      ! Calling arguments
      integer,intent(in) :: node, nodes, h_dim, nhmax, nvctr, BlockSize
      integer,dimension(h_dim),intent(in) :: numh_global
      integer,dimension(nhmax),intent(in) :: listh
      integer,dimension(nvctr),intent(out) :: row_ind_global

      ! Local variables
      integer :: ist, ist_global, jorb, jproc, jjproc, iblock

      ist = 0
      ist_global = 0
      jproc=0
      do jorb=1,h_dim,BlockSize
          jjproc = mod(jproc,nodes)
          do iblock=jorb,min(jorb+BlockSize-1,h_dim)
              if (node==jjproc) then
                  row_ind_global(ist_global+1:ist_global+numh_global(iblock)) = listh(ist+1:ist+numh_global(iblock))
                  ist = ist + numh_global(iblock)
              end if
              ist_global = ist_global + numh_global(iblock)
          end do
          jproc = jproc + 1
      end do
      call mpiallred(row_ind_global, mpi_sum, comm=mpi_comm_world)

    end subroutine gather_row_ind


    subroutine get_gol_ptr_global(h_dim, numh_global, col_ptr_global)
      implicit none

      ! Calling arguments
      integer,intent(in) :: h_dim
      integer,dimension(h_dim),intent(in) :: numh_global
      integer,dimension(h_dim),intent(out) :: col_ptr_global

      ! Local variables
      integer :: ii, ibasis
      ii = 1
      do ibasis=1,h_dim
          col_ptr_global(ibasis) = ii
          ii = ii + numh_global(ibasis)
      end do
    end subroutine get_gol_ptr_global


    subroutine order_matrices(nbasis, nhmax, numh, listh, s_sparse, h_sparse, conversion_lookup, s_work, h_work)
      use sparsematrix_base
      implicit none

      ! Calling arguments
      integer,intent(in) :: nbasis, nhmax
      integer,dimension(nbasis),intent(in) :: numh
      integer,dimension(nhmax),intent(in) :: listh
      real(dp),dimension(nhmax),intent(in) :: s_sparse
      real(dp),dimension(nhmax),intent(in) :: h_sparse
      integer,dimension(nhmax),intent(out) :: conversion_lookup
      real(dp),dimension(nhmax),intent(out) :: s_work, h_work

      ! Local variables
      integer :: ist, ii, i, ind, ibasis
      integer,dimension(:),allocatable :: lookup

      lookup = f_malloc(maxval(numh),id='lookup')
      ist = 0
      ii = 0
      do ibasis=1,nbasis
          lookup(1:numh(ibasis)) = listh(ist+1:ist+numh(ibasis))
          do i=1,numh(ibasis)
              ii = ii + 1
              ind = minloc(lookup(1:numh(ibasis)),1)
              s_work(ii) = s_sparse(ist+ind)
              h_work(ii) = h_sparse(ist+ind)
              conversion_lookup(ist+ind) = ii
              lookup(ind) = huge(ind)
          end do
          ist = ist + numh(ibasis)
      end do

      call f_free(lookup)

    end subroutine order_matrices


    subroutine get_matrices(node, nodes, h_dim, BlockSize, nhmax, nvctr, numh_global, s_work, h_work, s_global, h_global)
      use sparsematrix_base
      implicit none

      ! Calling arguments
      integer,intent(in) :: node, nodes, h_dim, BlockSize, nhmax, nvctr
      integer,dimension(h_dim),intent(in) :: numh_global
      real(dp),dimension(nhmax),intent(in) :: s_work, h_work
      real(dp),dimension(nvctr),intent(out) :: s_global, h_global

      ! Local variables
      integer :: ist, ist_global, jorb, jproc, jjproc, iblock 

      jproc = 0
      ist = 0
      ist_global = 0
      call f_zero(s_global)
      call f_zero(h_global)
      do jorb=1,h_dim,BlockSize
          jjproc = mod(jproc,nodes)
          do iblock=jorb,min(jorb+BlockSize-1,h_dim)
              if (node==jjproc) then
                  !write(*,*) 'node, jorb, iblock, jjproc, ist, ist_global', node, jorb, iblock, jjproc, ist, ist_global
                  s_global(ist_global+1:ist_global+numh_global(iblock)) = s_work(ist+1:ist+numh_global(iblock))
                  h_global(ist_global+1:ist_global+numh_global(iblock)) = h_work(ist+1:ist+numh_global(iblock))
                  ist = ist + numh_global(iblock)
              end if
              ist_global = ist_global + numh_global(iblock)
          end do
          jproc = jproc + 1
      end do
      call mpiallred(s_global, mpi_sum, comm=mpi_comm_world)
      call mpiallred(h_global, mpi_sum, comm=mpi_comm_world)
    end subroutine get_matrices


    subroutine set_matrices(node, nodes, h_dim, BlockSize, nhmax, nvctr, numh_global, k_global, k_work)
      implicit none

      ! Calling arguments
      integer,intent(in) :: node, nodes, h_dim, BlockSize, nhmax, nvctr
      integer,dimension(h_dim),intent(in) :: numh_global
      real(dp),dimension(nvctr),intent(out) :: k_global
      real(dp),dimension(nhmax),intent(out) :: k_work

      ! Local variables
      integer :: ist, ist_global, jorb, jproc, jjproc, iblock

      jproc = 0
      ist = 0
      ist_global = 0
      do jorb=1,h_dim,BlockSize
          jjproc = mod(jproc,nodes)
          do iblock=jorb,min(jorb+BlockSize-1,h_dim)
              if (node==jjproc) then
                  k_work(ist+1:ist+numh_global(iblock)) = k_global(ist_global+1:ist_global+numh_global(iblock))
                  ist = ist + numh_global(iblock)
              end if
              ist_global = ist_global + numh_global(iblock)
          end do
          jproc = jproc + 1
      end do

    end subroutine set_matrices


    subroutine unorder_matrices(nhmax, conversion_lookup, k_work, d_sparse)
      implicit none

      ! Calling arguments
      integer,intent(in) :: nhmax
      integer,dimension(nhmax),intent(in) :: conversion_lookup
      real(dp),dimension(nhmax),intent(in) :: k_work
      real(dp),dimension(nhmax),intent(out) :: d_sparse

      ! Local variables
      integer :: ibasis, ind
      do ibasis=1,nhmax
          ind = conversion_lookup(ibasis)
          d_sparse(ibasis) = k_work(ind)
      end do

    end subroutine unorder_matrices


    subroutine set_CheSS_parameter(key, val)
      implicit none

      ! Calling arguments
      character(len=*),intent(in) :: key
      real(mp),intent(in) :: val

      select case (key)
      case ('chess_buffer_kernel')
          chess_buffer_kernel = val
      case ('chess_buffer_mult')
          chess_buffer_mult = val
      case ('chess_betax')
          chess_betax = val
      case ('chess_fscale')
          chess_fscale = val
      case ('chess_fscale_lowerbound')
          chess_fscale_lowerbound = val
      case ('chess_fscale_upperbound')
          chess_fscale_upperbound = val
      case ('chess_evlow_h')
          chess_evlow_h = val
      case ('chess_evhigh_h')
          chess_evhigh_h = val
      case ('chess_evlow_s')
          chess_evlow_s = val
      case ('chess_evhigh_s')
          chess_evhigh_s = val
      case default
          call f_err_throw('wrong key in set_CheSS_parameter')
      end select

    end subroutine set_CheSS_parameter


    function get_CheSS_parameter(key) result(val)
      implicit none

      ! Calling arguments
      character(len=*),intent(in) :: key
      real(mp) :: val

      select case (key)
      case ('chess_buffer_kernel')
          val = chess_buffer_kernel
      case ('chess_buffer_mult')
          val = chess_buffer_mult
      case ('chess_betax')
          val = chess_betax
      case ('chess_fscale')
          val = chess_fscale
      case ('chess_fscale_lowerbound')
          val = chess_fscale_lowerbound
      case ('chess_fscale_upperbound')
          val = chess_fscale_upperbound
      case ('chess_evlow_h')
          val = chess_evlow_h
      case ('chess_evhigh_h')
          val = chess_evhigh_h
      case ('chess_evlow_s')
          val = chess_evlow_s
      case ('chess_evhigh_s')
          val = chess_evhigh_s
      case default
          call f_err_throw('wrong key in set_CheSS_parameter')
      end select

    end function get_CheSS_parameter

#endif

end module m_chess
