! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module kpoint_t_m

  use precision, only: dp

  implicit none

  public :: kpoint_t

  !< Method specification for non-specified k-list
  integer, parameter, public :: K_METHOD_NONE = 0
  !< Method specification for Monkhorst-Pack grid
  integer, parameter, public :: K_METHOD_MONKHORST_PACK = 1
  !< Method specification for cutoff designation
  integer, parameter, public :: K_METHOD_CUTOFF = 2
  !< Method specification for user-defined list
  integer, parameter, public :: K_METHOD_LIST = 3

  type :: kpoint_t
    
    !< Total number of k-points
    integer :: N = 0
    !< K-points
    real(dp), pointer :: k(:,:) => null()
    !< weights for k-points
    real(dp), pointer :: w(:) => null()
    !< whether these k-points are under time-reversal symmetry
    logical :: trs = .true.
    
    !< The method by which these k-points are generated
    integer :: method = K_METHOD_MONKHORST_PACK

    ! Specific elements for Monkhorst-Pack, cutoff grids
    integer :: k_cell(3,3) = 0
    real(dp) :: k_displ(3) = 0._dp

    ! The cutoff requested
    real(dp) :: cutoff = 0._dp

  end type kpoint_t

  public :: kpoint_associated
  public :: kpoint_associate
  public :: kpoint_nullify
  public :: kpoint_delete
  public :: kpoint_read
  
  public :: kpoint_read_file
  public :: kpoint_write_xml
  public :: kpoint_write_stdout
  public :: kpoint_write_file

contains

  !< Associate one k-point list with another (no copies made)
  subroutine kpoint_associate(out, in)
    type(kpoint_t), intent(inout) :: out
    type(kpoint_t), intent(in) :: in

    out%N = in%N
    out%k => in%k
    out%w => in%w
    out%trs = in%trs
    out%method = in%method
    out%k_cell = in%k_cell
    out%k_displ = in%k_displ

  end subroutine kpoint_associate

  !< Figure out if two types are associated with each other
  function kpoint_associated(a, b) result(assoc)
    type(kpoint_t), intent(in) :: a, b
    logical :: assoc

    assoc = associated(a%k, b%k) .and. &
        associated(a%w, b%w)

  end function kpoint_associated

  !< Nullify k-point list
  subroutine kpoint_nullify(this)
    type(kpoint_t), intent(inout) :: this

    nullify(this%k, this%w)
    call kpoint_delete(this)

  end subroutine kpoint_nullify

  !< Delete the k-point list
  subroutine kpoint_delete(this)
    use alloc, only: de_alloc
    type(kpoint_t), intent(inout) :: this

    call de_alloc(this%k, name='kpoint_t%k')
    call de_alloc(this%w, name='kpoint_t%w')
    this%N = 0
    this%trs = .true.
    this%method = K_METHOD_NONE
    this%k_cell = 0
    this%k_displ = 0
    
  end subroutine kpoint_delete

  !< Read settings for this k-point grid
  !!
  !! The order of reading the k-points are as follows:
  !!
  !! 1. kgrid.MonkhorstPack has precedence
  !!    If this block is found the block will be read as this:
  !!     %block kgrid.MonkhorstPack  # Defines k_cell and k_displ
  !!       4  0  0   0.50               # (k_cell(i,1),i=1,3), k_displ(1)
  !!       0  4  0   0.50               # (k_cell(i,2),i=1,3), k_displ(2)
  !!       0  0  4   0.50               # (k_cell(i,3),i=1,3), k_displ(3)
  !!     %endblock kgrid.MonkhorstPack
  !!    Note that the displacement defaults to zero and is an optional value.
  !!    One may also write the block like this:
  !!     %block kgrid.MonkhorstPack  # Defines k_cell and k_displ
  !!       4 0.50               # k_cell(1,1), k_displ(1)
  !!       4 0.50               # k_cell(2,2), k_displ(2)
  !!       4 0.50               # k_cell(3,3), k_displ(3)
  !!     %endblock kgrid.MonkhorstPack
  !! 2. kgrid.Cutoff
  !!    The k_cell is specified using a cutoff parameter and automatically
  !!    calculated.
  !! 3. kgrid.file
  !!    The k-points are user-supplied. In this case the k-points *must*
  !!    be provided in units of the reciprocal lattice vectors (i.e. in ]-0.5 ; 0.5])
  subroutine kpoint_read(this, prefix, cell, trs, process_k_cell)
    use fdf
    
    type(kpoint_t), intent(inout) :: this
    character(len=*), intent(in) :: prefix
    real(dp), intent(in) :: cell(3,3)
    logical, intent(in) :: trs
    abstract interface
      subroutine my_sub(k_cell, k_displ)
        use precision, only: dp
        integer, intent(inout) :: k_cell(3,3)
        real(dp), intent(inout) :: k_displ(3)
      end subroutine
    end interface
    procedure(my_sub), optional :: process_k_cell

    ! For reading the k-cell etc.
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline

    character(len=256) :: name_fdf, file
    real(dp) :: cutoff
    integer :: i, nvals
    logical :: spiral

    ! If this already has k-points associated we continue
    if ( this%N > 0 ) return

    ! Be sure to delete everything
    call kpoint_delete(this)
    this%method = K_METHOD_NONE

    spiral = fdf_get('Spin.Spiral', .false.)
    ! Allow the user to control the use of time-reversal-symmetry
    ! By default, it is on, except for "spin-spiral" calculations
    ! and/or non-collinear calculations
    this%trs = fdf_get("TimeReversalSymmetryForKpoints", (.not. spiral) .and. trs)


    call kpoint_fdf_name(prefix, 'kgrid.MonkhorstPack', name_fdf)
    if ( fdf_block(name_fdf, bfdf) ) then

      ! The method is a k-point MP grid
      this%method = K_METHOD_MONKHORST_PACK

      ! Now read data
      do i = 1, 3
        if ( .not. fdf_bline(bfdf,pline) ) then
          call die('kpoint_read: ERROR in ' // trim(name_fdf) // ' block')
        end if

        ! Read displacement
        nvals = fdf_bnvalues(pline)
        if ( nvals > 3 ) then
          this%k_displ(i) = mod(fdf_bvalues(pline,4), 1._dp)
        else
          this%k_displ(i) = 0._dp
        end if

        if ( nvals == 1 ) then
          this%k_cell(i,i) = fdf_bintegers(pline,1)
          
        else if ( nvals == 2 ) then
          this%k_cell(i,i) = fdf_bintegers(pline,1)
          this%k_displ(i) = mod(fdf_bvalues(pline,2), 1._dp)
          
        else
          this%k_cell(1,i) = fdf_bintegers(pline,1)
          this%k_cell(2,i) = fdf_bintegers(pline,2)
          this%k_cell(3,i) = fdf_bintegers(pline,3)
          
        end if
        
      end do

      call fdf_bclose(bfdf)

      call setup_mp()

    else if ( fdf_islist(name_fdf) ) then

      ! The method is a k-point MP grid
      this%method = K_METHOD_MONKHORST_PACK

      ! We have a list of integers
      i = -1
      call fdf_list(name_fdf, i, this%k_cell(:,1))
      if ( i /= 3 ) then
        call die('kpoint_read: ERROR in ' // trim(name_fdf) // ' list (have to have 3 numbers)')
      end if
      
      ! Initialize MP
      this%k_cell = 0
      this%k_displ = 0._dp
      call fdf_list(name_fdf, i, this%k_cell(:, 1))
      
      ! Re-arange the elements
      this%k_cell(2,2) = this%k_cell(2,1)
      this%k_cell(2,1) = 0
      this%k_cell(3,3) = this%k_cell(3,1)
      this%k_cell(3,1) = 0

      call setup_mp()

    end if

    ! Quick return
    if ( this%method /= K_METHOD_NONE ) return

    call kpoint_fdf_name(prefix, 'kgrid.Cutoff', name_fdf)
    if ( fdf_defined(name_fdf) ) then

      ! The method is cutoff based
      this%method = K_METHOD_CUTOFF
      cutoff = fdf_get(name_fdf, -1._dp, 'Bohr')
      if ( cutoff <= 0._dp ) then
        this%method = K_METHOD_NONE
      else
        call setup_cutoff()
      end if

    end if

    ! Quick return
    if ( this%method /= K_METHOD_NONE ) return

    call kpoint_fdf_name(prefix, 'kgrid.File', name_fdf)
    if ( fdf_defined(name_fdf) ) then

      ! A user-defined list of k-points
      this%method = K_METHOD_LIST
      file = fdf_get(name_fdf, '01234567890123456789')

      call setup_file()

    end if

    ! Quick return
    if ( this%method /= K_METHOD_NONE ) return

    call setup_gamma()

  contains

    subroutine setup_mp()
      use m_find_kgrid, only: find_kgrid
      
      if ( present(process_k_cell) ) then
        call process_k_cell(this%k_cell, this%k_displ)
      end if

      ! Find the grid
      call find_kgrid(cell, this%k_cell, this%k_displ, .true., &
          this%trs, &
          this%N, this%k, this%w, this%cutoff)
      
    end subroutine setup_mp
    
    subroutine setup_cutoff()
      use m_minvec, only : minvec
      use m_find_kgrid, only: find_kgrid

      real(dp) :: scmin(3,3), vmod
      real(dp) :: ctransf(3,3)
      integer :: factor, expansion_factor

      ! Find equivalent rounded unit-cell
      call minvec( cell, scmin, ctransf )
      
      expansion_factor = 1
      do i = 1 , 3
        
        vmod = sqrt( dot_product(scmin(:,i),scmin(:,i)) )
        factor = int(2.0_dp * cutoff / vmod) + 1
        expansion_factor = expansion_factor * factor
        
        ! Generate actual supercell skeleton
        this%k_cell(:,i) = ctransf(:,i) * factor
        
      end do

      ! Avoid confusing permutations, revert to identity
      if ( expansion_factor == 1 ) then
        
        this%k_cell(:,:) = 0
        do i = 1, 3
          this%k_cell(i,i) = 1
        end do
        
      end if

      if ( present(process_k_cell) ) then
        call process_k_cell(this%k_cell, this%k_displ)
      end if

      ! Find the grid
      call find_kgrid(cell,this%k_cell, this%k_displ, .false., &
          this%trs, &
          this%N, this%k, this%w, this%cutoff)
      
    end subroutine setup_cutoff

    subroutine setup_file()

      call kpoint_read_file(file, cell, this)

    end subroutine setup_file

    subroutine setup_gamma()
      use alloc, only: re_alloc

      ! We have a Gamma-only calculation
      this%method = K_METHOD_NONE
      this%N = 1
      call re_alloc(this%k, 1, 3, 1, 1, name='kpoint%k')
      call re_alloc(this%w, 1, 1, name='kpoint%w')
      this%k = 0._dp
      this%w = 1._dp
      this%k_cell = 0
      do i = 1, 3
        this%k_cell(i,i) = 1
      end do
      
    end subroutine setup_gamma

  end subroutine kpoint_read

  ! The user can specify their own k-points
  subroutine kpoint_read_file(fname, cell, this)
    use alloc, only: re_alloc
    use parallel, only : Node
    use m_os, only : file_exist
#ifdef MPI
    use mpi_siesta
#endif
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: cell(3,3)
    type(kpoint_t), intent(inout) :: this

    real(dp) :: rcell(3,3), k(3), wsum

#ifdef MPI
    integer :: MPIerror
#endif

    ! The user has requested to read in 
    ! k-points from a specific file...
    ! We will do that for him/her.
    integer :: iu, ik, ix, stat

    if ( Node == 0 ) then

      if ( .not. file_exist(trim(fname)) ) then
        call die('Could not locate file '//trim(fname)// &
            ' please ensure that the file exists.')
      end if

      call io_assign( iu )
      open( iu, file=trim(fname), form='formatted', status='old' )

      ! Read number of k-points
      read(iu,*,iostat=stat) this%N
      call kill_iokp(stat,0)

    end if

#ifdef MPI
    call MPI_Bcast(this%N,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

    call re_alloc(this%k, 1, 3, 1, this%N, name='kpoint%k')
    call re_alloc(this%w, 1, this%N, name='kpoint%w')

    if ( Node == 0 ) then

      call reclat(cell, rcell, 1)

      ! Read in the k-points
      wsum = 0._dp
      do ik = 1 , this%N
        
        ! read current k-point
        read(iu,*,iostat=stat) ix, k(:), this%w(ik) ! (i6,3f12.6,3x,f12.6)
        call kpoint_convert(rcell, k, this%k(:,ik), -2)
        call kill_iokp(stat,ik)
        wsum = wsum + this%w(ik)

      end do

      if ( abs(wsum - 1._dp) > 1.e-7_dp ) then
        write(*,'(a)')'WARNING: Weights for user specified k-points does &
            &not sum to 1.'
        call die('User specified k-points does not sum to 1.')
      end if

      call io_close( iu )

    end if

#ifdef MPI
    call MPI_Bcast(this%k(1,1),3*this%N,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(this%w(1),this%N,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
#endif

  contains

    subroutine kill_iokp(stat, line)
      integer, intent(in) :: stat, line

      if ( stat == 0 ) return

      write(*,*) 'Siesta kpoint_read could not read your input file'
      write(*,*) 'The k-points MUST be in units of reciprocal vectors!'
      write(*,*) 'Siesta will convert the unit to correct units.'
      write(*,*) 'Also the sum of weights MUST equal 1.'
      write(*,*) !
      if ( line == 0 ) then
        write(*,*) 'Error occured on reading number of k-points (first line)'
      else
        write(*,'(a,i0,a)') 'Error occured on reading the ',line,' kpoint.'
      end if
      write(*,*) 'Please format your file like this:'
      write(*,*) ' $> cat '//trim(fname)
      write(*,*) ' <nkpt>'
      write(*,*) '     1  <kpt-A1> <kpt-A2> <kpt-A3> <w-kpt>'
      write(*,*) '     2  <kpt-A1> <kpt-A2> <kpt-A3> <w-kpt>'
      write(*,*) ' ....'
      write(*,*) ' <nkpt> <kpt-A1> <kpt-A2> <kpt-A3> <w-kpt>'

      call die('Siesta reading user specified k-points')

    end subroutine kill_iokp

  end subroutine kpoint_read_file


  !< Write k-point list to the XML output file
  !!
  !! This will populate a property list with all information about the
  !! k-point generation.
  subroutine kpoint_write_xml(this, prefix)
    use siesta_cml
    use m_char, only: lcase
    use units, only: Ang
    
    type(kpoint_t), intent(inout) :: this
    character(len=*), intent(in), optional :: prefix
    character(len=64) :: lcase_prefix
    integer :: ik

    ! Quick return, if able
    if ( .not. cml_p ) return

    if ( present(prefix) ) then
      lcase_prefix = lcase(prefix)
      
      call cmlStartPropertyList(xf=mainXML, title=trim(prefix)//".k-points", &
          dictRef="siesta:"//trim(lcase_prefix)//".kpoints")
      call cmlAddProperty(xf=mainXML, value=this%N, &
          dictref='siesta:'//trim(lcase_prefix)//'.nkpnt', &
          units="cmlUnits:countable")
      do ik = 1, this%N
        call cmlAddKPoint(xf=mainXML, coords=this%k(:,ik), weight=this%w(ik))
      end do

      ! Add supplementary information for MP and Cutoff methods
      select case ( this%method )
      case ( K_METHOD_MONKHORST_PACK, K_METHOD_CUTOFF )
        call cmlAddProperty(xf=mainXML, value=this%cutoff/Ang, & 
            dictref='siesta:'//trim(lcase_prefix)//'.kcutoff', units='siestaUnits:angstrom')
        call cmlAddProperty(xf=mainXML, value=this%k_cell, &
            dictref='siesta:'//trim(lcase_prefix)//'.kscell', units="cmlUnits:countable")
        call cmlAddProperty(xf=mainXML, value=this%k_displ, &
            dictref='siesta:'//trim(lcase_prefix)//'.kdispl')
      end select
           
    else
      call cmlStartPropertyList(xf=mainXML, title="k-points", &
          dictRef="siesta:kpoints")
      call cmlAddProperty(xf=mainXML, value=this%N, dictref='siesta:nkpnt', &
          units="cmlUnits:countable")
      do ik = 1, this%N
        call cmlAddKPoint(xf=mainXML, coords=this%k(:,ik), weight=this%w(ik))
      end do

      ! Add supplementary information for MP and Cutoff methods
      select case ( this%method )
      case ( K_METHOD_MONKHORST_PACK, K_METHOD_CUTOFF )
        call cmlAddProperty(xf=mainXML, value=this%cutoff/Ang, &
            dictref='siesta:kcutoff', units='siestaUnits:angstrom')
        call cmlAddProperty(xf=mainXML, value=this%k_cell, &
            dictref='siesta:kscell', units="cmlUnits:countable")
        call cmlAddProperty(xf=mainXML, value=this%k_displ, &
            dictref='siesta:kdispl')
      end select

    end if

    call cmlEndPropertyList(mainXML)

  end subroutine kpoint_write_xml
  
  !< Write to std-out the k-points and some information regarding the generation of the k-list
  !!
  !! The k-points are only written if `all` is `.true.`.
  !! Otherwise only information regarding the generation will be written.
  subroutine kpoint_write_stdout(this, all, prefix)
    use units, only: Ang
    type(kpoint_t), intent(in) :: this
    logical, intent(in) :: all
    character(len=*), intent(in), optional :: prefix
    character(len=64) :: name
    integer :: ik

    if ( present(prefix) ) then
      name = trim(prefix) // ' k-'
    else
      name = 'k-'
    end if

    if ( all ) then
      write(*,'(/,3a)') 'siesta: ',trim(name), 'point coordinates (Bohr**-1) and weights:'
      do ik = 1, this%N
        write(*,'(a,i10,3(tr1,e13.6),tr3,e12.6)') 'siesta: ', ik, this%k(:,ik), this%w(ik)
      end do
    end if
    write(*,'(/3a,i10)')  'siesta: ', trim(name), 'grid: Number of k-points =', this%N

    select case ( this%method )
    case ( K_METHOD_NONE )
      write(*,'(3a)')  'siesta: ', trim(name), 'point is Gamma-only'
    case ( K_METHOD_MONKHORST_PACK )
      write(*,'(3a)')  'siesta: ', trim(name), 'points from Monkhorst-Pack grid'
      write(*,'(3a,f10.3,a)')  'siesta: ', trim(name), 'cutoff (effective) =', this%cutoff/Ang, ' Ang'
      write(*,'(3a)') 'siesta: ', trim(name), 'point supercell and displacements'
      do ik = 1, 3
        write(*,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ', this%k_cell(:,ik), this%k_displ(ik)
      end do
    case ( K_METHOD_CUTOFF )
      write(*,'(3a)')  'siesta: ', trim(name), 'points from cutoff'
      write(*,'(3a,f10.3,a)')  'siesta: ', trim(name), 'cutoff (effective) =', this%cutoff/Ang, ' Ang'
      write(*,'(3a)') 'siesta: ', trim(name), 'point supercell and displacements'
      do ik = 1, 3
        write(*,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ', this%k_cell(:,ik), this%k_displ(ik)
      end do
    case ( K_METHOD_LIST )
      write(*,'(3a)')  'siesta: ', trim(name), 'points from user-defined list'
    end select
    
  end subroutine kpoint_write_stdout

  !< Write to a file the k-point information
  subroutine kpoint_write_file(this, suffix)
    use files, only: slabel
    use m_io, only: io_assign, io_close
    type(kpoint_t), intent(in) :: this
    character(len=*), intent(in) :: suffix
    character(len=256) :: fname
    integer :: iu, ik

    fname = trim(slabel) // '.' // trim(suffix)

    call io_assign( iu )
    open( iu, file=fname, form='formatted', status='unknown' )      
    
    write(iu,'(i10)') this%N
    do ik = 1, this%N
      write(iu,'(i10,3(tr1,e13.6),tr3,e12.6)') ik, this%k(:,ik), this%w(ik)
    end do
      
    call io_close( iu )
    
  end subroutine kpoint_write_file
  
  subroutine kpoint_fdf_name(prefix, suffix, name)
    character(len=*), intent(in) :: prefix, suffix
    character(len=256), intent(out) :: name

    if ( len_trim(prefix) > 0 ) then
      name = trim(prefix) // '.' // trim(suffix)
    else
      name = trim(suffix)
    end if

  end subroutine kpoint_fdf_name

end module kpoint_t_m





