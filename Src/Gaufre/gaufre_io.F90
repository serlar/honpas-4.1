! *** Module: gaufre_io ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Collection of I/O routines for the fitting of Natural Atomic
!!        Orbitals to Gaussians
!!
!! This module reads and writes data and parameters related to the fitting
!! of Natural Atomic Orbitals by Gaussians.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 04.2018 Created [Yann Pouillon]
! *****************************************************************************
module gaufre_io

  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

  use gaufre_common

  implicit none

  private

  public :: &
&   gaufre_io_free, &
&   gaufre_io_init, &
&   gaufre_io_load_data, &
&   gaufre_io_save_data, &
&   gaufre_io_t

  ! I/O descriptor
  type :: gaufre_io_t
    integer :: errno = 0
    integer :: fd = -1
    integer :: fmt = GAUFRE_DATA_FMT_UNKNOWN
    integer :: mode = GAUFRE_IO_UNKNOWN
  end type

  ! Overloaded interfaces, to simplify the namespace
  interface gaufre_io_load_data
    module procedure :: gaufre_io_load_data_columns
  end interface
  interface gaufre_io_save_data
    module procedure :: gaufre_io_save_data_columns, gaufre_io_save_data_yaml
  end interface

contains

  ! ***************************************************************************
  !> \brief Destructor for the I/O descriptor
  !!
  !! This routine frees all memory allocated for the Gaussian fitting procedure
  !! configuration and should be called once the data structure is not needed
  !! anymore.
  !!
  !! \note This routine is idempotent and can be run as many times as desired.
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[in] iodesc: data structure containing I/O parameters
  ! ***************************************************************************
  subroutine gaufre_io_free(iodesc)

    implicit none

    ! Arguments
    type(gaufre_io_t), intent(inout) :: iodesc

    ! Local variables
    logical :: do_close

    ! -------------------------------------------------------------------------

    inquire(unit=iodesc%fd, opened=do_close)
    if ( do_close ) then
      close(unit=iodesc%fd)
    end if

    iodesc%fd = -1
    iodesc%fmt = GAUFRE_DATA_FMT_UNKNOWN
    iodesc%mode = GAUFRE_IO_UNKNOWN

  end subroutine gaufre_io_free

  ! ***************************************************************************
  !> \brief Constructor for the Gaussian expansion configuration
  !!
  !! This routine allocates all the necessary memory and initializes the
  !! corresponding parameters before running a Gaussian fitting procedure.
  !! It must be called before using the specified data structure.
  !!
  !! \note We want to be able to use subsets of mesh points and orbital values
  !!       not necessarily including boundaries and without interfering with
  !!       the original variables, this is why we explictly copy them for now.
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[out] iodesc: data structure containing I/O parameters
  !! \param[in] fd: file descriptor (Fortran unit) to use
  !! \param[in] iomode: GAUFRE_IO_LOAD or GAUFRE_IO_SAVE
  ! ***************************************************************************
  subroutine gaufre_io_init(iodesc, fd, iomode)

    implicit none

    ! Arguments
    type(gaufre_io_t), intent(out) :: iodesc
    integer, intent(in) :: fd
    integer, intent(in) :: iomode

    ! Local variables
    character(len=32) :: inq_access, inq_action, inq_formstr, run_action
    logical :: inq_opened, inq_formatted
    integer :: inq_number, inq_mode

    ! -------------------------------------------------------------------------

    ! Translate arguments and check their consistency
    select case ( iomode )
      case(GAUFRE_IO_LOAD)
        run_action = "READ"
      case(GAUFRE_IO_SAVE)
        run_action = "WRITE"
      case default
        write(unit=error_unit, fmt=*) "Error: unknown I/O mode ", iomode
        iodesc%errno = -1
        return
    end select

    ! Gather information about the file
    inquire(unit=fd, access=inq_access, action=inq_action, &
&     formatted=inq_formstr, number=inq_number, opened=inq_opened)
    if ( trim(inq_action) == "READ" ) then
      inq_mode = GAUFRE_IO_LOAD
    else if ( trim(inq_action) == "WRITE" ) then
      inq_mode = GAUFRE_IO_SAVE
    else
      inq_mode = iomode
    end if
    if ( trim(inq_formstr) == "YES" ) then
      inq_formatted = .true.
    else
      inq_formatted = .false.
    end if

    ! Perform necessary I/O operations
    if ( (.not. inq_opened) .or. (inq_number == -1) ) then
      open(unit=fd, action=trim(run_action), form="formatted", &
&       iostat=iodesc%errno)
      if ( iodesc%errno /= 0 ) return
    else
      if ( trim(inq_access) /= "SEQUENTIAL" ) then
        write(unit=error_unit, fmt=*) &
&         "Error: Fortran unit has been opened with incompatibe I/O access type"
        iodesc%errno = -2
        return
      end if
      if ( .not. inq_formatted ) then
        write(unit=error_unit, fmt=*) &
&         "Error: Fortran unit has been opened with incompatibe I/O format"
        iodesc%errno = -3
        return
      end if
      if ( inq_mode /= iomode ) then
        write(unit=error_unit, fmt=*) &
&         "Error: Fortran unit has been opened with incompatibe I/O mode"
        iodesc%errno = -4
        return
      end if
      rewind(unit=fd, iostat=iodesc%errno)
      if ( iodesc%errno /= 0 ) return
    end if

    ! Set I/O parameters
    iodesc%fd = fd
    iodesc%mode = iomode

  end subroutine gaufre_io_init

  ! ***************************************************************************
  !> \brief Loads a (X,F(X)) function from a simple text file
  !!
  !! This routine loads a (X,F(X)) pair written in columns in a text file,
  !! with a possible initial offset to ignore headers.
  !!
  !! It checks that X is a valid mesh of points (i.e. strictly monotone) and
  !! that the values F(X) are not set to zero beyond a certain point.
  !!
  !! \note If the destination arrays are already allocated, their sizes must
  !!       be equal and npts must be both strictly positive and less than or
  !!       equal to them.
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[inout] iodesc: data structure containing I/O parameters
  !! \param[inout] npts: number of points to read / actually read
  !! \param[out] xdata: array pointer of dimension 1
  !! \param[out] ydata: array pointer of dimension 1
  !! \param[in] io_fmt: file format to use (see \ref gaufre_common)
  !! \param[out] orb_data: optional data structure to store orbital information
  ! ***************************************************************************
  subroutine gaufre_io_load_data_columns(iodesc, npts, xdata, ydata, io_fmt, &
&   orb_data)

    use gaufre_orbital

    implicit none

    ! Arguments
    type(gaufre_io_t), intent(inout) :: iodesc
    integer, intent(inout) :: npts
    real(gfdp), dimension(:), pointer, intent(out) :: xdata, ydata
    integer, intent(in) :: io_fmt
    type(gaufre_orbital_t), optional, intent(out) :: orb_data

    ! Local variables
    character(len=3) :: orb_lbl
    character(len=80) :: garbage
    integer :: ioerr, ipts, maxpts, orb_l, orb_n, orb_z
    logical :: do_read, orb_pol, xalloc, yalloc
    real(gfdp) :: orb_pop, tol, xprec, xread, xdir, yprec, yread

    ! -------------------------------------------------------------------------

    ! Init
    do_read = .true.
    iodesc%errno = 0
    ioerr = 0
    tol = epsilon(1.0_gfdp)**2
    xprec = 0.0_gfdp

    ! Check the file format
    if ( io_fmt /= GAUFRE_IO_FMT_SIESTA ) then
      write(unit=error_unit, fmt='("Error: ",A,":",I2)') &
&       "File format not implemented", io_fmt
      return
    end if

    ! Check that the I/O descriptor allows data reading
    if ( (iodesc%fd == -1) .or. (iodesc%mode /= GAUFRE_IO_LOAD) ) then
      write(unit=error_unit, fmt=*) "Error: I/O descriptor not opened for reading"
      return
    end if

    ! Check that the number of points is positive
    if ( npts < 0 ) then
      write(unit=error_unit, fmt=*) "Error: invalid number of points (", &
&       npts, ")"
      return
    end if

    ! Manage destination arrays
    xalloc = associated(xdata)
    yalloc = associated(ydata)
    maxpts = 0
    if ( xalloc .and. yalloc ) then
      if ( size(xdata, 1) /= size(ydata, 1) ) then
        write(unit=error_unit, fmt=*) &
&         "Error: XDATA and YDATA must be in the same state"
        return
      end if
      maxpts = size(xdata, 1)
      if ( npts > 0 ) then
        maxpts = min(maxpts, npts)
      end if
    else if ( (.not. xalloc) .and. (.not. yalloc) ) then
      if ( npts == 0 ) then
        write(unit=error_unit, fmt=*) &
&         "Error: desired number of points must be strictly positive (got ", &
&         npts, ")"
        return
      end if
      maxpts = npts
      allocate(xdata(maxpts))
      allocate(ydata(maxpts))
      xdata(:) = 0.0_gfdp
      ydata(:) = 0.0_gfdp
    else
      write(unit=error_unit, fmt=*) &
&       "Error: XDATA and YDATA must be in the same state"
      write(unit=error_unit, fmt=*) "(XALLOC is ", xalloc, ", &
&       YALLOC is ", yalloc, ")"
      return
    end if

    ! Read header
    if ( io_fmt == GAUFRE_IO_FMT_SIESTA ) then
      read(unit=iodesc%fd, iostat=iodesc%errno, &
&       fmt='(A1,A3,4(1X,I2),1X,F9.4)') &
&       garbage, orb_lbl, orb_l, orb_n, orb_z, ipts, orb_pop
      if ( iodesc%errno /= 0 ) return

      if ( ipts == 0 ) then
        orb_pol = .false.
      else
        orb_pol = .true.
      end if

      read(unit=iodesc%fd, iostat=iodesc%errno, fmt='(A)') garbage
      if ( iodesc%errno /= 0 ) return

      if ( present(orb_data) ) then
        orb_data%species = adjustl(orb_lbl)
        orb_data%qn_n = orb_n
        orb_data%qn_l = orb_l
        orb_data%zeta = orb_z
        orb_data%polarized = orb_pol
        orb_data%population = orb_pop
      end if
    end if

    ! Read whatever data is requested and available
    ipts = 0
    do while ( do_read .and. (ioerr == 0) )
      read(unit=iodesc%fd, iostat=ioerr, fmt=*) xread, yread
      if ( ioerr /= 0 ) exit

      ! Determine the direction of X (forward or backward)
      if ( ipts == 1 ) then
        if ( (xread - xprec) > tol ) then
          xdir = 1.0_gfdp
        else
          xdir = -1.0_gfdp
        end if
      end if

      ! Enforce strictly monotonic X meshes
      if ( ipts > 0 ) then
        if ( (xdir * (xread - xprec)) < tol ) then
          write(unit=error_unit, fmt=*) "Error: mesh is not strictly monotonic"
          exit
        end if
      end if

      ! Stop if F(X) is too abruptly cut (better quality of the fitting)
      if ( ipts > 0 ) then
        if ( (abs(yread) < tol**2) .and. (abs(yprec) > abs(yread*1.0e3_gfdp)) ) then
          exit
        end if
      end if

      ipts = ipts + 1
      xdata(ipts) = xread
      ydata(ipts) = yread
      xprec = xread
      yprec = yread
      if ( ipts >= maxpts ) do_read = .false.
    end do

    ! Transmit final number of points
    npts = ipts

    ! Transmit errors other than end-of-file
    if ( ioerr > 0 ) iodesc%errno = ioerr

  end subroutine gaufre_io_load_data_columns

  ! ***************************************************************************
  !> \brief Saves a (X,F(X)) function to a simple text file
  !!
  !! This routine saves a (X,F(X)) pair of variables to columns in a
  !! text file.
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[inout] iodesc: data structure containing I/O parameters
  !! \param[in] xdata: X mesh points
  !! \param[out] yfunc: procedure pointer to a function calculating F(X)
  ! ***************************************************************************
  subroutine gaufre_io_save_data_columns(iodesc, xdata, yfunc)

    implicit none

    ! Arguments
    type(gaufre_io_t), intent(inout) :: iodesc
    real(gfdp), dimension(:), pointer, intent(in) :: xdata

    ! Interfaces of procedure arguments
    interface
      function yfunc(x)
        use gaufre_common, only: gfdp
        implicit none
        real(gfdp), intent(in) :: x
        real(gfdp) :: yfunc
      end function yfunc
    end interface

    ! Local variables
    integer :: ioerr, ipts

    ! -------------------------------------------------------------------------

    ! Check that the I/O descriptor allows data writing
    if ( (iodesc%fd == -1) .or. (iodesc%mode /= GAUFRE_IO_SAVE) ) then
      write(unit=error_unit, fmt=*) &
&       "Error: I/O descriptor not opened for writing"
      return
    end if

    ! Save the values of X and F(X) in columns
    do ipts = 1,size(xdata, 1)
      write(unit=iodesc%fd, fmt='(2(1X,E20.12))', iostat=ioerr) &
&       xdata, yfunc(xdata(ipts))
      if ( ioerr /= 0 ) exit
    end do

    ! Transmit errors
    if ( ioerr /= 0 ) iodesc%errno = ioerr

  end subroutine gaufre_io_save_data_columns

  ! ***************************************************************************
  !> \brief Saves a (X,F(X)) function to a file in YAML format
  !!
  !! This routine saves a (X,F(X)) pair of variables to columns in a
  !! text file in YAML format.
  !!
  !! \todo Allow indenting.
  !!
  !! \par History
  !!      - 04.2018 Created [Yann Pouillon]
  !!
  !! \param[inout] iodesc: data structure containing I/O parameters
  !! \param[in] varname: name of the YAML variable
  !! \param[in] xdata: X mesh points
  !! \param[out] yfunc: procedure pointer to a function calculating F(X)
  ! ***************************************************************************
  subroutine gaufre_io_save_data_yaml(iodesc, varname, xdata, yfunc)

    implicit none

    ! Arguments
    type(gaufre_io_t), intent(inout) :: iodesc
    character(len=*), intent(in) :: varname
    real(gfdp), dimension(:), pointer, intent(in) :: xdata

    ! Interfaces of procedure arguments
    interface
      function yfunc(x)
        use gaufre_common, only: gfdp
        implicit none
        real(gfdp), intent(in) :: x
        real(gfdp) :: yfunc
      end function yfunc
    end interface

    ! Local variables
    integer :: ioerr, ipts

    ! -------------------------------------------------------------------------

    ! Check that the I/O descriptor allows data writing
    if ( (iodesc%fd == -1) .or. (iodesc%mode /= GAUFRE_IO_SAVE) ) then
      write(unit=error_unit, fmt=*) "Error: I/O descriptor not opened for writing"
      return
    end if

    ! Save the values of X and F(X) as a YAML list
    write(unit=iodesc%fd, fmt='(A,":")', iostat=ioerr) varname
    do ipts = 1,size(xdata, 1)
      if ( ioerr /= 0 ) exit
      write(unit=iodesc%fd, fmt='(2X,"- [",E20.12,",",E20.12,"]")', &
&       iostat=ioerr) xdata, yfunc(xdata(ipts))
    end do
    write(unit=iodesc%fd, fmt='(A)', iostat=ioerr) ""

    ! Transmit errors
    if ( ioerr /= 0 ) iodesc%errno = ioerr

  end subroutine gaufre_io_save_data_yaml

end module gaufre_io
