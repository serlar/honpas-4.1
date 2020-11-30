program test_gaussian_fit

    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

    use nao2gto_fit

    implicit none

    ! Constants
    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: NGFS = 8
    integer, parameter :: NPRM = 3
    integer, parameter :: NPTS = 500

    ! Variable declarations
    logical :: converged
    integer :: idat, igfs, ios
    real(DP) :: rfit
    integer, dimension(NGFS) :: best_fits
    real(DP), dimension(NPRM,NGFS,NGFS) :: coeffs
    real(DP), dimension(NPTS) :: xdata
    real(DP), dimension(NPTS) :: ydata
    real(DP), dimension(NPTS,NGFS) :: ygauss
    real(DP), dimension(NGFS) :: residuals

    open(unit=10, file="test_gaussian_fit.dat", status="old", &
&       action="read", iostat=ios)
    do idat=1,NPTS
        if ( ios /= 0 ) exit
        read(unit=10, fmt=*, iostat=ios) xdata(idat), ydata(idat)
    end do
    close(unit=10)

    if ( ios /= 0 ) then
        write(unit=error_unit, fmt=*) "WARNING: problem reading input data"
    end if

    best_fits(:) = 0
    coeffs(:,:,:) = 0.0_DP
    residuals(:) = 0.0_DP
    converged = .true.
    do igfs=1,NGFS
        converged = converged .and. &
&           nao2gto_fit_find(igfs, NPTS, xdata, ydata, coeffs(:,1:igfs,igfs), &
&               residuals(igfs))
        call nao2gto_fit_eval(igfs, NPTS, coeffs(:,1:igfs,igfs), xdata, &
&           ygauss(:,igfs))
    end do

    if ( converged ) then
        rfit = huge(0.0_DP)
        idat = 1
        do igfs=1,NGFS
            if ( residuals(igfs) < rfit ) then
                best_fits(idat) = igfs
                idat = idat + 1
                rfit = residuals(igfs)
            end if
        end do
        write(unit=output_unit, fmt=*)
        write(unit=output_unit, fmt=*) "BEST FITS OBTAINED WITH:"
        write(unit=output_unit, fmt=*)
        do igfs=1,3
            if ( idat - igfs > 0 ) then
                write(unit=output_unit, &
&                   fmt='(3X,"- GAUSSIANS: ",I4,1X,"(RESIDUAL: ",E12.5,")")') &
&                       best_fits(idat-igfs), residuals(best_fits(idat-igfs))
            end if
        end do
        write(unit=output_unit, fmt=*)

        open(unit=11, file="test_gaussian_fit.out", status="replace", &
&           action="write", iostat=ios)
        if ( ios == 0 ) then
            do idat=1,NPTS
                if ( ios /= 0 ) exit
                write(unit=11, fmt='(10(1X,E20.8))') xdata(idat), ydata(idat), &
&                   ygauss(idat,1:min(8,NGFS))
            end do
            close(unit=11)
        else
            write(unit=error_unit, fmt=*) &
&               "WARNING: problem writing output data"
        end if

        open(unit=12, file="test_gaussian_fit.chk", status="replace", &
&           action="write", iostat=ios)
        if ( ios == 0 ) then
            do igfs=1,NGFS
                if ( ios /= 0 ) exit
                write(unit=12, fmt='(I4,1X,E20.8)') igfs, residuals(igfs)
            end do
            close(unit=12)
        else
            write(unit=error_unit, fmt=*) &
&               "WARNING: problem writing residual data"
        end if
    else
        write(unit=error_unit, fmt=*) &
&           "WARNING: the minimization did not converge"
    end if

end program test_gaussian_fit
