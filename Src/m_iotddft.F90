      MODULE m_iotddft

!     This module has different subroutines to write the total 
!     E_KS, the instantaneous so-called Eigen values and dipole moment
!     at every time step and time-dependent density (rho) after every 
!     given number of steps in case of TDDFT calculations depending 
!     on the user's choice.
!     Based on the modified version of Daniel Sanchez Portal's 
!     original  subroutines. 
!     Rafi Ullah November 2014.
!

      USE m_dipol,          ONLY: dipol
      USE m_steps,          ONLY: fincoor, istp
      USE files,            ONLY: slabel, label_length
      USE siesta_options,   ONLY: eigen_time, dip_time, etot_time, tdsaverho, &
                                  tdsaverho, ntdsaverho
      USE wavefunctions,    ONLY: wavef_ms
      USE parallel,         ONLY: IONode
      USE files,            ONLY: filesOut_t
      USE units,            ONLY: eV
      USE m_io,             ONLY: io_assign, io_close
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: write_tddft
      PUBLIC :: write_tdrho
      CHARACTER(LEN=15)  :: fform

      CONTAINS

      SUBROUTINE write_tddft(totime,istp,itd, ntd,rstart_time,         &
                             etot,eigen, maxo, nspin, nk)

      INTEGER, INTENT(IN)           :: istp, itd, ntd, maxo, nk, nspin
      DOUBLE PRECISION, INTENT(IN)  :: totime, rstart_time, etot
      DOUBLE PRECISION, INTENT(IN)  :: eigen(maxo,nspin,nk) 
      LOGICAL, SAVE                 :: laststp = .false.
      
      IF (istp .ge. fincoor .AND. itd .ge. ntd) THEN
         laststp = .true.
      END IF
      
      IF (dip_time) THEN
      CALL iodipole (totime, dipol, laststp, rstart_time)
      END IF 
      
      IF (etot_time) THEN
      CALL ioetot   (totime, etot, laststp, rstart_time)
      END IF
      
      IF (eigen_time) THEN
      CALL ioeigenvalues (totime, eigen, laststp, rstart_time, &
                                 maxo, nspin, nk)
      END IF


      END SUBROUTINE write_tddft
!----------------------------------------------------------------
      SUBROUTINE write_tdrho (filesOut)
      
      TYPE(filesOut_t), INTENT(INOUT)      :: filesOut
  
      IF (tdsaverho) THEN
         IF (mod(istp,ntdsaverho) .eq. 0) THEN
            write(filesOut%tdrho,"(i0,a)") istp, '.TDRho'
        ELSE
          filesOut%tdrho = ' '
        END IF
      END IF
      
      END SUBROUTINE write_tdrho

      SUBROUTINE  iodipole (totime, dipole,lastistp,rstart_time)
       
       
       CHARACTER(LEN=70 ) :: dipolefile
       DOUBLE PRECISION   :: dipole(3), extfield(3), totime, rstart_time
       INTEGER, SAVE      :: iu
       LOGICAL,    SAVE      :: frstme  = .true.
       LOGICAL, INTENT(IN)  :: lastistp

!      Only first node writes
       IF(IONode) THEN

       IF (frstme) THEN
         dipolefile = trim(slabel) // '.TDDIPOL'
         call io_assign( iu )
         fform='formatted'
         OPEN( iu, FILE=dipolefile, FORM=fform,POSITION='APPEND', STATUS='REPLACE' )
!        write(iu,'(a,3f15.6)') '#',extfield(1), extfield(2), extfield(3)
         frstme = .false.
       END IF
          WRITE(iu,'(4f15.6)')                                         &
          totime,                                                      &
          dipole(1),dipole(2), dipole(3)
       IF (lastistp) call io_close(iu)
        
       END IF ! IONode
      END SUBROUTINE iodipole
!----------------------------------------------------------------------
      SUBROUTINE ioetot (totime, etot, lastistp, rstart_time)
        
       DOUBLE PRECISION         ::totime, etot, rstart_time
       INTEGER                  :: iu, istp, itd, ntd
       LOGICAL                  :: lastistp
       LOGICAL, SAVE            :: frstme = .true. 
       SAVE                     :: iu
       CHARACTER(LEN=70)        :: etotfile
        
        IF(IONode) THEN

        IF (frstme) THEN
          etotfile = trim(slabel) // '.TDETOT'
          CALL io_assign ( iu )
          fform = 'formatted'
          OPEN (iu, FILE=etotfile, FORM=fform,POSITION='APPEND',STATUS='REPLACE')
          frstme = .false.
        END IF
           WRITE (iu, '(2f15.6)') totime, etot/eV
        IF (lastistp) CALL io_close(iu)
        END IF ! IONode
      END SUBROUTINE ioetot
!------------------------------------------------------------------------

SUBROUTINE ioeigenvalues (totime, eigen, lastistp, rstart_time, &
                           maxo, nspin, nk)

 INTEGER            :: maxo, nspin, nk, ik, ispin, ie
 INTEGER            :: nocc(nk,nspin)
 INTEGER, SAVE      :: iuu
 DOUBLE PRECISION   :: totime, rstart_time
 DOUBLE PRECISION   :: eigen(maxo,nspin,nk)
 LOGICAL            :: lastistp
 LOGICAL, SAVE      :: frstme = .true.
 CHARACTER(LEN=70)  :: eigenfile

 IF (IONode) THEN 
   IF (frstme) THEN
     eigenfile = trim(slabel) // '.TDEIG'
     fform = 'formatted'
     call io_newunit(iuu)
     OPEN (UNIT=iuu, FILE=eigenfile, FORM=fform, POSITION='APPEND',      &
          STATUS='REPLACE')
     WRITE(iuu,*) '#  ', nspin, nk
     frstme = .false.
   END IF
     WRITE(iuu,"(f12.8,/)") totime
     DO ik = 1, nk
       WRITE(iuu,"(i5,10f12.5,/,(5x,10f12.5))")               &
            ik, ((eigen(ie,ispin,ik)/eV,ie=1,(wavef_ms(ik,ispin)%dim2)),      &
            ispin=1,nspin)
     END DO
   IF (lastistp) CLOSE (iuu)
  END IF ! IONode
END SUBROUTINE ioeigenvalues 

  subroutine io_newunit(lun)
      integer, intent(out) :: lun
      logical used
      integer iostat

!     Looks for a free unit and assigns it to lun

      do lun= 10, 99
            inquire(unit=lun, opened=used, iostat=iostat)
            if (iostat .ne. 0) used = .true.
            if (.not. used) return
      enddo
  end subroutine io_newunit



       END MODULE m_iotddft
