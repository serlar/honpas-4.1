      module m_iodipol

!     This module implements a subroutine to write the dipole to
!     at every time step in case of TDDFT calculations.
!     Based on the modified version of Daniel Sanchez Portal's 
!     original 'iodipole' subroutine. 
!     Rafi Ullah November 2014.
!

      use m_dipol,   only: dipol
      use m_steps,   only: fincoor



      implicit none
     

      CONTAINS

      subroutine write_td_dipol(totime,istp,itd, ntd,rstart_time)


      integer, intent(in)           :: istp, itd, ntd
      double precision, intent(in)  :: totime, rstart_time
     
 
      logical, save         :: laststp = .false.
      logical, save         :: frstme  = .true.
      
      if (istp .gt. fincoor .and. itd .gt. ntd) then
         laststp = .true.
      end if

      call iodipole (totime, dipol, frstme, laststp, rstart_time)


      end subroutine write_td_dipol
!----------------------------------------------------------------
       
       subroutine  iodipole (totime, dipole,frstme, lastistp,rstart_time)
       
       use files,     only : slabel, label_length
       
       character(len=label_length+3) :: dipolefile
       external io_assign, io_close
       double precision dipole(3), extfield(3), totime, rstart_time
       integer iu
       character*15 fform
       logical, intent(inout) :: frstme
       logical, intent(in)    :: lastistp
       save iu 

      if(frstme) then
        dipolefile = trim(slabel)//'.dipol_vs_time'
        call io_assign( iu )
        fform='formatted'
        open( iu, file=dipolefile, form=fform, status='unknown' )
!        write(iu,'(a,3f15.6)') '#',extfield(1), extfield(2), extfield(3)
        frstme = .false.
      endif
       if (totime .gt. rstart_time) then
       write(iu,'(4f15.6)')                                             &
         totime,                                                        &
         dipole(1),dipole(2), dipole(3)
       end if
      if(lastistp) call io_close(iu)

      end subroutine iodipole

      end module m_iodipol
