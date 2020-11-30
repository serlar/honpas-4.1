! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_siesta_end
      private
      public :: siesta_end

      CONTAINS

      subroutine siesta_end( )
      use fdf,             only : fdf_shutdown
      use alloc,           only : alloc_report
      use siesta_cml,      only : cml_p, mainXML, cmlEndModule
      use siesta_cmlsubs,  only : siesta_cml_exit
      use siesta_master,   only : siesta_subroutine  ! Is siesta a subroutine?
#ifdef DEBUG_XC
      USE siestaXC,        only : closeDebugOutputFile
#endif /* DEBUG_XC */
      use m_timestamp,     only : timestamp
      use m_wallclock,     only : wallclock, wallclock_shutdown
      use parallel,        only : IOnode, SIESTA_worker, Node
#ifdef MPI
      use mpi_siesta
#endif
      use moreMeshSubs,    only : resetMeshDistr
      use densematrix,     only : resetDenseMatrix
      use sparse_matrices, only : resetSparseMatrices
      use m_fixed,         only : resetFixedPointers
      use band,            only : reset_bands
      use meshphi,         only : resetMeshPhi
      use meshdscf,        only : resetDscfPointers
      use m_dscfcomm,      only : resetdscfComm
      use m_rhog,          only : resetRhog
      ! NP, new mixing
      use m_mixing_scf,    only : mixers_scf_reset
      USE siesta_options,  only : isolve, SOLVE_CHESS

      use kpoint_scf_m, only: reset_kpoint_scf
      use kpoint_dos_m, only: reset_kpoint_dos
      use kpoint_pdos_m, only: reset_kpoint_pdos
      use kpoint_ldos_m, only: reset_kpoint_ldos
      use ts_kpoint_scf_m, only: reset_ts_kpoint_scf

#ifdef SIESTA__CHESS
      use m_chess, only: CheSS_finalize
#endif

#ifdef SIESTA__FLOOK
      use flook_siesta,    only : slua_close
      use siesta_options,  only : LUA
#endif

      use m_diag, only: diag_exit

#ifdef TRACING
      use extrae_eventllist, only : writeList, deleteList, eventlist
#endif
      use m_io,            only : io_assign, io_close

      use nao2gto_data, only: hfx_libint, hfx_system
      use nao2gto_data, only: nao2gto_libint_end, nao2gto_system_end
      use siesta_options,  only : hfx_wanted

      implicit none

      integer :: iu
#ifdef MPI
      integer MPIerror
#endif
!------------------------------------------------------------------------- BEGIN
      if (SIESTA_worker) then

#ifdef SIESTA__FLOOK
      ! Clean up the lua environment
      call slua_close(LUA)
#endif

!     Free allocated memory
      call resetSparseMatrices( )
      call resetMeshDistr( )
      call resetdscfComm( )
      call resetDenseMatrix( )
      call resetFixedPointers( )
      call reset_bands( )
      call resetMeshPhi( )
      call resetDscfPointers( )
      call resetRhog()
      call mixers_scf_reset()

      call reset_kpoint_ldos()
      call reset_kpoint_pdos()
      ! This *must* be called after LDOS+PDOS since
      ! they may be associated with DOS
      call reset_kpoint_dos()
      ! This *must* be called before SCF since
      ! they may be associated with SCF
      call reset_ts_kpoint_scf()
      call reset_kpoint_scf()

      ! Clean diagonalization
      call diag_exit()

#ifdef SIESTA__CHESS
      if (isolve == SOLVE_CHESS) then
          call CheSS_finalize()
      end if
#endif
      
!     Output memory use up to the end of the program
!!      call printmemory( 6, 1 )

!     Print allocation report
      call alloc_report( printNow=.true. , shutdown=.true.)

#ifdef DEBUG_XC
! Close debug outup file
      call closeDebugOutputFile()
#endif /* DEBUG_XC */

!     Stop time counter
      call timer( 'siesta', 2 )
      call timer( 'all', 3 )

#ifdef TRACING
!     Write extrae user function list
      if (IOnode) then
        call writeList(eventlist)
      end if
      call deleteList(eventlist)
#endif /* TRACING */


!     Print final date and time
      if (IOnode) then
        write(*,*) ! Newline
        call timestamp('End of run')
        call wallclock('End of run')
        ! Ensure clock file is closed
        call wallclock_shutdown()
      endif

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif

!     Free hybrid XC data structures
!     FIXME: commented for now because Libint causes a SIGABRT
!      if ( hfx_wanted ) then
!        call nao2gto_libint_end(hfx_libint)
!        call nao2gto_system_end(hfx_system)
!      endif

      endif  ! SIESTA_worker

! Finalize MPI, unless siesta is running as a subroutine
      !
      ! Add stamp to MESSAGES file
      call message('INFO','Job completed')

      ! Clean the fdf-files
      call fdf_shutdown()

#ifdef MPI
      call MPI_Barrier(MPI_Comm_World, MPIerror)
      MPIerror = MPI_SUCCESS
      if (.not.siesta_subroutine) call MPI_Finalize( MPIerror )
#endif

      if (IOnode) then
         call io_assign(iu)
         open(iu, file='0_NORMAL_EXIT', form='FORMATTED')
#ifdef MPI
         write(iu,'(a,i0)') 'SIESTA completed, MPI exit: ',
     &        MPIerror
#else
         write(iu,'(a)') 'SIESTA completed successfully'
#endif
         call io_close(iu)
      endif
!--------------------------------------------------------------------------- END
      END subroutine siesta_end
      END MODULE m_siesta_end
