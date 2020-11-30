BSC changes in SIESTA
=====================

Summary
-------

The Barcelona Computing Center has contributed code to SIESTA in the form of a
patch bomb which has been merged into the trunk at the beginning of 2011:

.. code::

   2011-02-03  11:15:00 GMT Alberto Garcia <albertog@icmab.es>   trunk-368
   Incorporation of BSC code: mesh optimization; parallel fdf library

The ongoing collaboration with the Barcelona Supercomputer Center
(through the department of Computer Applications in Science and
Engineering) has resulted in two new main enhancements:

- New code to improve the load-balance among processors in the
  operations carried out in the real-space mesh (involved in the
  computation of charge density, exchange-correlation potential, and
  matrix elements of the hamiltonian, among others).
- A new fdf library which avoids the 'master node bottleneck': each
  MPI processor now holds a copy of the relevant information, and can
  satisfy fdf queries without the need for interprocess communication or
  file access.


Load-balancing improvements
---------------------------

Since the different operations have different work-load patterns, it
has been necessary to define and handle three separate data
distributions, with appropriately optimized communication patterns.
Most of the new code is contained in new specialized modules, and
there are calls for data distribution handling in higher-level
modules, notably dhscf.

An overview of the rationale and implementation of the new functionality
can be found in:

    An efficient implementation of a QM-MM method in SIESTA,
    C. F. Sanz-Navarro, R. Grima, A. Garcia, E. A. Bea, A. Soba,
    J.M. Cela, and P. Ordejon, 
    Theoretical Chemistry Accounts. Online pub Sep 2010. 
    DOI: 10.1007/s00214-010-0816-5

We are working on a set of developer-oriented notes to facilitate
future enhancements.


Parallel FDF library
--------------------

The most noticeable change is in the processing of fdf blocks. The previous
version of the library provided a Fortran unit number from which the
user code could read the block contents. Now the user is given an
abstract handle with which to perform explicitly the operations of reading
lines, parse, and extract useful information. This paradigm was already
present as a wrapper in previous versions of Siesta, but it is now mandatory.

The two new developments, while formally independent, have been merged
at the same time due to historical and logistical reasons. Even though
the mesh enhancements affect directly only a relatively small number
of files, many others have been modified: some by changes in memory
allocation syntax and indentation, and some to clarify code patterns
or interfaces (see below). The new treatment of fdf blocks and the
removal of the explicit broadcast step in some existing fdf calls has
meant changes to many other files. In all, this patch affects several
dozen files in the Siesta distribution.

Other significant changes
-------------------------

The interface to the spherical-harmonics module in matel has been
changed to avoid the creation of temporary arrays at run time. Also,
some internal arrays in the spherical-harmonics module are now
pre-allocated to avoid overhead.

Throughout the code most 'allocatable' arrays have been turned into
pointers, to be able to exploit the functionality in the alloc module
to resize them concisely. It can be argued that we have gone too far
in some cases.

We have invested a lot of effort in this work, but it is very likely
that some issues have escaped us. A list of those known at this
time can be found in the file Docs/KNOWN.ISSUES.BSC-patch


Details of the commit
---------------------

Here is the full list of affected files:

.. code::

   removed:
     Src/FoX/arch.make
     Src/doping.F
     Src/fdf/fdf.f
     Src/fdf/fdf.h
     Src/fdf/fdf_mod.f
     Src/fdf/fdfdefs.h
     Src/fdf/io_for_fdf_tests.f
     Src/fdf/io_sample.f
     Src/fdf/parse.f
     Src/fdf/sample.f
     Src/m_efield.F
     Tests/bsc-compare.sh
     Tests/std-compare.sh
     Util/Denchar/Src/fdf/.dummy_dir
     Util/Gen-basis/fdf/.dummy_directory
     Util/TBTrans/Libs/
     Util/TBTrans/Libs/.dummy_directory
     Util/TBTrans/MPI/
     Util/TBTrans/MPI/.dummy_dir
     Util/TBTrans/fdf/
     Util/VCA/fdf/.dummy_directory
     Util/Vibra/Src/fdf/.dummy_directory
   added:
     .bzrignore
     Docs/BSC.CHANGES
     Docs/KNOWN.ISSUES.BSC-patch
     Src/FoX/arch.make
     Src/Libs/dgetrs_lapack.f
     Src/Libs/dspev_lapack.f
     Src/Sys/altix-32b-par.make
     Src/Sys/altix-32b-ser.make
     Src/Sys/atto-pgf95-openmpi.make
     Src/Sys/finisterrae.make.moved
     Src/Sys/g95-macosx-netcdf.make
     Src/Sys/gfortran-macosx-netcdf.make
     Src/Sys/intel-mpi-checks-metis.make
     Src/Sys/intel10-openmpi.make
     Src/Sys/intel11-openmpi.make
     Src/Sys/macosx-openmpi.make
     Src/Sys/marenostrum-lanczos.make
     Src/Sys/marenostrum-mpi-metis-64.make
     Src/Sys/mn-32b-par.make
     Src/Sys/mn-32b-ser.make
     Src/Sys/mn-openmp.make
     Src/Sys/nano-intel10-mvapich.make
     Src/Sys/pgf95-pgimpi.make
     Src/Sys/rogeli-intel-par.make
     Src/bsc_cellxc.F
     Src/bsc_xcmod.F
     Src/cellxc_mod.F
     Src/debugmpi.F
     Src/doping_uniform.F
     Src/fdf.Standard
     Src/fdf/XY.fdf
     Src/fdf/arch.make
     Src/fdf/fdf.F90
     Src/fdf/hostfile
     Src/fdf/io_fdf.F90
     Src/fdf/iso_fortran_env.F90
     Src/fdf/parse.F90
     Src/fdf/prec.F90
     Src/fdf/sample.F90
     Src/fdf/submit.sh
     Src/fdf/tags
     Src/fdf/utils.F90
     Src/final_H_f_stress.F
     Src/m_diagon.F
     Src/m_dscfcomm.F
     Src/m_efield.F
     Src/mesh.F
     Src/meshcomm.F
     Src/meshphi.F
     Src/moremeshsubs.F
     Src/qsort.F
     Src/schecomm.F
     Src/setup_H0.F
     Src/walltime.c
     Src/xc.f
     Tests/mgc-force/
     Tests/mgc-force/makefile
     Tests/mgc-force/mgc-force.fdf
     Tests/mgc-force/mgc-force.pseudos
     Util/Denchar/Src/timer_local.f
   renamed:
     Src/atom.f => Src/atom.F
     Src/forhar.f => Src/forhar.F
     Src/hsparse.f => Src/hsparse.F
     Src/meshmatrix.F => Src/meshdscf.F
     Src/parallelsubs.f => Src/parallelsubs.F
     Src/rhooda.f => Src/rhooda.F
     Src/rhoofd.f => Src/rhoofd.F
     Src/sparse_matrices.F90 => Src/sparse_matrices.F
     Src/timer.f90 => Src/timer.F90
     Src/vmat.f => Src/vmat.F
   modified:
     Docs/siesta.ind
     Docs/siesta.tex
     Src/MPI/generate.sh
     Src/MPI/mpi.F
     Src/MPI/mpi__include.f90
     Src/Makefile
     Src/SiestaXC/lib-makefile
     Src/SiestaXC/makefile
     Src/Sys/finisterrae.make
     Src/Sys/gfortran-netcdf.make
     Src/Sys/marenostrum-mpi-32.make
     Src/Sys/nano-intel-mpi-cdf.make
     Src/Sys/nano-intel-mpi.make
     Src/alloc.F90
     Src/arw.f
     Src/atm_transfer.f
     Src/atm_types.f
     Src/atmfuncs.f
     Src/atomlist.f
     Src/atomlwf.F
     Src/bands.F
     Src/basis_io.F
     Src/basis_specs.f
     Src/basis_types.f
     Src/bonds.f
     Src/born_charge.F
     Src/broadcast_basis.F
     Src/broyden_optim.F
     Src/cdiag.F
     Src/cell_broyden_optim.F
     Src/cell_fire_optim.F
     Src/cgvc.F
     Src/cgvc_zmatrix.F
     Src/cgwf.F
     Src/chemical.f
     Src/chempot.F
     Src/compute_dm.F
     Src/conjgr.f
     Src/constr.f
     Src/coor.F
     Src/denmat.F
     Src/denmatlomem.F
     Src/densematrix.f
     Src/detover.F
     Src/dfscf.f
     Src/dhscf.F
     Src/diag2g.F
     Src/diag2k.F
     Src/diagg.F
     Src/diagk.F
     Src/diagk_file.F
     Src/diagkp.F
     Src/diagon.F
     Src/diagpol.f
     Src/diagsprl.F
     Src/dipole.F
     Src/dnaefs.f
     Src/dynamics.f
     Src/egandd.F
     Src/eggbox.F
     Src/electrostatic.f
     Src/ener3.F
     Src/ener3lomem.F
     Src/fdf/Otherfile
     Src/fdf/README
     Src/fdf/coords.fdf
     Src/fdf/fdf.Standard
     Src/fdf/makefile
     Src/fdf/sample.fdf
     Src/fermid.F
     Src/fft.F
     Src/find_kgrid.F
     Src/fire_optim.F
     Src/fixed.F
     Src/get_target_stress.f
     Src/globalise.F
     Src/gradient.F
     Src/gradientlomem.F
     Src/grdsam.F
     Src/initatom.f
     Src/initparallel.F
     Src/iocg.f
     Src/iodm.F
     Src/iodm_netcdf.F90
     Src/iodmhs_netcdf.F90
     Src/ioeig.f
     Src/iofa.f
     Src/iogrid_netcdf.F90
     Src/iokp.f
     Src/iolwf.F
     Src/iomd.f
     Src/iopipes.F90
     Src/ioxv.F
     Src/iozm.F
     Src/kgrid.F
     Src/kgridinit.F
     Src/kinefsm.f
     Src/kpoint_grid.F90
     Src/kpoint_pdos.F90
     Src/ksv.f
     Src/ksvinit.F
     Src/linpack.F
     Src/listsc.f
     Src/local_DOS.F
     Src/m_broyden_mixing.f
     Src/m_fire_mixing.f
     Src/m_history.f90
     Src/m_iodm.F
     Src/m_iorho.F
     Src/m_iostruct.f
     Src/m_memory.F
     Src/m_mpi_utils.F
     Src/m_pulay.F90
     Src/m_sparse.F
     Src/m_spin.F90
     Src/m_ts_in_siesta.F
     Src/m_ts_kpoints.F90
     Src/m_ts_options.F90
     Src/matel.f
     Src/memory.F
     Src/memoryinfo.F
     Src/meshsubs.F
     Src/metaforce.F
     Src/mixer.F
     Src/mneighb.f
     Src/molecularmechanics.F90
     Src/moreParallelSubs.F90
     Src/mulliken.F
     Src/naefs.f
     Src/new_dm.F
     Src/nlefsm.f
     Src/normalize_dm.F
     Src/obc.f
     Src/old_atmfuncs.f
     Src/optical.F
     Src/ordern.F
     Src/outcoor.f
     Src/overfsm.f
     Src/overlap.f
     Src/pdos.F
     Src/pdosg.F
     Src/pdosk.F
     Src/pdoskp.F
     Src/phirphi.f
     Src/phirphi_opt.f
     Src/phonon.F
     Src/pixmol.f
     Src/plcharge.F
     Src/poison.F
     Src/post_scf_work.F
     Src/precision.F
     Src/projected_DOS.F
     Src/proximity_check.F
     Src/pseudopotential.f
     Src/pxf.F90
     Src/radfft.f
     Src/radial.f
     Src/rdiag.F
     Src/read_xc_info.F
     Src/readsp.F
     Src/redcel.F
     Src/reinit.F
     Src/reoptical.F
     Src/reord.f
     Src/rhoofdsp.f
     Src/savepsi.F
     Src/scfconvergence_test.F
     Src/setatomnodes.F
     Src/setspatial.f
     Src/setup_hamiltonian.F
     Src/setup_kscell.F
     Src/shaper.f
     Src/show_distribution.f
     Src/siesta.F
     Src/siesta_analysis.F
     Src/siesta_end.F
     Src/siesta_forces.F
     Src/siesta_init.F
     Src/siesta_move.F
     Src/siesta_options.F90
     Src/sorting.f
     Src/spher_harm.f
     Src/state_analysis.F
     Src/state_init.F
     Src/struct_init.F
     Src/vmatsp.f
     Src/vmb.F
     Src/write_md_record.F
     Src/write_subs.F
     Src/writewave.F
     Src/wxml/m_wxml_array_str.f90
     Src/wxml/m_wxml_buffer.f90
     Src/wxml/m_wxml_core.f90
     Src/wxml/m_wxml_dictionary.f90
     Src/wxml/m_wxml_elstack.f90
     Src/zm_broyden_optim.F
     Src/zm_fire_optim.F
     Src/zmatrix.F
     Tests/Makefile
     Util/Contrib/APostnikov/Makefile
     Util/Denchar/Src/Makefile
     Util/Denchar/Src/atompla.f
     Util/Denchar/Src/denchar.f
     Util/Denchar/Src/local_reinit.f
     Util/Denchar/Src/readpla.f
     Util/Denchar/Src/readsts.f
     Util/Gen-basis/Makefile
     Util/Gen-basis/gen-basis.F
     Util/HSX/iohs.F
     Util/Helpers/Makefile
     Util/Optimizer/Makefile
     Util/TBTrans/Makefile
     Util/TBTrans/green4.F
     Util/TBTrans/m_tbt_gf.F90
     Util/TBTrans/m_tbt_options.F90
     Util/TBTrans/mkqgrid.f
     Util/TBTrans/reinit_tb.F
     Util/TBTrans/tbtrans.F
     Util/VCA/Makefile
     Util/Vibra/Src/Makefile
     Util/Vibra/Src/fcbuild.f
     Util/Vibra/Src/klines.f
     Util/Vibra/Src/recoor.f
     Util/Vibra/Src/vibrator.f
     version.info
     Src/atom.F
     Src/forhar.F
     Src/hsparse.F
     Src/meshdscf.F
     Src/parallelsubs.F
     Src/rhooda.F
     Src/rhoofd.F
     Src/sparse_matrices.F
     Src/timer.F90
     Src/vmat.F

