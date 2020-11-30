# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#-------------------------------------------------------------------
# DOCUMENTED arch.make
#
# The most useful makefile symbols are explained. Use this file as
# a guide when you are looking at the .make files in this directory,
# or after 'configure' has produced a draft arch.make for you.
#
# This block tells make to consider only these suffixes in its operation
# It is included in most makefiles in the source tree, but it does not
# hurt to have it here too.
.SUFFIXES:
.SUFFIXES: .f .F .o .c .a .f90 .F90

# This string will be copied to the executable, you may write any
# descriptive statement regarding the compilation setup
SIESTA_ARCH = x86_64-honpas-linux-gnu

# CC is the C compiler.
# This currently need not be an MPI C-compiler for MPI compilation.
CC = mpicc

# CXX is the C++ compiler.
# This currently need not be an MPI C-compiler for MPI compilation.
CXX = mpic++

# In case your compiler does not understand the special meaning of 
# the .F and .F90 extensions ("files in need of preprocessing"), you
# will need to use an explicit preprocessing step.
# Typically this is sufficient to be the compiler with -E -P -x c
FPP = $(FC) -E -P -x c

# FC is typically the name of your fortran compiler. It is not always
# a good idea to add options here, except when they are essential for
# a proper operation.
# In case of MPI compilation this should be the mpi compiler (mpifort).
FC = mpifort

# The FC_SERIAL symbol is useful in at least two cases:
#   1. When the "MPI compiler environment" is so complex that it might
#      trick the configure scripts.
#   2. When executables compiled with a (parallel) FC are flagged by 
#      the computer centers as "queuing-system-only". 
# Most utilities are thus compiled with FC_SERIAL, which in practice
# defaults to FC if it is not defined.
FC_SERIAL = gfortran

# Here we should put mainly optimization flags for C
CFLAGS_OPTIM = -O2 -march=native
CFLAGS = -g $(CFLAGS_OPTIM)

# Here we should put mainly optimization flags for C++
CXXFLAGS_OPTIM = -O2 -march=native
CXXFLAGS = -g -std=c++11 $(CXXFLAGS_OPTIM)
CXX_LIBS =

# Here we should put mainly optimization flags for Fortran
FFLAGS_OPTIM = -O2 -march=native
FFLAGS = -g -fPIC $(FFLAGS_OPTIM)

# If you want to use a specific archive function you may specify it here
# Note that if you use LTO(gcc)/IPO(intel) specifying the correct AR
# is REQUIRED
AR = ar
# Extra flags for library creation by the 'ar' command
ARFLAGS_EXTRA =

# Some systems do not have 'ranlib'. If so, use "echo" instead of "ranlib"
RANLIB = :

# A compiler-specific file holding special versions of some routines
# For most f95 compilers, "nag" should work. (The name is historical)
# The currently supported systems are:
#   nag (generic)
#   bsd
#   ibm
#   ibm_pessl
#   sgi
#   t3e
#   xlf
SYS = nag

# These symbols should not need to be specified. They will be detected
# automatically at the time of compiling the MPI interface. Set them
# only if the automatic detection fails and you are sure of their values.
#SP_KIND = 4
#DP_KIND = 8
#KINDS = $(SP_KIND) $(DP_KIND)

# Some compilers (notably IBM's) are not happy with the standard 
# syntax for definition of preprocessor symbols ( -DSOME_SYMBOL),
# and thy need a prefix (i.e. -WF,-DSOME_SYMBOL). This is used
# in some utility makefiles.
# Typically this need not be defined
DEFS_PREFIX =

# Used only at the linking stage. For example, you might need "-static"
LDFLAGS =

# These symbols help to keep the building rules concise
# This enables specific compilation options for certain
# source files.
FCFLAGS_fixed_f =
FCFLAGS_free_f90 =
FPPFLAGS_fixed_F =
FPPFLAGS_free_F90 =

# This is the most installation-dependent part
# We can make things a bit easier by grouping symbols, and maybe
# using the -L flag to define search directories (see examples
# in this directory).
# For the most simplistic compilation one requires the following
# libraries:
#   BLAS
#   LAPACK
#   ScaLAPACK (only for MPI compilation)
BLAS_LIBS =
LAPACK_LIBS = -lopenblas
# The most recent ScaLAPACK versions have built in the BLACS library.
# If BLACS is not included in the ScaLAPACK library you are
# required to add it here as well:
BLACS_LIBS =
SCALAPACK_LIBS = -lscalapack

# If you do not have BLAS/LAPACK installed on your machine
# you may use the shipped BLAS/LAPACK versions.
#
# If you do not have the LAPACK library you may do:
# COMP_LIBS = libsiestaLAPACK.a
# If you have neither the BLAS, nor LAPACK library you may do:
# COMP_LIBS = libsiestaLAPACK.a libsiestaBLAS.a
#
# You are HIGHLY encouraged to use an optimized BLAS/LAPACK
# library as SIESTA performance is mainly governed by these
# libraries.
COMP_LIBS = 

# For netCDF support. Make sure you get a version compatible
# with the other options (for example, 32/64 bit). Don't forget
# to set -DCDF below.
# To use NetCDF define the installation directory of the library.
# By default NetCDF is not used, but it may become mandatory in
# later versions of SIESTA.
#NETCDF_ROOT = $(EBROOTNETCDFMINFORTRAN)
#NETCDF_INCFLAGS = -I$(EBROOTNETCDF)/include -I$(EBROOTNETCDFMINFORTRAN)/include
#NETCDF_LIBS = \
#  $(EBROOTNETCDFMINFORTRAN)/lib/libnetcdff.a \
#  $(EBROOTNETCDF)/lib/libnetcdf.a \
#  $(EBROOTHDF5)/lib/libhdf5hl_fortran.a \
#  $(EBROOTHDF5)/lib/libhdf5_fortran.a \
#  $(EBROOTHDF5)/lib/libhdf5_hl.a \
#  $(EBROOTHDF5)/lib/libhdf5.a \
#  $(EBROOTCURL)/lib/libcurl.a \
#  $(EBROOTSZIP)/lib/libsz.a \
#  -lssl -lcrypto -lz -ldl
#DEFS_NETCDF = $(DEFS_PREFIX)-DCDF

# To fit natural orbitals on Gaussians for hybrid XC functionals.
#GAUFRE_ROOT = $(HONPAS_HOME)/workdir/gaufre/gnu/check/$(GAUFRE_VERSION)
#GAUFRE_INCFLAGS = -I$(GAUFRE_ROOT)/include
#GAUFRE_LIBS = $(GAUFRE_ROOT)/lib/libgaufre.a
COMP_LIBS += libgaufre.a

# For hybrid XC functionals. Make sure you get a version compatible
# with the other options (for example, 32/64 bit).
#LIBINT_ROOT = $(HONPAS_HOME)/workdir/libint/gnu/check/$(LIBINT_VERSION)
#LIBINT_INCFLAGS = -I$(LIBINT_ROOT)/include
#LIBINT_LIBS = \
#  $(LIBINT_ROOT)/lib/libint1_glue.a \
#  $(LIBINT_ROOT)/lib/libderiv.a \
#  $(LIBINT_ROOT)/lib/libint.a
DEFS_LIBINT = $(DEFS_PREFIX)-DHAVE_LIBINT
COMP_LIBS += libint1_glue.a
COMP_LIBS += libint.a

# This (as well as the -DMPI definition) is essential for MPI support
MPI_INTERFACE = libmpi_f90.a
MPI_INCLUDE = .
DEFS_MPI = $(DEFS_PREFIX)-DMPI

# Preprocessor definitions or flags.
# Here we use FPPFLAGS (as 'configure' calls them), but historically
# it was very common to use DEFS. Try to use only FPPFLAGS from now on,
# converting any old arch.make files you might have lying around, and
# remember that you have to change the final building rules at the end
# to use only FPPFLAGS. DEFS is deprecated.
DEFS_GRID = $(DEFS_PREFIX)-DGRID_SP
DEFS_FORTRAN = $(DEFS_PREFIX)-DFC_HAVE_ABORT

# CDF and MPI are self-explanatory
# Other definitions might be needed to work around some glitch in the compiler
# See the manual for more details.
FPPFLAGS = $(DEFS_FORTRAN) $(DEFS_GRID) $(DEFS_LIBINT) $(DEFS_MPI) $(DEFS_NETCDF)

# We put here all the neeeded libraries.
# Sometimes the BLAS are included in LAPACK (or it could be that everything
# is included in SCALAPACK...). You might need to experiment if you find 
# duplicate symbols. See examples in this directory.
#
LIBS = $(NETCDF_LIBS) $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(COMP_LIBS)

# Dependency rules ---------

# Some compilers are not able to compile certain files with full optimization,
# or they produce wrong results if they do. For example, the PGI compiler 
# has trouble with atom.f and electrostatic.f. In these cases, we need to
# insert extra lines. Use exactly the format shown, as it is general enough
# to work with VPATH.
#
FFLAGS_DEBUG = -g3 -ggdb -fPIC -fbacktrace -Og

# The atom.f code is very vulnerable. Particularly the Intel compiler
# will make an erroneous compilation of atom.f with high optimization
# levels.
atom.o: atom.F
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) $<

# Finally, the default building rules which will be used everywhere,
# unless overriden.
# These were created by a former run of 'configure'.
# See other examples in this directory. If you cut and paste, 
# MAKE SURE that there are TABS, not spaces, at the beginning.
#
# Important points to note:
#  - INCFLAGS must be present. It is used in several utility makefiles
#  - Either FPPFLAGS (preferred) or DEFS (deprecated) must be present
#    (see above) -- Note that the use of DEFS might break Util compilations.
#  - If your compiler does not recognize .F and .F90 extensions as in
#    need of preprocessing, you will need to use an intermediate
#    preprocessing step (see above about FPP). For example:
##
#.F90.o:
#        $(FPP) $(FPPFLAGS) $< > tmp_$*.f90
#        $(FC) -c $(FFLAGS) $(INCFLAGS) tmp_$*.f90
#        @mv tmp_$*.o $*.o
#        @rm -f tmp_$*.f90
#
#---------------- Example of actual rules  
.c.o:
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(CPPFLAGS) $< 
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< 
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<
