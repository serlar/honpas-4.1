#
#  Makefile for the GAUFRE library.
#  
.PHONY: check_flags
#
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location 
# of the arch.make file.
#
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
# Add items to the INCFLAGS as set by the arch.make file
#
INCFLAGS:=-I$(VPATH) $(INCFLAGS)    # For VPATH operation
#
# This is needed to pick up top-level module files, such as mpi_siesta
#
LOCAL_INCFLAGS=-I../
INCFLAGS:=$(INCFLAGS) $(LOCAL_INCFLAGS)
#
# List of library objects
#
minpack_objs = \
  chkder.o \
  dogleg.o \
  enorm.o \
  enorm2.o \
  fdjac1.o \
  fdjac2.o \
  hybrd.o \
  hybrd1.o \
  hybrj.o \
  hybrj1.o \
  lmder.o \
  lmder1.o \
  lmdif.o \
  lmdif1.o \
  lmpar.o \
  lmstr.o \
  lmstr1.o \
  qform.o \
  qrfac.o \
  qrsolv.o \
  r1mpyq.o \
  r1updt.o \
  r8_uniform_01.o \
  r8mat_print.o \
  r8mat_print_some.o \
  r8vec_print.o \
  rwupdt.o \
  timestamp.o
gaufre_objs = \
  gaufre.o \
  gaufre_common.o \
  gaufre_config.o \
  gaufre_data.o \
  gaufre_driver.o \
  gaufre_io.o \
  gaufre_minpack.o \
  gaufre_orbital.o
#
# List of library modules
#
gaufre_mods = \
  gaufre.mod \
  gaufre_common.mod \
  gaufre_config.mod \
  gaufre_data.mod \
  gaufre_driver.mod \
  gaufre_io.mod \
  gaufre_minpack.mod \
  gaufre_orbital.mod
#
# Include copying operations in rule to make sure that
# they are always performed.
#
# Library module libgaufre.a
default: module
module: libgaufre.a
	cp libgaufre.a $(gaufre_mods) ..
#
libgaufre.a: check_flags $(minpack_objs) $(gaufre_objs)
	$(AR) $(ARFLAGS_EXTRA) cru libgaufre.a $(minpack_objs) $(gaufre_objs)
	-$(RANLIB) libgaufre.a

check_flags:
	@echo "In gaufre, INCFLAGS is: $(INCFLAGS)"
#
# Cleaning
#
clean:
	rm -f *.o *.mod libgaufre.a
#
# Explicit dependencies for MINPACK
#
dogleg.o: \
  enorm.o
hybrd.o: \
  dogleg.o \
  enorm.o \
  fdjac1.o \
  qform.o \
  qrfac.o \
  r1mpyq.o \
  r1updt.o
hybrd1.o: \
  hybrd.o
hybrj.o: \
  dogleg.o \
  enorm.o \
  qrfac.o \
  qform.o \
  r1mpyq.o \
  r1updt.o
hybrj1.o: \
  hybrj.o
lmder.o: \
  enorm.o \
  lmpar.o \
  qrfac.o
lmder1.o: \
  lmder.o
lmdif.o: \
  enorm.o \
  fdjac2.o \
  lmpar.o \
  qrfac.o
lmdif1.o: \
  lmdif.o
lmpar.o: \
  enorm.o \
  qrsolv.o
lmstr.o: \
  enorm.o \
  lmpar.o \
  qrfac.o \
  rwupdt.o
lmstr1.o: \
  lmstr.o
qrfac.o: \
  enorm.o
r8mat_print.o: \
  r8mat_print_some.o
#
# Explicit dependencies for GAUFRE
#
gaufre.o: \
  gaufre_common.o \
  gaufre_config.o \
  gaufre_data.o \
  gaufre_driver.o \
  gaufre_io.o
gaufre_config.o: \
  gaufre_common.o
gaufre_data.o: \
  gaufre_common.o \
  gaufre_config.o
gaufre_driver.o: \
  gaufre_common.o \
  gaufre_config.o \
  gaufre_minpack.o
gaufre_io.o: \
  gaufre_common.o \
  gaufre_orbital.o
gaufre_minpack.o: \
  gaufre_common.o \
  gaufre_data.o
gaufre_orbital.o: \
  gaufre_common.o
