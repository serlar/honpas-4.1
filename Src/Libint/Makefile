#
#  Makefile for the Libint library.
#  
.PHONY: check_flags libint_core
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
# Select versions of Libint and Libint1-Glue
#
LIBINT_VERSION = 1.1.4
#
# Include copying operations in rule to make sure that
# they are always performed.
#
# Library module libgaufre.a
default: module
module: check_flags libint_core
	cp -a include/* ..
	cp lib/*.a ..
#
libint_core:
	rm -rf libint-$(LIBINT_VERSION)
	tar xzf $(LIBINT_TARDIR)/libint-$(LIBINT_VERSION).tar.gz
	cd libint-$(LIBINT_VERSION) && \
	  ./configure \
	    --prefix="$(PWD)" \
	    --enable-static --disable-shared \
	    --with-cc-optflags="$(CFLAGS_OPTIM)" \
	    --with-cxx-optflags="$(CXXFLAGS_OPTIM)" \
	    CC="$(CC)" CFLAGS="$(CFLAGS)" \
	    CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" \
	    F77="$(FC)" FFLAGS="$(FCFLAGS)" \
	    F90="$(FC)" F90FLAGS="$(FCFLAGS)" \
	    FC="$(FC)" FCFLAGS="$(FCFLAGS)" && \
	  make -j8 SHELL="/bin/bash" && \
	  make install SHELL="/bin/bash"
#
check_flags:
	@echo "In libint, CC is       : $(CC)"
	@echo "In libint, CFLAGS is   : $(CFLAGS)"
	@echo "In libint, CXX is      : $(CXX)"
	@echo "In libint, CXXFLAGS is : $(CXXFLAGS)"
#
# Cleaning
#
clean:
	rm -rf libint-$(LIBINT_VERSION)
	rm -rf bin include lib share
