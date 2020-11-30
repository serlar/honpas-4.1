#
#  Makefile for the Libint1-Glue library.
#  
.PHONY: check_flags libint_glue
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
# Select version of Libint1-Glue
#
LIBINT1_GLUE_VERSION = 0.1.0
#
# Include copying operations in rule to make sure that
# they are always performed.
#
# Library module libgaufre.a
default: module
module: check_flags libint_glue
	cp include/*.h ..
	cp lib/*.a ..
#
libint_glue:
	rm -rf libint1-glue-$(LIBINT1_GLUE_VERSION)
	tar xzf $(LIBINT_GLUE_TARDIR)/libint1-glue-$(LIBINT1_GLUE_VERSION).tar.gz
	mkdir libint1-glue-$(LIBINT1_GLUE_VERSION)/_siesta_build
	cd libint1-glue-$(LIBINT1_GLUE_VERSION)/_siesta_build && \
	  ../configure \
	    --prefix="$(PWD)" \
	    --enable-static --disable-shared \
	    --with-libint1="$(PWD)/../Libint" \
	    CC="$(CC)" CFLAGS="$(CFLAGS)" \
	    CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" && \
	  make && \
	  make install
#
check_flags:
	@echo "In libint_glue, CC is       : $(CC)"
	@echo "In libint_glue, CFLAGS is   : $(CFLAGS)"
	@echo "In libint_glue, CXX is      : $(CXX)"
	@echo "In libint_glue, CXXFLAGS is : $(CXXFLAGS)"
#
# Cleaning
#
clean:
	rm -rf libint1-glue-$(LIBINT1_GLUE_VERSION)
	rm -rf bin include lib share
