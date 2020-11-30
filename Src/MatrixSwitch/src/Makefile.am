# -*- Automake -*-
#
# Makefile for the MatrixSwitch package
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the MatrixSwitch software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

                    # ------------------------------------ #

#
# Main source files
#

# Common source files
msw_srcs = \
  MatrixSwitch.F90 \
  MatrixSwitch_m_copy.F90 \
  MatrixSwitch_m_register.F90 \
  MatrixSwitch_ops.F90 \
  MatrixSwitch_wrapper_params.F90 \
  MatrixSwitch_m_add.F90 \
  MatrixSwitch_mm_multiply.F90 \
  MatrixSwitch_m_set.F90 \
  MatrixSwitch_wrapper.F90

# Fortran modules
# TODO: support both lower-case and upper-case file names
msw_f03_mods = \
  matrixswitch.$(MODEXT) \
  matrixswitch_m_copy.$(MODEXT) \
  matrixswitch_m_register.$(MODEXT) \
  matrixswitch_ops.$(MODEXT) \
  matrixswitch_wrapper_params.$(MODEXT) \
  matrixswitch_m_add.$(MODEXT) \
  matrixswitch_mm_multiply.$(MODEXT) \
  matrixswitch_m_set.$(MODEXT) \
  matrixswitch_wrapper.$(MODEXT)

# Libraries to install
lib_LTLIBRARIES = libMatrixSwitch.la

libMatrixSwitch_la_SOURCES = $(msw_srcs)
libMatrixSwitch_la_LDFLAGS = -version-info 0:0:0

# Fortran modules to install
f03moddir = $(includedir)
install-data-local:
	$(INSTALL) -d -m 755 $(DESTDIR)$(f03moddir)
	$(INSTALL) -m 644 $(msw_f03_mods) $(DESTDIR)$(f03moddir)

uninstall-local:
	cd $(DESTDIR)$(f03moddir) && rm -f $(msw_f03_mods)

# Local cleaning
clean-local:
	rm -f *.mod *.MOD *.obj

# Explicit dependencies
include $(srcdir)/dependencies.mf

# Files to mark explicitly for distribution
EXTRA_DIST = dependencies.mf
