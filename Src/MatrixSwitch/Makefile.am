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

# Autotools parameters
ACLOCAL_AMFLAGS = -I config/m4

# Subdirectories containing makefiles
SUBDIRS = src examples doc

# Support for pkg-config (see http://pkg-config.freedesktop.org/wiki/)
pkgconfigdir = $(prefix)/lib/pkgconfig
pkgconfig_DATA = config/data/MatrixSwitch.pc

# Files to explicitly distribute in the source package
EXTRA_DIST = \
  config/data/MatrixSwitch.pc.in \
  config/data/MatrixSwitch-config.yaml.in \
  config/scripts/markdown

# Custom targets
distclean-local:
	rm -f $(distdir).tar.bz2 $(distdir).tar.gz $(distdir).tar.xz
