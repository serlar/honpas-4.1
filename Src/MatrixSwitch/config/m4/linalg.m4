# -*- Autoconf -*-
#
# M4 macros for MatrixSwitch
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the MatrixSwitch software package. For license information,
# please see the COPYING file in the top-level directory of the MatrixSwitch source
# distribution.
#

#
# Linear algebra support
#



# MSW_LINALG_DETECT()
# -------------------
#
# Checks that the selected linear algebra libraries properly work.
#
AC_DEFUN([MSW_LINALG_DETECT],[
  dnl Init
  msw_linalg_has_lapack="unknown"
  msw_linalg_has_scalapack="unknown"
  msw_linalg_ok="unknown"

  dnl Prepare environment
  saved_CPPFLAGS="${CPPFLAGS}"
  saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${msw_linalg_incs}"
  LIBS="${msw_linalg_libs} ${LIBS}"

  dnl Check BLAS and LAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether linear algebra libraries work])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
      call zhpev
    ]])], [msw_linalg_ok="yes"; msw_linalg_has_lapack="yes"], [msw_linalg_ok="no"])
  AC_MSG_RESULT([${msw_linalg_ok}])
  AC_LANG_POP([Fortran])

  dnl Check ScaLAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether linear algebra libraries have ScaLAPACK])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call pzheevx
    ]])], [msw_linalg_has_scalapack="yes"], [msw_linalg_has_scalapack="no"])
  AC_MSG_RESULT([${msw_linalg_ok}])
  AC_LANG_POP([Fortran])

  dnl Restore environment
  CPPFLAGS="${saved_CPPFLAGS}"
  LIBS="${saved_LIBS}"
]) # MSW_LINALG_DETECT
