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
# pspBLAS support
#



# MSW_PSP_DETECT()
# -------------------
#
# Checks that the selected linear algebra libraries properly work.
#
AC_DEFUN([MSW_PSP_DETECT],[
  dnl Init
  msw_psp_ok="unknown"

  dnl Prepare environment
  saved_CPPFLAGS="${CPPFLAGS}"
  saved_FCFLAGS="${FCFLAGS}"
  saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${msw_psp_incs}"
  FCFLAGS="${FCFLAGS} ${msw_psp_incs}"
  LIBS="${msw_psp_libs} ${LIBS}"

  dnl Check pspBLAS routine
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether the pspBLAS library works])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use pspblas
      integer :: M,N,K
      type(psp_matrix_spm) :: A
      character(1) :: opA, opB
      double precision :: B(2,2),C(2,2),alpha,beta
      call psp_gespmm(M,N,K,A,opA,B,opB,C,alpha,beta)
    ]])], [msw_psp_ok="yes"], [msw_psp_ok="no"])
  AC_MSG_RESULT([${msw_psp_ok}])
  AC_LANG_POP([Fortran])

  dnl Restore environment
  CPPFLAGS="${saved_CPPFLAGS}"
  FCFLAGS="${saved_FCFLAGS}"
  LIBS="${saved_LIBS}"
]) # MSW_PSP_DETECT
