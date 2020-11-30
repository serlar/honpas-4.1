# -*- Autoconf -*-
#
# M4 macros for MatrixSwitch
#
# Copyright (C) 2015 Yann Pouillon
#
# This file is part of the MatrixSwitch software package. For license information,
# please see the COPYING file in the top-level directory of the MatrixSwitch source
# distribution.
#

#
# MPI support
#



# MSW_MPI_DETECT()
# ----------------
#
# Checks whether the configured MPI implementation is working.
#
AC_DEFUN([MSW_MPI_DETECT], [
  dnl Init
  msw_mpi_ok="unknown"

  dnl Display current MPI status
  AC_MSG_CHECKING([how MPI parameters have been set])
  AC_MSG_RESULT([${msw_mpi_type}])
  AC_MSG_CHECKING([whether the MPI C compiler is set])
  AC_MSG_RESULT([${msw_mpi_cc_set}])
  AC_MSG_CHECKING([whether the MPI C compiler is wrapped])
  AC_MSG_RESULT([${msw_mpi_cc_wrap}])
  AC_MSG_CHECKING([whether the MPI Fortran compiler is set])
  AC_MSG_RESULT([${msw_mpi_fc_set}])
  AC_MSG_CHECKING([whether the MPI Fortran compiler is wrapped])
  AC_MSG_RESULT([${msw_mpi_fc_wrap}])

  dnl Warn if serial component of wrapped compilers supports MPI
  if test "${msw_mpi_cc_wrap}" = "yes"; then
    AC_MSG_NOTICE([validating that '${msw_sercc}' is indeed serial])
    AC_MSG_NOTICE([please ignore possible warnings about mpi.h not found])
    _MSW_MPI_CHECK_CC([${msw_sercc}])
    if test "${msw_mpi_cc_ok}" = "yes"; then
      AC_MSG_WARN([the serial C compiler is MPI-aware
                    Your current configuration is probably ill-defined.
                    The build will likely fail.])
      sleep 5
    fi
  fi
  if test "${msw_mpi_fc_wrap}" = "yes"; then
    AC_MSG_NOTICE([validating that '${msw_serfc}' is indeed serial])
    _MSW_MPI_CHECK_FC([${msw_serfc}])
    if test "${msw_mpi_fc_ok}" = "yes"; then
      AC_MSG_WARN([the serial Fortran compiler is MPI-aware
                    Your current configuration is probably ill-defined.
                    The build will likely fail.])
      sleep 5
    fi
  fi

  dnl Test MPI compilers
  _MSW_MPI_CHECK_CC([${CC}])
  if test "${msw_mpi_cc_ok}" = "yes"; then
    _MSW_MPI_CHECK_FC([${FC}])
  fi

  dnl Look for mpirun
  dnl FIXME: hard-coded command-line options
  if test "${MPIRUN}" = ""; then
    AC_CHECK_PROGS([MPIRUN], [mpirun mpiexec])
  fi
  if test "${MPIRUN}" != ""; then
    MPIRUN="${MPIRUN} -np 4"
  fi

  dnl Take final decision
  AC_MSG_CHECKING([whether we have a full MPI support])
  if test "${msw_mpi_cc_ok}" = "yes" -a \
          "${msw_mpi_fc_ok}" = "yes"; then
    msw_mpi_ok="yes"
  else
    msw_mpi_ok="no"
  fi
  AC_MSG_RESULT([${msw_mpi_ok}])
]) # MSW_MPI_DETECT



# MSW_MPI_INIT()
# --------------
#
# Initializes MPI parameters.
#
AC_DEFUN([MSW_MPI_INIT], [
  if test "${msw_mpi_enable}" != "no"; then
    AC_MSG_CHECKING([how MPI parameters have been set])
    AC_MSG_RESULT([${msw_mpi_type}])
    if test "${msw_mpi_type}" = "env"; then
      _AC_SRCDIRS(["."])
    fi
    _MSW_MPI_INIT_CC
    _MSW_MPI_INIT_FC
  fi
]) # MSW_MPI_INIT



                    ########################################



# _MSW_MPI_CHECK_CC(CC)
# ---------------------
#
# Check whether the MPI C compiler is working.
#
AC_DEFUN([_MSW_MPI_CHECK_CC], [
  dnl Init
  msw_mpi_cc_ok="unknown"
  msw_mpi_cc_has_funs="unknown"
  msw_mpi_cc_has_hdrs="unknown"

  dnl Prepare environment
  msw_saved_CC="${CC}"
  msw_saved_CC="${CC}"
  CC="$1"
  tmp_mpi_header=mpi.h
  tmp_mpi_cache=AS_TR_SH([ac_cv_header_${tmp_mpi_header}])
  ${as_unset} ${tmp_mpi_cache}

  dnl Look for C includes
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([mpi.h],
    [msw_mpi_cc_has_hdrs="yes"], [msw_mpi_cc_has_hdrs="no"])
  AC_LANG_POP([C])

  dnl Look for C functions
  if test "${msw_mpi_cc_has_hdrs}" = "yes"; then
    AC_CHECK_FUNC([MPI_Init], [msw_mpi_cc_has_funs="yes"],
      [msw_mpi_cc_has_funs="no"])
  fi

  dnl Validate C support
  AC_MSG_CHECKING([whether the MPI C compiler works])
  if test "${msw_mpi_cc_has_funs}" = "yes" -a \
          "${msw_mpi_cc_has_hdrs}" = "yes"; then
    msw_mpi_cc_ok="yes"
  else
    msw_mpi_cc_ok="no"
  fi
  AC_MSG_RESULT([${msw_mpi_cc_ok}])

  dnl Restore environment
  CC="${msw_saved_CC}"
  unset tmp_mpi_cache
  unset tmp_mpi_header
]) # _MSW_MPI_CHECK_CC



# _MSW_MPI_CHECK_FC(FC)
# ---------------------
#
# Check whether the MPI Fortran compiler is working.
#
AC_DEFUN([_MSW_MPI_CHECK_FC], [
  dnl Init
  msw_mpi_fc_ok="unknown"
  msw_mpi_fc_has_funs="unknown"
  msw_mpi_fc_has_mods="unknown"

  dnl Prepare environment
  msw_saved_FC="${FC}"
  FC="$1"

  dnl Look for Fortran modules
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([for a Fortran MPI module])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use mpi
    ]])], [msw_mpi_fc_has_mods="yes"], [msw_mpi_fc_has_mods="no"])
  AC_MSG_RESULT([${msw_mpi_fc_has_mods}])
  AC_LANG_POP([Fortran])

  dnl Look for Fortran functions
  if test "${msw_mpi_fc_has_mods}" = "yes"; then
    AC_LANG_PUSH([Fortran])
    AC_MSG_CHECKING([for a Fortran MPI_Init])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use mpi
        integer :: ierr
        call mpi_init(ierr)
      ]])], [msw_mpi_fc_has_funs="yes"], [msw_mpi_fc_has_funs="no"])
    AC_MSG_RESULT([${msw_mpi_fc_has_funs}])
    AC_LANG_POP([Fortran])
  fi

  dnl Validate Fortran support
  AC_MSG_CHECKING([whether the MPI Fortran compiler works])
  if test "${msw_mpi_fc_has_funs}" = "yes" -a \
          "${msw_mpi_fc_has_mods}" = "yes"; then
    msw_mpi_fc_ok="yes"
  else
    msw_mpi_fc_ok="no"
  fi
  AC_MSG_RESULT([${msw_mpi_fc_ok}])

  dnl Restore environment
  FC="${msw_saved_FC}"
]) # _MSW_MPI_CHECK_FC



# _MSW_MPI_INIT_CC()
# ------------------
#
# Initializes MPI parameters related to the C compiler.
#
AC_DEFUN([_MSW_MPI_INIT_CC], [
  dnl Init
  msw_sercc="${CC}"
  msw_mpicc=""
  msw_mpi_cc_set="no"
  msw_mpi_cc_wrap="unknown"

  dnl Look for a MPI C compiler
  case "${msw_mpi_type}" in

    def)
      msw_mpi_cc_wrap="no"
      ;;

    dir)
      msw_mpicc="${with_mpi}/bin/mpicc"
      if test -x "${msw_mpicc}"; then
        AC_MSG_CHECKING([for an executable MPI C compiler])
        AC_MSG_RESULT([${msw_mpicc}])
        if test "${msw_sercc}" = ""; then
          AC_MSG_NOTICE([setting CC to '${msw_mpicc}'])
          CC="${msw_mpicc}"
          msw_mpi_cc_set="yes"
          msw_mpi_cc_wrap="no"
        else
          msw_mpi_cc_wrap="yes"
        fi
      else
        AC_MSG_ERROR([MPI C compiler not found in ${with_mpi}/bin])
      fi
      ;;

    env|yon)
      if test -n "${MPICC}"; then
        msw_mpicc="${MPICC}"
      else
        AC_CHECK_PROGS([msw_mpicc], [mpicc])
      fi
      if test -n "${msw_sercc}" -a -n "${msw_mpicc}"; then
        msw_mpi_cc_wrap="yes"
      elif test -n "${msw_mpicc}"; then
        AC_MSG_NOTICE([setting CC to '${msw_mpicc}'])
        CC="${msw_mpicc}"
        msw_mpi_cc_set="yes"
        msw_mpi_cc_wrap="no"
      fi
      ;;

  esac

  if test "${msw_mpi_cc_wrap}" = "yes"; then
    _MSW_MPI_CREATE_WRAPPER([CC], [${msw_sercc}], [${msw_mpicc}])
    msw_mpi_cc_set="yes"
  fi
]) # _MSW_MPI_INIT_CC



# _MSW_MPI_INIT_FC()
# ------------------
#
# Initializes MPI parameters related to the Fortran compiler.
#
AC_DEFUN([_MSW_MPI_INIT_FC], [
  dnl Init
  msw_serfc="${FC}"
  msw_mpifc=""
  msw_mpi_fc_set="no"
  msw_mpi_fc_wrap="unknown"

  dnl Look for a MPI Fortran compiler
  case "${msw_mpi_type}" in

    def)
      msw_mpi_fc_wrap="no"
      ;;

    dir)
      msw_mpifc="${with_mpi}/bin/mpif90"
      if test -x "${msw_mpifc}"; then
        AC_MSG_CHECKING([for an executable MPI Fortran compiler])
        AC_MSG_RESULT([${msw_mpifc}])
        if test "${msw_serfc}" = ""; then
          AC_MSG_NOTICE([setting FC to '${msw_mpifc}'])
          FC="${msw_mpifc}"
          msw_mpi_fc_set="yes"
          msw_mpi_fc_wrap="no"
        else
          msw_mpi_fc_wrap="yes"
        fi
      else
        AC_MSG_ERROR([MPI Fortran compiler not found in ${with_mpi}/bin])
      fi
      ;;

    env|yon)
      if test -n "${MPIFC}"; then
        msw_mpifc="${MPIFC}"
      else
        AC_CHECK_PROGS([msw_mpifc], [mpif90 mpif95])
      fi
      if test -n "${msw_serfc}" -a -n "${msw_mpifc}"; then
        msw_mpi_fc_wrap="yes"
      elif test -n "${msw_mpifc}"; then
        AC_MSG_NOTICE([setting FC to '${msw_mpifc}'])
        FC="${msw_mpifc}"
        msw_mpi_fc_set="yes"
        msw_mpi_fc_wrap="no"
      fi
      ;;

  esac

  if test "${msw_mpi_fc_wrap}" = "yes"; then
    _MSW_MPI_CREATE_WRAPPER([FC], [${msw_serfc}], [${msw_mpifc}])
    msw_mpi_fc_set="yes"
  fi
]) # _MSW_MPI_INIT_FC



# _MSW_MPI_CREATE_WRAPPER(COMPILER_TYPE, SERIAL_COMPILER, MPI_COMPILER)
# ---------------------------------------------------------------------
#
# Creates a wrapper for MPI compilers when they can be configured to
# accept different serial compilers.
#
# Note: it is impossible to set two compiler levels with the Autotools,
#       because Automake requires CC, CXX, and FC to be set to
#       the actual compilers.
#
AC_DEFUN([_MSW_MPI_CREATE_WRAPPER], [
  dnl Init
  tmp_comp_name=`echo "$1" | sed -e 's/.*/\L&/'`
  ${MKDIR_P} config/wrappers

  dnl Create file
  cat >config/wrappers/wrap-mpi${tmp_comp_name} <<EOF
#!/bin/sh

$1="$2"
export $1

$3 \[$]{*}
EOF

  dnl Fix permissions
  chmod u+x config/wrappers/wrap-mpi${tmp_comp_name}

  dnl Overwrite compiler setting
  eval tmp_wrapper_path="${ac_abs_top_builddir}/config/wrappers/wrap-mpi${tmp_comp_name}"
  tmp_wrapper_name=`basename "${tmp_wrapper_path}"`
  AC_MSG_NOTICE([wrapping serial and MPI compilers into ${tmp_wrapper_name}])
  $1="${tmp_wrapper_path}"

  dnl Clean-up
  unset tmp_comp_name
  unset tmp_wrapper_name
  unset tmp_wrapper_path
]) # _MSW_MPI_CREATE_WRAPPER
