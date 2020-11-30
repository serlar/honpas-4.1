# Copyright (C) 2016 Yann Pouillon <notifications@materialsevolution.es>
#
# This file is part of the MatrixSwitch software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Configure script for the MatrixSwitch package
#

                    # ------------------------------------ #

#
# Autotools startup
#

# Init Autoconf
AC_PREREQ(2.69)
AC_INIT([MatrixSwitch], [0.0.0], [http://www.e-cam2020.eu:10080/ESL/omm/issues])
AC_REVISION([Autotools support for MatrixSwitch])
AC_CONFIG_AUX_DIR(config/gnu)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_SRCDIR([src/MatrixSwitch.F90])
_AC_SRCDIRS([.])

# Init Automake
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([1.14 parallel-tests color-tests])
AM_CONFIG_HEADER([config.h])

# Generate timestamp
msw_timestamp=`date '+%Y%m%dT%H%M%S%z'`
AC_SUBST(msw_timestamp)

                    # ------------------------------------ #

#
# System utilities
#

# Check for common programs
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_SED
AC_PROG_AWK
AC_PROG_GREP

                    # ------------------------------------ #

#
# Default settings
#

# Note for developers: you may edit these parameters to tune the behaviour
# of the build system. Please do it with care.

# Lists of options to manage (msw_*_enable_def variables)
#
# Note:
#   * switches (--enable-*) must be set to "yes" or "no"
#   * external dependencies (--with-*) must be set to "yes", "no", or "auto"
#
msw_cfg_switches="conv debug"
msw_cfg_extdeps="linalg mpi psp"

# Optional features
msw_conv_enable_def="no"
msw_debug_enable_def="yes"

# Linear algebra
# Note: default libraries must be empty because some compilers have
#       built-in support.
msw_linalg_enable_def="no"
msw_linalg_incs_def=""
msw_linalg_libs_def=""

# MPI
msw_mpi_enable_def="auto"
msw_mpi_incs_def=""
msw_mpi_libs_def=""

# pspBLAS
msw_psp_enable_def="no"
msw_psp_incs_def=""
msw_psp_libs_def=""

                    # ------------------------------------ #

#
# Check consistency of defaults
#

# Default triggers for optional features must be yes or no, and not empty
for msw_cfg_default in ${msw_cfg_switches}; do
  tmp_default=`eval echo \$\{msw_${msw_cfg_default}_enable_def\}`
  AC_MSG_CHECKING([default value of the ${msw_cfg_default} switch])
  AC_MSG_RESULT([${tmp_default}])
  if test "${tmp_default}" != "no" -a \
          "${tmp_default}" != "yes"; then
    AC_MSG_ERROR([invalid default value: msw_${msw_cfg_default}_enable_def = ${tmp_default}])
  fi
done

# Default triggers for external dependencies must be yes, no, or auto
for msw_cfg_default in ${msw_cfg_extdeps}; do
  tmp_default=`eval echo \$\{msw_${msw_cfg_default}_enable_def\}`
  AC_MSG_CHECKING([default detection behavior for ${msw_cfg_default}])
  AC_MSG_RESULT([${tmp_default}])
  if test "${tmp_default}" != "" -a \
          "${tmp_default}" != "auto" -a \
          "${tmp_default}" != "no" -a \
          "${tmp_default}" != "yes"; then
    AC_MSG_ERROR([invalid default value: msw_${msw_cfg_default}_enable_def = ${tmp_default}])
  fi
done

# Clean-up
unset tmp_default

                    # ------------------------------------ #

#
# Optional features
#

# Conditional compilation: --enable-* options (alphabetical order)
AC_ARG_ENABLE([conv],
  [AS_HELP_STRING([--enable-conv],
    [Enable real <-> complex conversion])],
  [msw_conv_enable="${enableval}"; msw_conv_type="yon"],
  [msw_conv_enable="${msw_conv_enable_def}"; msw_conv_type="def"])
AC_SUBST(enable_conv)
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],
    [Enable extra debug information])],
  [msw_debug_enable="${enableval}"; msw_debug_type="yon"],
  [msw_debug_enable="${msw_debug_enable_def}"; msw_debug_type="def"])
AC_SUBST(enable_debug)

                    # ------------------------------------ #

#
# External dependencies
#

# MPI (optional)
AC_ARG_WITH([mpi],
  [AS_HELP_STRING([--with-mpi],
    [Install prefix of MPI (e.g. /usr/local). The default behaviour is to detect whether the specified compilers properly support MPI and to fall back to serial mode if not. You may use --with-mpi without argument to force MPI detection, in which case detection failures will result in errors, and --without-mpi to disable MPI support completely.])],
  [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
      msw_mpi_enable="${withval}"
      msw_mpi_type="yon"
    else
      msw_mpi_enable="yes"
      msw_mpi_type="dir"
    fi],
  [msw_mpi_enable="${msw_mpi_enable_def}"; msw_mpi_type="def"])
AC_SUBST(with_mpi)

AC_ARG_VAR([MPICC], [MPI-aware C compiler (usually called 'mpicc'). Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPICC}"; then
  msw_mpi_enable="yes"
  msw_mpi_type="env"
fi
AC_ARG_VAR([MPIFC], [MPI-aware Fortran compiler (usually called 'mpif90'. Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPIFC}"; then
  msw_mpi_enable="yes"
  msw_mpi_type="env"
fi
AC_ARG_VAR([MPI_INCLUDES],
  [Include flags for MPI. Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPI_INCLUDES}"; then
  msw_mpi_enable="yes"
  msw_mpi_type="env"
  msw_mpi_incs="${MPI_INCLUDES}"
fi
AC_ARG_VAR([MPI_LIBS],
    [Library flags for MPI. Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPI_LIBS}"; then
  msw_mpi_enable="yes"
  msw_mpi_type="env"
  msw_mpi_libs="${MPI_LIBS}"
fi

# Linear algebra (optional)
AC_ARG_WITH([linalg],
  [AS_HELP_STRING([--with-linalg],
    [Install prefix of the linear algebra libraries (e.g. /usr/local). Linear algebra support will be enabled only if you use this option. You may use --with-linalg without argument if linear algebra libraries are available system-wide, and --without-linalg to disable linear algebra support completely.])],
  [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
      msw_linalg_enable="${withval}"
      msw_linalg_type="yon"
    else
      msw_linalg_enable="yes"
      msw_linalg_type="dir"
      msw_linalg_incs="-I${withval}/include"
      msw_linalg_libs="-L${withval}/lib"
    fi],
  [msw_linalg_enable="${msw_linalg_enable_def}"; msw_linalg_type="def";
    msw_linalg_incs=""; msw_linalg_libs=""])
AC_SUBST(with_linalg)

AC_ARG_VAR([LINALG_INCLUDES],
  [Include flags for linear algebra. Please note that this variable conflicts with the --with-linalg option.])
if test -n "${LINALG_INCLUDES}"; then
  msw_linalg_enable="yes"
  msw_linalg_type="env"
fi
AC_ARG_VAR([LINALG_LIBS],
    [Library flags for linear algebra. Please note that this variable conflicts with the --with-linalg option.])
if test -n "${LINALG_LIBS}"; then
  msw_linalg_enable="yes"
  msw_linalg_type="env"
fi

# pspBLAS (optional)
AC_ARG_WITH([psp],
  [AS_HELP_STRING([--with-psp],
    [Install prefix of the pspBLAS library (e.g. /usr/local). pspBLAS support will be enabled only if you use this option. You may use --with-psp without argument if pspBLAS is installed as a system library, and --without-psp to disable pspBLAS support completely.])],
  [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
      msw_psp_enable="${withval}"
      msw_psp_type="yon"
    else
      msw_psp_enable="yes"
      msw_psp_type="dir"
      msw_psp_incs="-I${withval}/include"
      msw_psp_libs="-L${withval}/lib -lpspBLAS"
    fi],
  [msw_psp_enable="${msw_psp_enable_def}"; msw_psp_type="def";
    msw_psp_incs=""; msw_psp_libs=""])
AC_SUBST(with_psp)

AC_ARG_VAR([PSP_INCLUDES],
  [Include flags for pspBLAS. Please note that this variable conflicts with the --with-psp option.])
if test -n "${PSP_INCLUDES}"; then
  msw_psp_enable="yes"
  msw_psp_type="env"
fi
AC_ARG_VAR([PSP_LIBS],
    [Library flags for pspBLAS. Please note that this variable conflicts with the --with-psp option.])
if test -n "${PSP_LIBS}"; then
  msw_psp_enable="yes"
  msw_psp_type="env"
fi


                    # ------------------------------------ #

#
# Check option consistency
#

# All --enable-* options must be yes or no
for msw_cfg_option in ${msw_cfg_switches}; do
  tmp_option=`eval echo \$\{enable_${msw_cfg_option}\}`
  if test "${tmp_option}" != "" -a \
          "${tmp_option}" != "no" -a \
          "${tmp_option}" != "yes"; then
    AC_MSG_ERROR([--enable-${msw_cfg_option} must be "yes" or "no"])
  fi
done

# Linear algebra must be enabled
#if test "${msw_linalg_enable}" != "yes"; then
#  AC_MSG_ERROR([linear algebra support has been disabled
#                  while it is mandatory])
#fi

# All MPI wrappers must be provided at once
tmp_mpi_vars="${MPICC}${MPIFC}"
if test -n "${tmp_mpi_vars}" -a \
        \( -z "${MPICC}" -o -z "${MPIFC}" \); then
  AC_MSG_ERROR([MPICC and MPIFC must both be set])
fi

# All pspBLAS environment variables must be provided at once
tmp_psp_vars="${PSP_INCLUDES}${PSP_LIBS}"
if test -n "${tmp_psp_vars}" -a \
        \( -z "${PSP_INCLUDES}" -o -z "${PSP_LIBS}" \); then
  AC_MSG_ERROR([PSP_INCLUDES and PSP_LIBS must both be set])
fi

# Environment variables conflict with --with-* options
# Note: for linear algebra, only LINALG_LIBS is relevant
if test -n "${tmp_mpi_vars}" -a -n "${with_mpi}"; then
  AC_MSG_ERROR([conflicting option settings
                  Please use MPICC + MPIFC or --with-mpi,
                  not both.])
fi
if test -n "${LINALG_LIBS}" -a -n "${with_linalg}"; then
  AC_MSG_ERROR([conflicting option settings
                  Please use LINALG_LIBS or --with-linalg,
                  not both.])
fi
if test -n "${tmp_psp_vars}" -a -n "${with_psp}"; then
  AC_MSG_ERROR([conflicting option settings
                  Please use PSP_INCLUDES + PSP_LIBS or --with-psp,
                  not both.])
fi

# When using environment variables, triggers must be set to yes
if test -n "${tmp_mpi_vars}"; then
  msw_mpi_enable="yes"
  msw_mpi_type="env"
fi
if test -n "${LINALG_LIBS}"; then
  msw_linalg_enable="yes"
  msw_linalg_type="env"
fi
if test -n "${tmp_psp_vars}"; then
  msw_psp_enable="yes"
  msw_psp_type="env"
fi

# pspBLAS requires MPI
if test "${msw_mpi_enable}" = "no" -a "${msw_psp_enable}" = "yes"; then
  AC_MSG_ERROR([pspBLAS requires MPI
                  Please adjust your parameters and re-run configure.])
fi

# Clean-up
unset tmp_option
unset tmp_mpi_vars
unset tmp_psp_vars

                    # ------------------------------------ #

#
# Substitute build-system variables
#

# Defaults
AC_SUBST(msw_conv_enable_def)
AC_SUBST(msw_debug_enable_def)
AC_SUBST(msw_mpi_enable_def)
AC_SUBST(msw_mpi_incs_def)
AC_SUBST(msw_mpi_libs_def)
AC_SUBST(msw_linalg_enable_def)
AC_SUBST(msw_linalg_incs_def)
AC_SUBST(msw_linalg_libs_def)
AC_SUBST(msw_psp_enable_def)
AC_SUBST(msw_psp_incs_def)
AC_SUBST(msw_psp_libs_def)

# Triggers
AC_SUBST(msw_conv_enable)
AC_SUBST(msw_debug_enable)
AC_SUBST(msw_mpi_enable)
AC_SUBST(msw_linalg_enable)
AC_SUBST(msw_psp_enable)

# Initialization types
AC_SUBST(msw_conv_type)
AC_SUBST(msw_debug_type)
AC_SUBST(msw_mpi_type)
AC_SUBST(msw_linalg_type)
AC_SUBST(msw_psp_type)

# Internal MPI parameters
AC_SUBST(msw_mpi_cc_set)
AC_SUBST(msw_mpi_cc_wrap)
AC_SUBST(msw_mpi_fc_set)
AC_SUBST(msw_mpi_fc_wrap)

# Build flags for external dependencies
AC_SUBST(msw_mpi_incs)
AC_SUBST(msw_mpi_libs)
AC_SUBST(msw_linalg_incs)
AC_SUBST(msw_linalg_libs)
AC_SUBST(msw_linalg_has_lapack)
AC_SUBST(msw_linalg_has_scalapack)
AC_SUBST(msw_psp_incs)
AC_SUBST(msw_psp_libs)

                    # ------------------------------------ #

#
# Architecture setup
#

# Improve defaults using information from EasyBuild
if test "${msw_mpi_enable}" != "no" -a \
        -n "${EBROOTOPENMPI}" -a \
        "${msw_mpi_type}" = "def"; then
  if test "${CC}" = "" -a "${FC}" = ""; then
    CC="mpicc"
    FC="mpifort"
  else
    MPICC="mpicc"
    MPIFC="mpifort"
  fi
fi

# Init MPI parameters, possibly setting and/or wrapping compilers
MSW_MPI_INIT

                    # ------------------------------------ #

#
# C language support
#

# Init C compiler and preprocessor
AC_PROG_CC
AC_PROG_CPP

# C compiler peculiarities (for Libtool)
AM_PROG_CC_C_O

                    # ------------------------------------ #

#
# Fortran language support
#

# Look for the Fortran compiler
if test "${FC}" != "" -a ! -x "${FC}"; then
  msw_fc_probe=`echo "${FC}" | sed -e 's/ .*//'`
  if test ! -x "${msw_fc_probe}"; then
    AC_PATH_PROG([msw_fc_path],[${msw_fc_probe}])
    if test "${msw_fc_path}" = ""; then
      AC_MSG_ERROR([could not run Fortran compiler "${FC}"])
    fi
  fi
fi
AC_PROG_FC

# Fail if no Fortran compiler is available
if test "${FC}" = ""; then
  AC_MSG_ERROR([no Fortran compiler available])
fi

# Set default file extensions
MSW_FC_EXTENSIONS

# Get module file case
MSW_FC_MOD_CASE

                    # ------------------------------------ #

#
# Libtool configuration
#

# Init Libtool (must be done once compilers are fully set)
LT_INIT
LT_PREREQ([2.4.2])
LT_LANG([Fortran])
LTOBJEXT="lo"
AC_SUBST(LTOBJEXT)

                    # ------------------------------------ #

#
# Multicore architecture support
#

# Look for MPI
if test "${msw_mpi_enable}" != "no"; then
  MSW_MPI_DETECT
  if test "${msw_mpi_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI], 1, [Define to 1 to enable MPI support.])
    msw_mpi_enable="yes"
  else
    if test "${msw_mpi_enable}" = "auto"; then
      AC_MSG_NOTICE([disabling MPI support])
      msw_mpi_enable="no"
    else
      AC_MSG_FAILURE([MPI support is broken - please check your configuration])
    fi
  fi
fi

# From now on, the MPI trigger must be yes or no
if test "${msw_mpi_enable}" != "no" -a "${msw_mpi_enable}" != "yes"; then
  AC_MSG_ERROR([unexpected MPI trigger value: '${msw_mpi_enable}'])
fi

# Let the test suite know about mpirun
TESTS_CMD="${MPIRUN}"
AC_SUBST(TESTS_CMD)

                    # ------------------------------------ #

#
# External dependencies
#

# Set linear algebra includes
if test "${msw_linalg_enable}" != "no"; then
  if test "${msw_linalg_type}" != "dir"; then
    if test -n "${LINALG_INCLUDES}"; then
      msw_linalg_incs="${LINALG_INCLUDES}"
    else
      msw_linalg_incs="${msw_linalg_incs_def}"
    fi
  fi
fi

# Set linear algebra libraries
if test "${msw_linalg_enable}" != "no"; then
  if test "${msw_linalg_type}" != "dir"; then
    if test -n "${LINALG_LIBS}"; then
      msw_linalg_libs="${LINALG_LIBS}"
    else
      msw_linalg_libs="${msw_linalg_libs_def}"
      if test "${msw_linalg_type}" = "def" -o "${msw_linalg_type}" = "yon"; then
        if test -n "${EBROOTOPENBLAS}"; then
          AC_MSG_NOTICE([using EasyBuild configuration to set linear algebra libraries])
          msw_linalg_type="ebf"
          msw_linalg_libs="-lopenblas"
        fi
        if test -n "${EBROOTSCALAPACK}" -a "${msw_mpi_enable}" = "yes"; then
          msw_linalg_libs="-lscalapack ${msw_linalg_libs}"
        fi
      fi
    fi
  fi
fi

# Look for linear algebra
AC_MSG_CHECKING([whether to look for linear algebra])
AC_MSG_RESULT([${msw_linalg_enable}])
if test "${msw_linalg_enable}" = "yes"; then
  AC_MSG_CHECKING([how linear algebra parameters have been set])
  AC_MSG_RESULT([${msw_linalg_type}])
  AC_MSG_CHECKING([for linear algebra include flags])
  AC_MSG_RESULT([${msw_linalg_incs}])
  AC_MSG_CHECKING([for linear algebra library flags])
  AC_MSG_RESULT([${msw_linalg_libs}])
  MSW_LINALG_DETECT
  if test "${msw_linalg_ok}" != "yes" -a "${msw_linalg_type}" = "def"; then
    AC_MSG_NOTICE([retrying with generic linear algebra library flags (Netlib)])
    LIBS="-lscalapack -lblacs -lblacsCinit -lblacsF77init -llapack -lblas ${LIBS}"
    MSW_LINALG_DETECT
  fi
  if test "${msw_linalg_ok}" = "yes"; then
    CPPFLAGS="${CPPFLAGS} ${msw_linalg_incs}"
    LIBS="${msw_linalg_libs} ${LIBS}"
  else
    AC_MSG_FAILURE([linear algebra is missing or incomplete])
  fi
fi

                              # ---------------- #

# Check that MPI is still enabled
if test "${msw_mpi_enable}" = "no" -a "${msw_psp_enable}" = "yes"; then
  AC_MSG_ERROR([pspBLAS requires MPI
                  Please adjust your parameters and re-run configure.])
fi

# Set pspBLAS includes
if test "${msw_psp_enable}" != "no"; then
  if test "${msw_psp_type}" != "dir"; then
    if test -n "${PSP_INCLUDES}"; then
      msw_psp_incs="${PSP_INCLUDES}"
    else
      msw_psp_incs="${msw_psp_incs_def}"
      if test -d "../../tmp-pspblas"; then
        AC_MSG_NOTICE([using OMM Bundle configuration to set PspBLAS includes])
        msw_psp_type="omm"
        msw_psp_incs="-I../../tmp-pspblas/include"
      fi
    fi
  fi
fi

# Set pspBLAS libraries
if test "${msw_psp_enable}" != "no"; then
  if test "${msw_psp_type}" != "dir"; then
    if test -n "${PSP_LIBS}"; then
      msw_psp_libs="${PSP_LIBS}"
    else
      msw_psp_libs="${msw_psp_libs_def}"
      if test -d "../../tmp-pspblas"; then
        AC_MSG_NOTICE([using OMM Bundle configuration to set PspBLAS libraries])
        msw_psp_type="omm"
        msw_psp_libs="-L../../tmp-pspblas/lib -lpspBLAS"
      fi
    fi
  fi
fi

# Look for pspBLAS
AC_MSG_CHECKING([whether to look for pspBLAS])
AC_MSG_RESULT([${msw_psp_enable}])
if test "${msw_psp_enable}" = "yes"; then
  AC_MSG_CHECKING([how pspBLAS parameters have been set])
  AC_MSG_RESULT([${msw_psp_type}])
  AC_MSG_CHECKING([for pspBLAS include flags])
  AC_MSG_RESULT([${msw_psp_incs}])
  AC_MSG_CHECKING([for pspBLAS library flags])
  AC_MSG_RESULT([${msw_psp_libs}])
  MSW_PSP_DETECT
  if test "${msw_psp_ok}" = "yes"; then
    CPPFLAGS="${CPPFLAGS} ${msw_psp_incs}"
    LIBS="${msw_psp_libs} ${LIBS}"
  else
    AC_MSG_FAILURE([pspBLAS is missing or incomplete])
  fi
fi

                    # ------------------------------------ #

#
# Adjust build parameters according to configure options
#

# Reset option counter
msw_cpp_options=""

# ISO C Bindings are always available nowadays
AC_DEFINE([CBIND], 1, [Define to 1 to use ISO C Bindings])
msw_cpp_options="${msw_cpp_options} CBIND"

# Conditional compilation (alphabetical order)
if test "${msw_conv_enable}" = "yes"; then
  AC_DEFINE([CONV], 1,
    [Define to 1 if you want to use real <-> complex conversion.])
  msw_cpp_options="${msw_cpp_options} CONV"
fi
if test "${msw_debug_enable}" = "yes"; then
  AC_DEFINE([DEBUG_MODE], 1,
    [Define to 1 if you want to enable enhanced debugging information.])
  msw_cpp_options="${msw_cpp_options} DEBUG_MODE"
fi

# External dependencies (alphabetical order)
if test "${msw_linalg_has_lapack}" = "yes"; then
  AC_DEFINE([HAVE_LAPACK], 1, [Define to 1 if you have LAPACK libraries available.])
  msw_cpp_options="${msw_cpp_options} HAVE_LAPACK"
fi
if test "${msw_linalg_has_scalapack}" = "yes"; then
  AC_DEFINE([HAVE_SCALAPACK], 1, [Define to 1 if you have ScaLAPACK libraries available.])
  msw_cpp_options="${msw_cpp_options} HAVE_SCALAPACK"
fi
if test "${msw_psp_ok}" = "yes"; then
  AC_DEFINE([HAVE_PSPBLAS], 1, [Define to 1 if you have the pspBLAS library available.])
  msw_cpp_options="${msw_cpp_options} HAVE_PSPBLAS"
fi

                    # ------------------------------------ #

# Write down YAML configuration now (allows full report with deferred errors)
AC_OUTPUT([config/data/MatrixSwitch-config.yaml])

                    # ------------------------------------ #

# Report configuration
AC_MSG_NOTICE([])
AC_MSG_NOTICE([Final build parameters])
AC_MSG_NOTICE([----------------------])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([TSTAMP     = ${msw_timestamp}])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([CONVERT  = ${msw_conv_enable} (init: ${msw_conv_type})])
AC_MSG_NOTICE([DEBUG    = ${msw_debug_enable} (init: ${msw_debug_type})])
AC_MSG_NOTICE([MPI      = ${msw_mpi_enable} (init: ${msw_mpi_type})])
AC_MSG_NOTICE([LINALG   = ${msw_linalg_enable} (init: ${msw_linalg_type})])
AC_MSG_NOTICE([PSP      = ${msw_psp_enable} (init: ${msw_psp_type})])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([CPP      = ${CPP}])
AC_MSG_NOTICE([CPPFLAGS = ${CPPFLAGS}])
AC_MSG_NOTICE([FC       = ${CC}])
AC_MSG_NOTICE([MPIFC    = ${MPICC}])
AC_MSG_NOTICE([FC       = ${FC}])
AC_MSG_NOTICE([MPIFC    = ${MPIFC}])
AC_MSG_NOTICE([FCFLAGS  = ${FCFLAGS}])
AC_MSG_NOTICE([LDFLAGS  = ${LDFLAGS}])
AC_MSG_NOTICE([LIBS     = ${LIBS}])
AC_MSG_NOTICE([OPTIONS  = ${msw_cpp_options}])
AC_MSG_NOTICE([])

                    # ------------------------------------ #

#
# Output configuration
#

AC_CONFIG_FILES([
  config/data/MatrixSwitch.pc
  Makefile
  doc/Makefile
  src/Makefile
  examples/Makefile
])

AC_OUTPUT
