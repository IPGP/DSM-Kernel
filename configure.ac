#
#
#    Process this file with autoconf to produce the configure script
# 
#       You must have recent versions of Autoconf and Automake installed.
#
#
#      Nobuaki Fuji, Institut de Physique du Globe de Paris
#
#
#                    NF referred to configure.ac of specfem2d


AC_PREREQ(2.61)
AC_INIT([DSM Kernel],[1.0.0],[nobuaki@ipgp.fr])
AC_CONFIG_SRCDIR([README.md])
AC_CONFIG_HEADER([etc/config_calcul/config.h])
AC_CONFIG_MACRO_DIR([m4])

# 'configure'options

AC_ARG_WITH([mpi],
    [AC_HELP_STRING([--with--mpi],
        [build parallel version @<:@default=yes@:>@])],
    [want_cuda="$withval"],
    [want_cuda=no])
AM_CONDITIONAL([COND_CUDA],[test "$want_cuda" !=no])
AM_CONDITIONAL([COND_CUDA5],[test "$want_cuda" = cuda5])

# debugging

AC_ARG_ENABLE([debug],
    [AC_HELP_STRING([--enable-debug],
        [build with debugging options enabled @<:@default=no@:>@])],
    [want_debug="$enableval"],
    [want_debug=no])
AM_CONDITIONAL([COND_DEBUG], [test x"$want_debug" != xno])
export COND_DEBUG_FALSE
export COND_DEBUG_TRUE



AC_ARG_ENABLE([double-precision],
    [AC_HELP_STRING([--enable-double-precision],
        [solver in double precision @<:@default=no@:>@])],
    [want_double_precision="$enableval"],
    [want_double_precision=no])
if test "$want_double_precision" = no; then
    CUSTOM_REAL=SIZE_REAL
    CUSTOM_MPI_TYPE=MPI_REAL
else
    CUSTOM_REAL=SIZE_DOUBLE
    CUSTOM_MPI_TYPE=MPI_DOUBLE_PRECISION
fi
AC_SUBST([CUSTOM_REAL])
AC_SUBST([CUSTOM_MPI_TYPE])

