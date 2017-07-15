#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])


# Checks for programs.
AC_PROG_CC
AC_PROG_FC([mpiifort mpif90])
AC_SUBST([AC_SAVE_FCFLAGS],[-O4 -check nobounds -xAVX -qopt-report -ftz -assume buffered_io -assume byterecl -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 examples/multiple_kernels/Makefile
                 src/KernelMaker/Makefile
                 src/Makefile
                 src/ModelDrawer/Makefile
                 src/ReadKernel/Makefile
                 src/SGTpsv/Makefile
                 src/SGTsh/Makefile])
AC_OUTPUT
