#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([DSM Kernel], [1.0.0], [nobuaki@ipgp.fr])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([README.md])
AC_CONFIG_HEADERS([./config.h])


# Checks for programs.


AC_PROG_CC

AC_PROG_FC
#AM_INIT_AUTOMAKE([-O4 -check nobounds -xAVX -ftz -assume buffered_io -assume byterecl -vec-report3 -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel])

#AC_SUBST([AM_MAKEFLAGS],[-O4 -check nobounds -xAVX -ftz -assume buffered_io -assume byterecl -vec-report3 -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel])

# Checks for programs.

# a courtesy to the installed base of users

#AC_FC_MODULE_FLAG([-check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -assume byterecl])


# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/SGTpsv/Makefile
                 src/SGTsh/Makefile
		 src/KernelMaker/Makefile])


AC_OUTPUT
