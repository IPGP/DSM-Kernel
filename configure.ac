#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([DSM-Kernel], [1.0.1], [nobuaki@ipgp.fr])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])


# Checks for programs.
# AC_PROG_CC
#AC_PROG_FC([mpiifort mpif90])

if [ [  -z "${FCFLAGS+set}" ] ]
 then
 FCFLAGS='-O4 -check nobounds -xAVX -qopt-report -ftz -assume buffered_io -assume byterecl -implicitnone -warn truncated_source -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel' 
 FCFLAGS='-O4 -check nobounds -xAVX -qopt-report -ftz -assume buffered_io -assume byterecl -implicitnone -warn truncated_source -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel'
 #FCFLAGS='-std=gnu -fimplicit-none -frange-check -O2 -pedantic -pedantic-errors -Waliasing -Wampersand -Wline-truncation -Wsurprising -Wunderflow -fbounds-check -assume byterecl'
fi 
#AC_SUBST([ac_test_FCFLAGS],[set])
#AC_SUBST([ac_save_FCFLAGS],['-O4 -check nobounds -xAVX -qopt-report -ftz -assume buffered_io -assume byterecl -implicitnone -warn truncated_source -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel'])


AC_PROG_FC([mpiifort mpif90])
#AC_PROG_FC([mpif90 mpiifort])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.


AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/ModelDrawer/Makefile
                 src/ReadKernel/Makefile
                 src/SGTpsv/Makefile
                 src/SGTsh/Makefile
                 src/KernelMaker/Makefile
                 examples/global_kernels/database.inf
                 examples/global_kernels/kernel.inf
                 examples/innercorekernels/database.inf
                 examples/innercorekernels/kernel.inf
                 examples/multiple_kernels/database.inf
                 examples/multiple_kernels/kernel_config/ABC.20200101.Z.inf
                 examples/single_kernel/database.inf
                 examples/single_kernel/kernel.inf
])


AC_OUTPUT


#AC_CONFIG_FILES([Makefile
#                 examples/Makefile
#                 examples/multiple_kernels/Makefile
#                 src/KernelMaker/Makefile
#                 src/Makefile
#                 src/ModelDrawer/Makefile
#                 src/ReadKernel/Makefile
#                 src/SGTpsv/Makefile
#                 src/SGTsh/Makefile])
#AC_OUTPUT
