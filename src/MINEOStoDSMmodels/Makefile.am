FC=gfortran
#AM_FCFLAGS = -L /cm/shared/modulefiles/lapack/intel/64/3.4.2 
#AM_FCFLAGS = -llapack --ffixed-line-length-none
# this flag for malbec IPGP
AM_FCFLAGS = -O2 -ftz -fpe0 -traceback -ftrapuv -implicitnone -check all -assume byterecl -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -L/cm/shared/apps/lapack/intel/64/3.4.2 -llapack -L/cm/shared/apps/blas/intel/1/lib64 -lblas
bin_PROGRAMS = $(top_builddir)/bin/MINEOStoDSMmodels  
__top_builddir__bin_MINEOStoDSMmodels_SOURCES =   modules.f90 others.f90 sub.inverse.f90 MINEOStoDSMmodels.f90 
