include ../../etc/config_calcul/config.h
bin_PROGRAMS = $(top_builddir)/bin/mpiSGTsh 
__top_builddir__bin_mpiSGTsh_SOURCES =   glu2.f90 dclib.f90 dclib3.f90 others.f90 trial.f90 calmat.f90 SGTsh-0.1.0-MPI.f90
#clean-local: clean-local-SGTsh
.PHONY: clean-local-SGTsh
clean-local-SGTpsv:
	rm $(bin_PROGRAMS) *.o *.mod *.optrpt $(top_builddir)/bin/*.optrpt
