BIN = ../../bin
program = $(BIN)/mpiSGTpsv 
include ../../etc/config_calcul/config.h
objs =   trial.o rk3.o glu2.o dcsymbdl.o dcsymbdl3.o others.o calmat.o SGTpsv-0.3.3-MPI.o

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
clean:
	rm $(program) $(objs) *.mod *.optrpt 
