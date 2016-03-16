BIN = ../../bin
program = $(BIN)/mpiSGTsh
include ../../etc/config-calcul/config.h
objs =    glu2.o  dclib.o dclib3.o others.o trial.o calmat.o SGTsh-0.1.0-MPI.o
.SUFFIXES: .o .f90

$(program): $(objs)
	        $(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	        $(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
	clean:
	        rm $(program) $(objs)  


