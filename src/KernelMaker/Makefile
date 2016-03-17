BIN = ../../bin/
include ../../etc/config_calcul/config.h
program = $(BIN)/mpiKernelMaker
objs = bwfilt.o module.o others.o fourier.o calculateK.o KernelMaker-1.0.0.o

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
clean:
	rm $(program) $(objs)  *.lst *.o
