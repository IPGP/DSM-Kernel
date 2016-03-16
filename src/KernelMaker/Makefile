BIN = ../../bin
program = $(BIN)/mpiKernelMaker
include ../../etc/config-calcul/config.h
objs =   bwfilt.o module.o others.o fourier.o calculateK.o KernelMaker-0.1.0.o

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(option)
.f90.o:
	$(FC) -o $@ -c $< $(option)

.PHONY: clean
clean:
	rm $(program) $(objs)  
