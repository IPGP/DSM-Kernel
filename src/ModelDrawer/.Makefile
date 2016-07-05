program = ../../bin/ModelDrawer
objs =    module.o others.o ModelDrawer-0.1.o
#FC=mpif90
FC=ifort
option = -O2 -ftz -fpe0 -traceback -ftrapuv -vec-report0 -implicitnone -check all -assume byterecl -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage  

#option = -check nobounds -O3 -ftz -fpe0 -vec-report0 -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage -assume byterecl 
# -L//usr/local/ansys121/v121/ansys/intelmpi/linx64/lib64/  
# -warn unused has been removed 
.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(option)
.f90.o:
	$(FC) -o $@ -c $< $(option)

.PHONY: clean
clean:
	rm $(program) $(objs) *.mod 
