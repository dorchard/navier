FC:=gfortran
GFLAGS:=-O2
OUTPUT:=navier
SOURCES:=helpers_mod.f90 init_mod.f90 boundary_mod.f90 simulation_mod.f90 output_mod.f90 navier.f90
OBJS:=$(patsubst %.f90, %.o, $(SOURCES))

.PHONY: all 
all: $(MODS) $(OBJS) output
	$(FC) -o $(OUTPUT) $(OBJS)

%.o: %.f90
	$(FC) $(GFLAGS) -c $< 

output: 
	mkdir output

.PHONY: clean
clean:
	rm -f *.o *.mod $(OUTPUT) 

