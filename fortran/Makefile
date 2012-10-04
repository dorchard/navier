GFLAGS = -O2

all: 	
	gfortran $(GFLAGS) -c helpers.f90
	gfortran $(GFLAGS) -c init.f90
	gfortran $(GFLAGS) -c boundary.f90
	gfortran $(GFLAGS) -c simulation.f90
	gfortran $(GFLAGS) -c output.f90
	gfortran $(GFLAGS) -o navier main.f90 helpers.f90 init.f90 boundary.f90 simulation.f90 output.f90