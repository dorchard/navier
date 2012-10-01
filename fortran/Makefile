all: 	
	gfortran -c helpers.f90
	gfortran -c init.f90
	gfortran -c boundary.f90
	gfortran -o navier main.f90 helpers.f90 init.f90 boundary.f90