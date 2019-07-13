# Makefile for D Wave internship code

#compiler
gfortErr = $(shell which gfortran >/dev/null; echo $$?)
ifeq "$(gfortErr)" "0"
F90C = gfortran
#flags
F90FLAGS = -O3 -fopenmp -ffast-math -fmax-errors=4 -fbounds-check
endif

OBJECTS = montecarlo.o problems.o precision.o graphs.o model.o
DRIVER = driver.F90

default: sampler

sampler: $(OBJECTS) $(DRIVER)
	$(F90C) $(F90FLAGS) $(OBJECTS) $(DRIVER) -o sampler

clean :
	rm -f *.o *.a *.d *.mod

# OBJECTS instructions

precision.o : precision.f90
	$(F90C) $(F90FLAGS) -c precision.f90 -o precision.o

montecarlo.o: precision.o montecarlo.f90
	$(F90C) $(F90FLAGS) -c montecarlo.f90 -o montecarlo.o

problems.o : precision.o problems.f90
	$(F90C) $(F90FLAGS) -c problems.f90 -o problems.o

graphs.o : precision.o graphs.f90
	$(F90C) $(F90FLAGS) -c graphs.f90 -o graphs.o

model.o : precision.o model.f90
	$(F90C) $(F90FLAGS) -c model.f90 -o model.o

