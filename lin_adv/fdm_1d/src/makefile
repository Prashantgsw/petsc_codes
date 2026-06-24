.SUFFIXES:
FC      = mpif90
GIT_VERSION = $(shell git describe --abbrev=8 --dirty --always)
FFLAGS  = -O3 -I$(PETSC_DIR)/include \
		-DVERSION=\"$(GIT_VERSION)\"
LDFLAGS = -lpetsc -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/usr/lib64

SRC  = $(wildcard *.F90)
OBJ  = $(patsubst %.F90,%.o,$(SRC))
MOD  = double.mod variables.mod read_write.mod auxillary_conditions.mod \
	fdm.mod solve.mod

TARGET  = exe

$(TARGET): $(OBJ) $(MOD)
	$(FC) -o $(TARGET) $(OBJ) $(LDFLAGS)

double.mod: double.F90
	$(FC) -c $(FFLAGS) double.F90

variables.mod: variables.F90 double.mod
	$(FC) -c $(FFLAGS) variables.F90

read_write.mod: read_write.F90 variables.mod double.mod
	$(FC) -c $(FFLAGS) read_write.F90

auxillary_conditions.mod: auxillary_conditions.F90 variables.mod double.mod
	$(FC) -c $(FFLAGS) auxillary_conditions.F90

fdm.mod: fdm.F90 variables.mod double.mod
	$(FC) -c $(FFLAGS) fdm.F90

solve.mod: solve.F90 variables.mod double.mod
	$(FC) -c $(FFLAGS) solve.F90

%.o: %.F90 $(MOD)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.mod *.o $(OBJ) $(TARGET) *.dat verify_convergence