include ./Makefile.inc
EXEC= ../BIN

all: eikonal_isotrope stencil_isotropic_test


clean:
	rm *.o *.mod *~ 

######################## NETLIB TOMS package 954

SUB = Polynomial2RootSolvers.o neighbors.o stencil_triangle.o local_solver.o eikonal_solver.o read_inputs.o

############################################################################################################

############################################################################################################
stencil_test: $(SUB) stencil_test.o
	$(F90) $(OPTF) $(SUB) stencil_test.o -o $(EXEC)/stencil_test

stencil_isotropic_test: $(SUB) stencil_isotropic_test.o
	$(F90) $(OPTF) $(SUB) stencil_isotropic_test.o -o $(EXEC)/stencil_isotropic_test



eikonal_isotrope: $(SUB) EIKONAL.o
	$(F90) $(OPTF) $(SUB) EIKONAL.o -o $(EXEC)/eikonal_isotrope

