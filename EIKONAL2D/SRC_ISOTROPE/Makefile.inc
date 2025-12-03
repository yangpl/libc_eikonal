##################################################################
CC       = gcc
F77      = gfortran -Wunused -Wuninitialized -D -g -fbounds-check -Wall -fbacktrace -finit-real=nan
F90      = gfortran -Wunused -Wuninitialized -D -g -fbounds-check -Wall -fbacktrace -finit-real=nan

F77      = gfortran -Wunused -Wuninitialized -D -g -fbounds-check -fbacktrace -finit-real=nan
F90      = gfortran -Wunused -Wuninitialized -D -g -fbounds-check -fbacktrace 

OPTF = -O3
OPTC =  

%.o : %.mod 

.SUFFIXES: .f .o
.f.o:
	$(F77) $(OPTF) -c $*.f

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(OPTF) -c $*.f90

.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) -c $*.c

