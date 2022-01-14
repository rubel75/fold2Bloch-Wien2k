# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Intel fortran compiler
FC = ifort
FLFLAGS = # none
FCFLAGS = -free #-g -traceback -check all -debug all

# gfortran compiler
#FC = gfortran
#FLFLAGS = # none
#FCFLAGS = -ffree-form


# ~~~ Do not edit after that line ~~~

PROGRAM = fold2Bloch

# source files and objects
SRCS = $(patsubst %.F, %.o, $(wildcard *.F)) \
    $(patsubst %.F90, %.o, $(wildcard *.F90))

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FLFLAGS) $(FLINK) -o $@ $^

%.o: %.F
	$(FC) $(FCFLAGS) $(FOPT) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) $(FOPT) -c $<

#%.mod: %.h
#	$(FC) $(FCFLAGS) -o $@ $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)
