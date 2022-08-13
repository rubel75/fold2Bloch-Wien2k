# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Intel fortran compiler
FC = ifort
FLFLAGS = # none
#FCFLAGS = -free # performance
FCFLAGS = -free -g -traceback -check all -debug all # debug

# gfortran compiler
#FC = gfortran
#FLFLAGS = # none
#FCFLAGS = -ffree-form # performance
#FCFLAGS = -ffree-form -g -fbacktrace # debug


# ~~~ Do not edit after that line ~~~

PROGRAM = fold2Bloch

# source files and objects
SRCS = util.o $(patsubst %.F, %.o, $(wildcard *.F)) \
    $(patsubst %.F90, %.o, $(wildcard *.F90))

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FLFLAGS) $(FLINK) -o $@ $^

util.o: util.F90
	$(FC) $(FCFLAGS) $(FOPT) -c util.F90

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
