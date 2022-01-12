# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Intel fortran compiler (gfortran options are in brackets)
FC = ifort # (gfortran)
FCFLAGS = #-free -g -traceback -check all -debug all # none
FLFLAGS = -free #-g -traceback -check all -debug all # (-ffree-form)

# gfortran compiler


# ~~~ Do not edit after that line ~~~

PROGRAM = fold2Bloch

# source files and objects
SRCS = $(patsubst %.F, %.o, $(wildcard *.F)) \
    $(patsubst %.F90, %.o, $(wildcard *.F90))

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) $(FLINK) -o $@ $^

%.o: %.F
	$(FC) $(FLFLAGS) $(FOPT) -c $<

%.o: %.F90
	$(FC) $(FLFLAGS) $(FOPT) -c $<

#%.mod: %.h
#	$(FC) $(FLFLAGS) -o $@ $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)
