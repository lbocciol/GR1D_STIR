USE_SCREENING ?= TRUE
DEFINES :=
ifeq ($(USE_SCREENING),TRUE)
   DEFINES += -DSCREENING
endif

FC  := gfortran
CXX := g++
FFLAGS   := -O2 -I./pynucnet
CXXFLAGS := -O2 -std=c++20 -I./pynucnet $(DEFINES)

# TODO(build): replace with -lmkl_rt or point to a local LAPACK build if needed
LIBS := -llapack -lblas -lstdc++

.PHONY: all clean

all: burn_test

pynucnet/wrapper.o: pynucnet/wrapper.cpp pynucnet/*.H
	$(CXX) $(CXXFLAGS) -c $< -o $@

# pynet.mod lands in pynucnet/ (already built); recompile only if .f90 changes
pynucnet/fortran_interface.o: pynucnet/fortran_interface.f90
	$(FC) $(FFLAGS) -J./pynucnet -c $< -o $@

burn.o: burn.F90 pynucnet/fortran_interface.o
	$(FC) $(FFLAGS) -c $< -o $@

burn_test.o: burn_test.f90 burn.o
	$(FC) $(FFLAGS) -c $< -o $@

burn_test: burn_test.o burn.o pynucnet/fortran_interface.o pynucnet/wrapper.o
	$(FC) -o $@ $^ $(LIBS)

clean:
	rm -f burn_test burn.o burn_test.o
