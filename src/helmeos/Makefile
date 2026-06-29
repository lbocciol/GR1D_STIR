# Standalone build for the Helmholtz EOS unit tests.
#
#   make                  # build temp_inversion_test
#   make test             # build and run the pure-Helmholtz (rho,T,Ye) sweep
#   make helm_test        # build the composite nuc_eos+Helmholtz blend test
#   make helm_test_run    # build and run the composite blend test
#   make clean
#
# temp_inversion_test exercises only wlHelmholtzEOS (self-contained, builds
# with plain gfortran).  helm_test exercises the nuc_eos + Helmholtz blend, so
# it links nuc_eos.a and therefore needs mpifort + HDF5 (set HDF5_DIR), the
# same toolchain as the main GR1D build.

FC      ?= gfortran
FFLAGS  ?= -O2

# composite-test toolchain (matches make.inc); override on the command line
HELMFC    ?= mpifort
NUCEOS     = ../nuc_eos/nuc_eos.a
HDF5DIR   ?= $(HDF5_DIR)
HDF5INCS  ?= -I$(HDF5DIR)/include
HDF5LIBS  ?= -L$(HDF5DIR)/lib -lhdf5 -lhdf5_fortran -lhdf5 -lz
# nuc_eos table for the composite test (override: make helm_test_run TABLE=...)
TABLE     ?= ../../Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5

EOSOBJ  = helmholtzEOS.o

all: temp_inversion_test

$(EOSOBJ): helmholtzEOS.F90
	$(FC) $(FFLAGS) -c $< -o $@

# helmholtzEOS.o for the composite test must be built with the same compiler
# (mpifort) as nuc_eos.a so the .mod / objects are compatible.
helmholtzEOS_helm.o: helmholtzEOS.F90
	$(HELMFC) $(FFLAGS) -c $< -o $@

temp_inversion_test: temp_inversion_test.f90 $(EOSOBJ)
	$(FC) $(FFLAGS) $^ -o $@

$(NUCEOS):
	$(MAKE) -C ../nuc_eos

helm_test: helm_test.f90 helmholtzEOS_helm.o $(NUCEOS)
	$(HELMFC) $(FFLAGS) -I../nuc_eos $(HDF5INCS) \
		helm_test.f90 helmholtzEOS_helm.o $(NUCEOS) $(HDF5LIBS) -o helm_test

# run from this directory so helm_table.dat is found
test: temp_inversion_test
	./temp_inversion_test

helm_test_run: helm_test
	./helm_test $(TABLE)

clean:
	rm -f temp_inversion_test helm_test $(EOSOBJ) helmholtzEOS_helm.o wlhelmholtzeos.mod

.PHONY: all test helm_test_run clean
