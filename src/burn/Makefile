include ../../make.inc

# Compiler settings (use make.inc defaults, but allow override)
F90 ?= gfortran
F90FLAGS ?= -g -O2
CXX ?= g++
CXXFLAGS ?= -O2 -std=c++20 -DSCREENING

# Screening can be toggled: USE_SCREENING=FALSE to disable
USE_SCREENING ?= TRUE
ifeq ($(USE_SCREENING),FALSE)
    # Remove -DSCREENING from CXXFLAGS if present
    CXXFLAGS := $(filter-out -DSCREENING,$(CXXFLAGS))
endif

# Debug mode
DEBUG ?= 0
ifeq ($(DEBUG),1)
    F90FLAGS = -g -O0 -Wall -Wextra -fcheck=all -fsanitize=address,undefined
    CXXFLAGS := -g -O0 -std=c++20 $(filter -DSCREENING,$(CXXFLAGS))
endif

# Include paths
MODINC ?= -I./
F90FLAGS += $(MODINC)

# LAPACK
LAPACKDIR ?=
LAPACKLIBS ?= -llapack -lblas

# Source files and objects
BURN_SOURCES = burn.F90 composition.F90 nse.F90
BURN_OBJECTS = burn.o composition.o nse.o

# pynucnet network objects
NETWORK_SUBDIR = pynucnet
NETWORK_OBJECTS = $(NETWORK_SUBDIR)/fortran_interface.o $(NETWORK_SUBDIR)/wrapper.o

# Test executable
TEST_SOURCE = burn_test.f90
TEST_OBJECT = burn_test.o
TEST_EXE = burn_test

# Default target
.PHONY: all
all: $(TEST_EXE)

# Build pynucnet objects (with header dependencies)
$(NETWORK_SUBDIR)/%.o: $(NETWORK_SUBDIR)/%.cpp $(wildcard $(NETWORK_SUBDIR)/*.H)
	$(CXX) $(CXXFLAGS) -I. -I$(NETWORK_SUBDIR) -c $< -o $@

$(NETWORK_SUBDIR)/%.o: $(NETWORK_SUBDIR)/%.f90
	$(F90) $(F90FLAGS) -c $< -o $@

# Fortran compilation rules for burn modules
%.o: %.F90
	$(F90) $(F90FLAGS) -c $< -o $@

%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

# Module build order (composition uses burn+pynet; nse uses composition)
burn.o: $(NETWORK_SUBDIR)/fortran_interface.o
composition.o: burn.o $(NETWORK_SUBDIR)/fortran_interface.o
nse.o: composition.o
burn_test.o: $(BURN_OBJECTS)

# Executable (depends on network objects)
$(TEST_EXE): $(NETWORK_OBJECTS) $(BURN_OBJECTS) $(TEST_OBJECT)
	$(F90) $(F90FLAGS) -o $@ $^ $(LAPACKLIBS) -lstdc++

# Clean up
.PHONY: clean
clean:
	rm -f *.o *.mod $(TEST_EXE)
	rm -f pynucnet/*.o

.PHONY: distclean
distclean: clean
	rm -f *.mod
	rm -f pynucnet/*.mod

# Help
.PHONY: help
help:
	@echo "burn_test Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all      - build burn_test executable (default)"
	@echo "  clean    - remove object files and executable"
	@echo "  distclean - also remove .mod files"
	@echo "  help     - show this message"
	@echo ""
	@echo "Variables:"
	@echo "  DEBUG=0/1        - enable debug flags (default: 0)"
	@echo "  USE_SCREENING=TRUE/FALSE - toggle nuclear screening (default: TRUE)"
	@echo "  F90=...          - Fortran compiler (default: gfortran)"
	@echo "  CXX=...          - C++ compiler (default: g++)"
	@echo ""
	@echo "Examples:"
	@echo "  make              # build with default settings"
	@echo "  make DEBUG=1      # build with debug flags"
	@echo "  make USE_SCREENING=FALSE  # build without screening"
	@echo "  make F90=ifort    # use Intel Fortran compiler"
