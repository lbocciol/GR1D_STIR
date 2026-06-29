include ../make.inc

SOURCES=GR1D_module.F90 \
	grad_module.F90 \
	timers.F90 \
	GR1D.F90 \
	atmosphere.F90 \
	eos.F90 \
	ppm.F90 \
	ye_of_rho.F90 \
	allocate_vars.F90 \
	initialize_vars.F90 \
	input_parser.F90 \
	problem.F90 \
	outinfo.F90 \
  output.F90 \
  output_control.F90 \
	grid.F90 \
	map_profile.F90 \
	mass_interior.F90 \
	shocktube.F90 \
	sedov.F90 \
	collapse.F90 \
	OSC.F90 \
	boundaries.F90 \
	driver.F90 \
	tools.F90 \
	Step.F90 \
	prim2con.F90 \
	reconstruct.F90 \
	flux_differences_HLLE.F90 \
	gravity.F90 \
	TVD_reconstruction.F90 \
	con2prim.F90 \
	press_sources.F90 \
	restart_H5.F90 \
	analysis.F90 \
	neutrino_pressure.F90 \
	start.F90 \
	turbulence.F90 \
	M1/M1_init.F90 \
	M1/M1_reconstruct.F90 \
	M1/M1_closure.F90 \
	M1/M1_explicitterms.F90 \
	M1/M1_implicitstep.F90 \
	M1/M1_conservativeupdate.F90 \
	M1/M1_control.F90 \
	M1/M1_updateeas.F90 \
	M1/M1_testcases.F90

MODULES=GR1D_module.mod
OBJECTS=$(SOURCES:.F90=.o )

EXTRAINCS += $(HDF5INCS)
EXTRAOBJECTS += $(HDF5LIBS)

ifeq ($(HAVE_NUC_EOS),1)
	DEFS += -DHAVE_NUC_EOS
	EXTRAINCS += -Inuc_eos
	EXTRADEPS += nuc_eos/nuc_eos.a
	EXTRAOBJECTS += nuc_eos/nuc_eos.a
endif
ifeq ($(HAVE_LEAK_ROS),1)
	DEFS += -DHAVE_LEAK_ROS
	EXTRAINCS += -Ileakage_rosswog
	EXTRAOBJECTS += leakage_rosswog/leakage_rosswog.a 
	EXTRADEPS += leakage_rosswog/leakage_rosswog.a 
endif	
ifeq ($(HAVE_HDF5_OUTPUT),1)
  DEFS += -DHAVE_HDF5_OUTPUT
  EXTRAINCS += -Ihdf5_output
  EXTRAOBJECTS += hdf5_output/hdf5_output.a
  EXTRADEPS += hdf5_output/hdf5_output.a
endif

EXTRAINCS += -Inulibtable
EXTRAOBJECTS += nulibtable/nulibtable.a
EXTRADEPS += nulibtable/nulibtable.a

ifeq ($(HAVE_BURN),1)
	DEFS += -DHAVE_BURN
	CXX ?= g++
	CXXFLAGS ?= -O2 -std=c++20 -DSCREENING
	EXTRAINCS += -Iburn -Iburn/pynucnet -Ihelmeos
	BURNOBJECTS = burn/pynucnet/fortran_interface.o burn/pynucnet/wrapper.o \
		burn/burn.o burn/composition.o burn/nse.o helmeos/helmholtzEOS.o
	# listed before LAPACK below so burn.o (needs dgesv) precedes -llapack at link
	EXTRAOBJECTS += $(BURNOBJECTS) -lstdc++
	EXTRADEPS += $(BURNOBJECTS)
endif

ifeq ($(HAVE_LAPACK),1)
	EXTRAOBJECTS += $(LAPACKLIBS)
	DEFS += -DHAVE_LAPACK
endif	

ifeq ($(HAVE_OMP),1)
  DEFS += -DHAVE_OMP
endif

GR1D:  GR1D_module.o $(EXTRADEPS) $(OBJECTS)  
	$(F90) $(LDFLAGS) -o ../GR1D $(OBJECTS) $(EXTRAOBJECTS)

$(OBJECTS): %.o: %.F90 
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

GR1D.o: GR1D_module.o

nuc_eos/nuc_eos.a: nuc_eos/*.F90 nuc_eos/*.f
	$(MAKE) -C nuc_eos

leakage_rosswog/leakage_rosswog.a: leakage_rosswog/leak_rosswog.F90
	$(MAKE) -C leakage_rosswog

nulibtable/nulibtable.a: nulibtable/*.F90
	$(MAKE) -C nulibtable

hdf5_output/hdf5_output.a: hdf5_output/*.F90 GR1D_module.o nulibtable/nulibtable.a
	$(MAKE) -C hdf5_output

ifeq ($(HAVE_BURN),1)
# pynucastro-generated interface (lowercase .f90 -> no cpp); produces pynet.mod
burn/pynucnet/fortran_interface.o: burn/pynucnet/fortran_interface.f90
	$(F90) $(F90FLAGS) -Jburn/pynucnet -Iburn/pynucnet -c $< -o $@

# pynucastro-generated C++ rate wrapper (do not edit)
burn/pynucnet/wrapper.o: burn/pynucnet/wrapper.cpp
	$(CXX) $(CXXFLAGS) -Iburn/pynucnet -c $< -o $@

# burn integrator; produces burn.mod, needs pynet.mod
burn/burn.o: burn/burn.F90 burn/pynucnet/fortran_interface.o
	$(F90) $(F90FLAGS) $(DEFS) -Iburn/pynucnet -Jburn -c $< -o $@

# GR1D-level composition (network species, plus inert free n,p); produces
# composition.mod, needs pynet.mod and burn.mod
burn/composition.o: burn/composition.F90 burn/burn.o burn/pynucnet/fortran_interface.o
	$(F90) $(F90FLAGS) $(DEFS) -Iburn -Iburn/pynucnet -Jburn -c $< -o $@

# NSE composition solver; produces nse.mod, needs composition.mod (and pynet.mod)
burn/nse.o: burn/nse.F90 burn/composition.o burn/pynucnet/fortran_interface.o
	$(F90) $(F90FLAGS) $(DEFS) -Iburn -Iburn/pynucnet -Jburn -c $< -o $@

# Helmholtz EOS module; produces wlhelmholtzeos.mod
helmeos/helmholtzEOS.o: helmeos/helmholtzEOS.F90
	$(F90) $(F90FLAGS) $(DEFS) -Ihelmeos -Jhelmeos -c $< -o $@

# main objects that `use burn` / `use composition` / `use wlHelmholtzEOS` need those .mod files first
start.o Step.o: burn/burn.o burn/nse.o
start.o eos.o: helmeos/helmholtzEOS.o
# eos.o uses the nuc_eos `eosmodule` (offset-table grids); needs its .mod built first
eos.o: nuc_eos/nuc_eos.a
start.o Step.o eos.o allocate_vars.o map_profile.o: burn/composition.o
endif

clean:
	rm -rf ../GR1D
	rm -rf *.o
	rm -rf *.mod
	rm -rf nuc_eos/*.o
	rm -rf nuc_eos/*.mod
	rm -rf nuc_eos/*.a
	rm -rf leakage_rosswog/*.a
	rm -rf leakage_rosswog/*.o
	rm -rf leakage_rosswog/*.mod
	rm -rf nulibtable/*.a
	rm -rf nulibtable/*.o
	rm -rf nulibtable/*.mod
	rm -rf M1/*.o
	rm -rf hdf5_output/*.a
	rm -rf hdf5_output/*.o
	rm -rf hdf5_output/*.mod
	rm -rf burn/*.o
	rm -rf burn/*.mod
	rm -rf burn/pynucnet/fortran_interface.o
	rm -rf burn/pynucnet/wrapper.o
	rm -rf burn/pynucnet/*.mod
	rm -rf helmeos/*.o
	rm -rf helmeos/*.mod
