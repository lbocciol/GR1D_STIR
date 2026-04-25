include ../../make.inc

SOURCES=eosmodule.F90 readtable.F90 nuc_eos.F90 bisection.F90 findtemp.F90 linterp_many.F90
FSOURCES=linterp.f

CLEANSTUFF=rm -rf *.o *.mod *.a

OBJECTS=$(SOURCES:.F90=.o)
FOBJECTS=$(FSOURCES:.f=.o)

EXTRADEPS=timers.o

MODINC=$(HDF5INCS)

all: nuc_eos.a

nuc_eos.a: $(OBJECTS) $(FOBJECTS) timers.o
	ar r nuc_eos.a *.o

$(OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) -c $< -o $@

$(FOBJECTS): %.o: %.f $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) -c $< -o $@

timers.o: ../timers.F90
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) -c $< -o $@

clean:
	$(CLEANSTUFF)
