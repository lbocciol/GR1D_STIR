include ../../make.inc

# List of source files
FFC_SOURCES = FFC.F90 \
              tools_FFC.F90

# List of object files (compiled .o files)
FFC_OBJECTS = $(FFC_SOURCES:.F90=.o) ../GR1D_module.o

# Extra include directories and object files if needed
EXTRAINCS = "-I../nulibtable" "-I../"
EXTRAOBJECTS = 

# All target that will be built
all: FFC.a

# Compile each Fortran source file into an object file
$(FFC_OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

# Ensure FFC.a is built but does not contain nulibtable.a
FFC.a: $(FFC_OBJECTS)
	ar rcs FFC.a $(FFC_OBJECTS)

# Clean up object files, modules, and static libraries
clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
