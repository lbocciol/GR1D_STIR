include ../../make.inc

H5_SOURCES=hdf5_output_utils.F90 \
           output_HDF5.F90

H5_OBJECTS=$(H5_SOURCES:.F90=.o)

EXTRAINCS = $(HDF5INCS) -I.. -I../nulibtable
# Note: You don't usually need EXTRAOBJECTS (like HDF5LIBS) when creating the .a archive,
# those are resolved later when linking the final executable in the main directory.

all: hdf5_output.a

# --- IMPORTANT: Internal Dependencies ---
# Even in a separate directory, you still must define the internal compilation order
# so the .mod files are generated in the right sequence.
output_HDF5.o: hdf5_output_utils.o

$(H5_OBJECTS): %.o: %.F90
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

hdf5_output.a: $(H5_OBJECTS)
	ar r hdf5_output.a $(H5_OBJECTS)

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
