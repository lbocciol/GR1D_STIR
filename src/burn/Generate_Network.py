import pynucastro as pyna

# Load the ReacLib library
reaclib_lib = pyna.ReacLibLibrary()

# Combine them into a single master library
full_lib = reaclib_lib

# approx 21
isotopes = ["n", "h1", "he4", "c12", "n14", "o16", "ne20", "mg24",
            "si28", "s32", "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
            "fe54", "cr56", "fe56"]

# approx 13
isotopes = ["he4", "c12" , "o16" , "ne20", "mg24", "si28", \
            "s32", "ar36", "ca40", "ti44", "cr48", "fe52", "ni56"]

# Create the linked network using the combined library
linked_lib = full_lib.linking_nuclei(isotopes)

# Export to Fortran for GR1D
net = pyna.FortranNetwork(libraries=[linked_lib])
net.write_network('src/burn/pynucnet')