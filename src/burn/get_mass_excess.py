#!/usr/bin/env python3
"""
Extract mass excess values from pynucastro for the approx21 network.
"""

import pynucastro as pyna

# Load the ReacLib library
reaclib_lib = pyna.ReacLibLibrary()

# approx21 isotopes
isotopes = ["n", "h1", "he4", "c12", "n14", "o16", "ne20", "mg24",
            "si28", "s32", "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
            "fe54", "cr56", "fe56"]

print("Extracting mass excess values from pynucastro for approx21 network:")
print(f"Number of isotopes: {len(isotopes)}")
print()

# Extract mass excess for each isotope
mass_excess = []
for iso_name in isotopes:
    # Get the nucleus from pynucastro
    nucleus = pyna.Nucleus(iso_name)
    # Extract mass excess in MeV
    # pynucastro stores this in nucleus.Q_value or similar
    # Let's print available attributes
    print(f"{iso_name:5s}: {nucleus}")

print()
print("To get the actual mass excess values, we need to use the nuclear data library.")
print("The values are typically retrieved via pynucastro's interface to JINA ReacLib.")
