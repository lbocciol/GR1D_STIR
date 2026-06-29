import pynucastro as pyna
from pynucastro.rates.aprox_family_rates import make_CO_approx_rates

net = pyna.network_helper(["p", "he4",
                           "c12", "o16", "ne20", "na23",
                           "mg24", "al27", "si28", "p31", "s32",
                           "cl35", "ar36", "k39", "ca40",
                           "sc43", "ti44", "v47", "cr48",
                           "mn51", "fe52", "co55", "ni56"])

approx_net = pyna.PythonNetwork(rates=net.get_rates())

approx_net.make_CO_burning_approx("C")
approx_net.remove_nuclei(["na23"])

approx_net.make_CO_burning_approx("CO")
approx_net.remove_nuclei(["al27"])

approx_net.make_CO_burning_approx("O")
approx_net.remove_nuclei(["p31"])

approx_net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47", "mn51", "co55"])
approx_net.remove_nuclei(["cl35", "k39", "sc43", "v47", "mn51", "co55"])

print(approx_net.summary())

# Export to Fortran for GR1D
# Extract the finalized, approximated rates from approx_net
final_rates = approx_net.get_rates()

# # Initialize the FortranNetwork using those extracted rates
fortran_net = pyna.FortranNetwork(rates=final_rates)
fortran_net.write_network('pynucnet/')
