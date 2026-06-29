# This follows in detail https://pynucastro.github.io/pynucastro/he-burning-example.html
import pynucastro as pyna

nuclei = ["p",
          "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
          "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
          "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
          "n13", "n14", "f18", "ne21", "na22", "na23"]

reaclib_lib = pyna.ReacLibLibrary()
core_lib = reaclib_lib.linking_nuclei(nuclei)

other_rates = [("c12(c12,n)mg23", "mg24"),
               ("o16(o16,n)s31", "s32"),
               ("o16(c12,n)si27", "si28")]

for r, mp in other_rates:
    _r = reaclib_lib.get_rate_by_name(r)
    new_rate = pyna.ModifiedRate(_r, new_products=[mp])    
    core_lib += pyna.Library(rates=[new_rate])

iron_peak = ["n", "p", "he4",
             "mn51",
             "fe52", "fe53", "fe54", "fe55", "fe56",
             "co55", "co56", "co57",
             "ni56", "ni57", "ni58"]

iron_reaclib = reaclib_lib.linking_nuclei(iron_peak)

weak_lib = pyna.TabularLibrary()
iron_weak_lib = weak_lib.linking_nuclei(iron_peak)

all_lib = core_lib + iron_reaclib + iron_weak_lib

rates_to_derive = all_lib.backward().get_rates()

# This is only supported for python
# # now for each of those derived rates, look to see if the pair exists
# for r in rates_to_derive:
#     fr = all_lib.get_rate_by_nuclei(r.products, r.reactants)
#     if fr:
#         all_lib.remove_rate(r)
#         d = pyna.DerivedRate(source_rate=fr, use_pf=True, use_unreliable_spins=True)
#         all_lib.add_rate(d)

all_lib.eliminate_duplicates()

net = pyna.PythonNetwork(libraries=[all_lib])
net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47"])
net.remove_nuclei(["cl35", "k39", "sc43", "v47"])

net.make_nn_g_approx(intermediate_nuclei=["fe53", "fe55", "ni57"])
net.remove_nuclei(["fe53", "fe55", "ni57"])

print(net.summary())
print(net.get_nuclei())

# Export to Fortran for GR1D
# Extract the finalized, approximated rates from approx_net
final_rates = net.get_rates()

# # Initialize the FortranNetwork using those extracted rates
fortran_net = pyna.FortranNetwork(rates=final_rates)
fortran_net.write_network('pynucnet/')

