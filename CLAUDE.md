# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

GR1D is a general-relativistic, spherically-symmetric, neutrino-radiation-hydrodynamics
code for stellar collapse and core-collapse supernovae, written in free-form Fortran
(`.F90`, preprocessed). See `docs/README.pdf` and the O'Connor & Ott 2010 / O'Connor 2015
papers referenced in `README.md` for the physics.

This checkout is on the `burn` branch: work in progress to add a **nuclear reaction-network
burning** stage and a **Helmholtz EOS** to the existing hydro+transport code (see
"Work in progress" below).

## Build

The build is `make.inc`-driven. There is no `make.inc` in the repo — create one from the
template and fill in the empty fields before building:

```sh
cp make.inc.template make.inc      # then edit: set F90, HDF5DIR, LAPACKDIR
make                               # top-level Makefile recurses into src/ -> ./GR1D
make clean
make DEBUG=1                       # -O0 -fcheck=all -fsanitize=address,undefined
```

`make.inc` selects the compiler (`F90=`), external library paths (`HDF5DIR`, `LAPACKDIR`),
and the **feature flags** that toggle `-D` defines and link which sub-library archives:

- `HAVE_OMP` — OpenMP (`-fopenmp`).
- `HAVE_NUC_EOS` — tabulated hot nuclear EOS (`src/nuc_eos/`, builds `nuc_eos.a`).
- `HAVE_LAPACK` — LAPACK/BLAS.
- `HAVE_HDF5_OUTPUT` — HDF5 output module (`src/hdf5_output/`).
- `HAVE_LEAK_ROS` — Rosswog leakage scheme (`src/leakage_rosswog/`).
- `HAVE_RESTART` — restart support.

Each sub-library (`nuc_eos`, `nulibtable`, `leakage_rosswog`, `hdf5_output`) builds its own
`.a` via its own Makefile and is linked into `../GR1D`. `nulibtable` (neutrino opacity tables)
is always linked.

## Run

```sh
./GR1D            # reads a file literally named "parameters" from the current directory
```

`./GR1D` does **not** take command-line arguments. Set up a run directory containing a
`parameters` file, a stellar `profile` (see `profiles/`), and the required EOS / opacity
HDF5 tables (download links in `README.md`). Copy and edit one of the ready-made parameter
sets in `sample_parameter_files/` (e.g. `latest_recommended_params`,
`oconnor2015_taggedversion_params/`, `turbulence_recommended_params`). `RunDir` and `GR1D`
are git-ignored, so runs are meant to happen in a scratch directory.

There is no automated test suite. Validation is done by running the sample parameter files
(several reproduce published test cases, e.g. the `TC*` cases in
`sample_parameter_files/oconnor2015_taggedversion_params/`).

## Architecture

**Entry point** `src/GR1D.F90`: calls `start` (setup), then an unbounded time-integration
loop: `SetTimeStep` → `handle_output` → `Step(dt)` → `postStep_analysis`, until `tend`/`ntmax`.

**Global state** lives in the `GR1D_module` module (`src/GR1D_module.F90`): every field
variable is a module-level `allocatable, save` array (`rho`, `temp`, `ye`, `eps`, `v`, etc.),
plus all run parameters. Nearly every routine does `use GR1D_module`. `src/timers.F90` holds
the `timer_*` accumulators printed at shutdown.

**Setup** `start` (`src/start.F90`) orchestrates: `input_parser` (reads the `parameters`
file via `get_*_parameter` helpers keyed by name), `grid`, `allocate_vars`,
`initialize_vars`, problem setup (`collapse`/`shocktube`/`sedov`/`OSC`/M1 test), and
`map_profile` (interpolates a 1D stellar profile onto the grid).

**`Step(dts)`** (`src/Step.F90`) is the operator-split heart of one timestep, in order:
1. **Hydro** — RK update (`iorder_hydro`) of finite-volume GR hydro: `reconstruct`
   (`tvd`/`ppm`/pc) → `prim2con`/`con2prim` → `flux_differences_HLLE` → `gravity`/metric →
   source terms. Timed in `timer_hydro`.
2. **Burn** — nuclear network update (work in progress, `HAVE_BURN`, `timer_burn`).
3. **Neutrinos** — `do_M1` two-moment transport (`src/M1/`, uses `nulibtable` opacities) or
   the leakage scheme.

**EOS** is dispatched through `eos(...)` in `src/eos.F90`, selected by `eoskey`
(1 = hybrid analytic, 2 = polytrope, 3 = tabulated hot nuclear `nuc_eos`, 4 = ideal gas) and
an `eosflag`/`keytemp` argument convention that says which thermodynamic variable is the
input vs. to be solved for. Read the header comments in `eos.F90` before touching call sites.

## Conventions

- Free-form Fortran preprocessed by the compiler (`.F90`). Feature code is guarded by
  `#ifdef HAVE_*` matching the `make.inc` flags — keep new optional features behind a flag and
  wire the flag through `src/Makefile` (`DEFS`, `EXTRA*`).
- `real*8` / `real(8)` throughout; `implicit none` in new code.
- Loops over the radial grid run `i = ghosts1+1, n1-ghosts1` (interior); `ghosts1` ghost
  zones on each side, total `n1 = radial_zones + 2*ghosts1`.

## Work in progress: burning network + Helmholtz EOS (`burn` branch)

Not yet wired into `make.inc`/`src/Makefile` — guarded by a new `HAVE_BURN` define. The main-code
integration hooks (in `Step.F90`, `start.F90`, `allocate_vars.F90`, `initialize_vars.F90`,
`map_profile.F90`, `GR1D_module.F90`, `timers.F90`) are **incomplete and currently contain
placeholder/pseudo-code** (e.g. a trailing-comma `use` statement in `Step.F90`, a prose
placeholder for the composition read in `map_profile.F90`); they will not compile under
`HAVE_BURN` until finished. Treat the diffs as a sketch, not working code.

- `src/burn/burn.F90` — module `burn`: a stiff backward-Euler + Newton integrator
  (`burn_state`, `burn_newton`, `burn_rhs`, `burn_init`) that wraps a **pynucastro**-generated
  reaction network. Works in molar abundances `Y = X/A`, CGS units, `dgesv` (LAPACK) for the
  Newton linear solve. `src/burn/CLAUDE.md` is the original task spec for this integrator.
- `src/burn/pynucnet/` — the **generated** network (Fortran module `pynet` + C++ rate headers).
  Do **not** hand-edit; regenerate via `src/burn/Generate_Network.py` (pynucastro), which writes
  this directory. Edit the isotope list there to change the network.
- `src/burn/GNUmakefile` — standalone build of `burn_test` (driver `burn_test.f90`) that links
  `burn.F90` against the generated network and LAPACK, independent of the main GR1D build. Use
  this to develop/test the integrator in isolation. `USE_SCREENING=TRUE` toggles `-DSCREENING`.
- `src/helmeos/helmholtzEOS.F90` — module `wlHelmholtzEOS`: a Helmholtz free-energy EOS
  (electron-positron + ions + radiation + Coulomb) with multi-mode root finding
  (`eos_input_*` selects which pair of `rho/T/p/e/h/s` are inputs). Intended for the burning
  regime; not yet hooked into `eos.F90`.
- New `GR1D_module` burn state: `nspec`, `aion`, `zion`, `Yion(:,:)` (per-zone composition),
  and thresholds `T_NSE` (assume NSE above) / `T_interp` (EOS interpolation).
