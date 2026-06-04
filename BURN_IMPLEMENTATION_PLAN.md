# Burning + Helmholtz EOS — Implementation Plan (plain-English)

This file explains, in easy steps, what it takes to finish wiring the **nuclear
burning network** and the **Helmholtz EOS** into GR1D on the `burn` branch.
It is a roadmap, not code. Read it top to bottom.

The two governing specs are already written:
- `CLAUDE.md` (root) — the big picture of GR1D.
- `src/CLAUDE.md` — the detailed rules for the composite EOS (the hard part).
- `src/burn/CLAUDE.md` — the rules for the network integrator itself.

---

## 1. What we are actually building

Three new physics pieces, glued onto the existing hydro+neutrino code:

1. **A reaction network** (`src/burn/`) — pynucastro generated a 13-isotope alpha
   network (`he4 … ni56`). `burn.F90` integrates it per zone over a hydro substep
   and returns the energy released. **This code is essentially written.**
2. **A Helmholtz EOS** (`src/helmeos/helmholtzEOS.F90`) — an analytic
   electron/ion/radiation EOS for the *low-temperature, non-NSE* regime where the
   network governs composition. **The math routine is written; the table loader is not.**
3. **A composite EOS** — the glue that makes `eoskey==3` use the tabulated nuclear
   EOS (`nuc_eos`) when hot (T ≥ `T_NSE`) and Helmholtz when cool (T ≤ `T_interp`),
   smoothly blending in between. **This is the real work and does not exist yet.**

Everything is gated behind a new `HAVE_BURN` preprocessor flag. **Rule #1: with
`HAVE_BURN=0` the binary must be byte-for-byte the old GR1D.**

---

## 2. Current state — what works and what is broken

The diffs on this branch are a **sketch**. They do not compile under `HAVE_BURN`.
Here is the honest inventory.

### Already working (leave alone)
- `src/burn/pynucnet/` — generated network (module `pynet`). Exposes `nspec=13`,
  `aion`, `zion`, `rhs_f`, `jac_f`, `ener_gener_f`, `network_init`. **Never hand-edit.**
- `src/burn/burn.F90` — the integrator (`burn_init`, `burn_rhs`, `burn_state`,
  `burn_newton`). Backward-Euler + Newton, LAPACK `dgesv`.
- `src/burn/burn_test.f90` + `GNUmakefile` — standalone test, builds `burn_test`.
- `src/helmeos/helmholtzEOS.F90` — `FullHelmEOS(input, HelmTable, State)` works
  *if* you hand it a filled `HelmTableType`.

### Broken / missing (the to-do list)
There are **5 compile-blocking bugs** in the current sketch and **6 missing pieces**:

| # | Problem | File |
|---|---------|------|
| B1 | `use burn, only: burn_newton,` — trailing comma | `src/Step.F90:12` |
| B2 | `e_step` is used but never declared | `src/Step.F90` |
| B3 | `burn_state` calls `be_newton` — that routine doesn't exist (it's `burn_newton`) | `src/burn/burn.F90:67` |
| B4 | `read(67,*) read composition for each zone…` — English prose, not Fortran | `src/map_profile.F90` |
| B5 | `lprofile_comp_name` is used but never declared | `src/map_profile.F90` |
| M1 | `nspec` is declared in `GR1D_module` but **never set** (arrays size to 0) | setup |
| M2 | `Yion` shape is inconsistent: allocated `(n1,nspec)`, used `(:,i)` and `(i,j)` | several |
| M3 | No `make.inc`/`Makefile` wiring for `HAVE_BURN`, helmeos, or the C++ objects | build |
| M4 | Helmholtz EOS has **no table loader** and there is no `helm_table.dat` | `src/helmeos/` |
| M5 | The composite blend in `eos.F90` (`eoskey==3`) **does not exist** | `src/eos.F90` |
| M6 | Composition input file format + reader is undefined | `src/map_profile.F90` |

---

## 3. The plan, in four phases

Do them **in this order**. Each phase ends at a state you can compile and check.

### Phase 0 — Make it build (no new physics) ≈ half a day
Goal: `make HAVE_BURN=1` compiles and links; `HAVE_BURN=0` is unchanged.

1. **Build wiring (M3).**
   - Add `HAVE_BURN=0` to `make.inc.template` (and your `make.inc`).
   - In `src/Makefile`: under `ifeq ($(HAVE_BURN),1)` add `DEFS += -DHAVE_BURN`,
     compile `helmeos/helmholtzEOS.F90` and the burn objects
     (`burn/burn.F90` + `burn/pynucnet/fortran_interface.f90` +
     `burn/pynucnet/wrapper.o` from g++), and link `-lstdc++` plus LAPACK
     (LAPACK is already available via `HAVE_LAPACK`). Mirror how `nuc_eos.a` is built.
2. **Fix the 5 compile bugs (B1–B5):**
   - B1: `use burn, only: burn_newton` (drop comma) — or import `burn_state` instead
     (see decision D1 below).
   - B2: declare `real*8 :: e_step` in `Step.F90`.
   - B3: rename the call in `burn.F90` from `be_newton` to `burn_newton`.
   - B4/B5: handled in Phase 2 — for now stub the composition read so it compiles.
3. **Set `nspec` and fix `Yion` shape (M1, M2).** Decide the array layout once
   (recommend `Yion(nspec, n1)` so `Yion(:,i)` is one zone's vector — matches
   `Step.F90`). Set `nspec` from `pynet`'s `nspec` in `start.F90` *before*
   `allocate_vars`, and fix the `allocate(Yion(...))` and every index accordingly.

**Checkpoint:** builds both ways; a collapse run with `HAVE_BURN=1` but burning
effectively idle behaves like before.

### Phase 1 — The composite EOS (the hard part) ≈ 2–3 days
This is where `src/CLAUDE.md` is your bible. Thermodynamic consistency first.

4. **Helmholtz table loader (M4).** Write one new routine (e.g. `ReadHelmTable`)
   in `helmholtzEOS.F90` that reads `helm_table.dat` once into a `save`d
   `HelmTableType`. Call it from `start.F90` under `#ifdef HAVE_BURN`, next to the
   other table reads. **Prerequisite:** obtain `helm_table.dat` (the standard
   Timmes Helmholtz table, 271×101 grid — same file MESA/flash use).
5. **Composite branch in `eos.F90` (M5).** In the `eoskey==3` block of **both**
   `eos_full` (used by `Step.F90`) and `eos` (single-quantity), add, under
   `#ifdef HAVE_BURN`:
   - T ≥ `T_NSE` → call `nuc_eos_short` only (today's behaviour).
   - T ≤ `T_interp` → call `FullHelmEOS` only (needs `abar,zbar,ye` from `Yion`:
     `abar = 1/ΣY_i`, `zbar = abar·ΣZ_iY_i`).
   - In between → call **both** and blend **every** returned quantity with the
     **same** smooth weight `w(T)` (a smoothstep between `T_interp` and `T_NSE`).
   - Keep the **CGS-in / `*_gf`-scale-out** pattern identical to the existing branch.
   - Honour `keytemp` (0 = solve T from `eps`, 1 = T known) and propagate `keyerr`.

**Checkpoint:** at fixed composition, pressure/`cs2`/`dedt` are continuous across
`T_interp…T_NSE`; `cs2>0` and `dedt>0` everywhere.

### Phase 2 — Composition I/O + coupling ≈ 1 day
6. **Composition input (M6, B4, B5).** Declare `lprofile_comp_name` in
   `GR1D_module`, add it to `input_parser.F90`, define the file format (e.g. one
   row per profile zone, `nspec` mass-fraction columns), write the real read loop,
   convert X→Y into `Yion`, and map onto the grid with the existing `map_map`.
7. **Burn call in `Step.F90` (D1).** After the hydro update, loop interior zones,
   and **only where `T < T_NSE`** call the integrator, then `eps += e_step`.
   Decide D1 (below) — `burn_state` (does X↔Y + Ye check) vs `burn_newton` (raw,
   on molar `Yion`). Recommend `burn_state` for safety; keep `Yion` in molar `Y`.

**Checkpoint:** a test run burns He/C/O at low T, releases energy, `Ye` stays put.

### Phase 3 — Validation & polish ≈ 1 day
8. Run `src/burn/burn_test` (already builds) to sanity-check the network alone.
9. Run a small collapse/heating problem with `HAVE_BURN=1`; confirm energy
   conservation and that the EOS handoff doesn't kick `cs2` negative.
10. Diff all non-`#ifdef HAVE_BURN` regions against `master` to prove the
    default-off no-op rule (src/CLAUDE.md principle #3).

---

## 4. Decisions you need to make (small but real)

- **D1 — Step coupling entry point.** Call `burn_state` (handles X↔Y, floors
  negatives, Ye assertion, returns integrated energy) or `burn_newton` (raw molar
  step)? *Recommendation: `burn_state`.* This changes the `use` line in `Step.F90`
  and the per-zone call.
- **D2 — `Yion` units.** Store molar abundances `Y` (matches `burn_newton`/network)
  or mass fractions `X`? *Recommendation: molar `Y`*, convert to `X`/`abar`/`zbar`
  only at the EOS and output boundaries.
- **D3 — Blend weight `w(T)`.** Pick the smoothstep form and document it once.
- **D4 — `helm_table.dat` provenance.** Where the table comes from and where it
  lives at runtime (alongside the EOS/opacity tables in the run dir).

---

## 5. How big is this, really?

**Scope: moderate.** No architectural change to GR1D — burning is one more
operator-split stage and the composite EOS lives entirely inside the existing
`eoskey==3` branch. The bulk of the new physics code (network integrator +
Helmholtz math) is already written; what remains is **glue + one table loader +
one blend function**.

| Area | Files touched | New code | Risk |
|------|--------------|----------|------|
| Build wiring | `make.inc.template`, `src/Makefile` | small | low |
| Compile-bug fixes | `Step.F90`, `burn.F90`, `map_profile.F90`, `GR1D_module.F90` | tiny | low |
| Setup (nspec, alloc, init) | `start.F90`, `allocate_vars.F90`, `initialize_vars.F90` | small | low |
| **Composite EOS** | `eos.F90`, `helmholtzEOS.F90` (+table loader) | **medium** | **high** |
| Composition I/O | `input_parser.F90`, `map_profile.F90` (+format) | small–medium | medium |
| Data | `helm_table.dat`, a composition profile | (obtain) | — |

**Totals:** ~11 source files to edit (8 already partly sketched), **2 genuinely new
routines** (Helmholtz table loader + the EOS blend), **1 data file to obtain**, and
**5 one-line bug fixes**. Roughly **4–6 focused days**, with ~70% of the effort and
~90% of the *risk* concentrated in Phase 1 (the composite EOS), because that is the
only part that must be thermodynamically consistent rather than merely correct.

The smallest path to "it runs": **Phase 0** alone (≈ half a day) gets a compiling,
default-off-safe binary with burning idle — a good first commit before touching the EOS.
