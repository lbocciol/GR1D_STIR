# GR1D + burning network — repository guide

GR1D-style 1D (spherical) GR/Newtonian supernova hydro code. A pynucastro
burning network is being added and coupled in operator-split fashion. The active
work is the **composite equation of state** that hands off between the finite-T
`nuc_eos` (NSE) at high T and a Helmholtz EOS at low T.

## Guiding principles (in priority order)

1. **Thermodynamic consistency first.** Pressure, energy, entropy and their
   derivatives (`dpdrhoe`, `dpderho`, `cs2`, `dedt`) must come from a single
   consistent free energy in each regime, and must vary *continuously* across the
   transition. Never blend quantities in a way that breaks `cs2 > 0` or makes
   `de/dT < 0`. When interpolating between the two EOS, interpolate the **same
   set of returned quantities** with the **same weight**, so the result stays a
   valid thermodynamic state.
2. **Clarity / readability over cleverness.** The EOS is cheap — do not optimize
   it. Prefer an obvious, well-commented blend over a fast opaque one. Someone
   reading `eos.F90` should see *exactly* which branch ran and why.
3. **Minimal, gated changes.** With `HAVE_BURN = 0` the code must be identical to
   today. All new EOS-blending logic lives behind `#ifdef HAVE_BURN` (and the
   existing `#if HAVE_NUC_EOS`). Do not perturb the non-burning code paths.

## Preprocessor flags

- `HAVE_BURN`   — enables network + composite EOS. **Default-off must be a no-op.**
- `HAVE_NUC_EOS`— enables the finite-T table EOS (`eoskey == 3`).
Files are `.F90` (capital F → cpp runs). Guard new code; never assume a flag is set.

## EOS architecture (where to work)

`eos.F90` is the single entry point for all thermodynamics. Two interfaces:

- `eos_full(i, rho, temp, ye, eps, press, pressth, ent, cs2, dedt, dpderho,
  dpdrhoe, xa,xh,xn,xp, abar,zbar, mu_e,mu_n,mu_p,munu, keytemp,keyerr,eoskey,
  rfeps)` — computes the **whole** thermo vector at once. **This is what
  `Step.F90` calls.** Do the composite logic here.
- `eos(i, ri, tio, y, eio, xx, keytemp,keyerr, eosflag, eoskey, rfeps)` — returns
  **one** quantity selected by `eosflag` (1=p, 2=dpdrhoe, 3=dpderho, 4=eps,
  5=temp, 6=cs2, 7=gamma, 8=entropy, 9=munu). Mirror any composite logic here too.

`eoskey` selects the backend: 1=hybrid, 2=poly, 3=finite-T `nuc_eos`, 4=gamma-law.
The composite (nuc_eos ↔ Helmholtz) belongs to the **`eoskey == 3`** branch, only
when `HAVE_BURN` is set.

### Units convention (do not break this)

Every sub-EOS works in **CGS**, then the caller scales to geometric units with
`rho_gf, press_gf, eps_gf` (from `GR1D_module`). The Helmholtz path must follow
the *identical* convert-in/convert-out pattern as the existing `nuc_eos` branch,
so the blend happens in CGS and the `*_gf` scaling is applied once, uniformly.

### `keytemp` (root-find mode), preserve semantics

- `keytemp == 0`: T unknown, solve from `eps` (and ρ, Ye, composition).
- `keytemp == 1`: T known, return everything (energy set from T).
Both backends must honor the same `keytemp` so the composite is transparent to
callers. `keyerr /= 0` signals EOS failure (see the bounce-recovery loop in
`Step.F90`); propagate it, don't swallow it.

## The two backends

**nuc_eos** (`nuc_eos/`): called via `nuc_eos_short(rho,temp,ye,eps,prs,ent,cs2,
dedt,dpderho,dpdrhoe,munu,keytemp,keyerr,rfeps)`, CGS in/out. Valid for T ≥ T_NSE.

**Helmholtz** (`helmeos/helmholtzEOS.F90`, `MODULE wlHelmholtzEOS`):
`FullHelmEOS(input, HelmTable, HelmholtzState)`.
- `input`: `eos_input_rt` (ρ,T → everything) or `eos_input_re` (ρ,e → solve T).
- Needs `abar, zbar, ye` set on `HelmholtzState` (from the network composition).
- Returns CGS; entropy already in `k_B`/baryon; provides `dpdr_e`, `dpde`, `cs`,
  `gam1`, `cv` (= `dedT`) — the consistent set needed to fill `eos_full` outputs.
- **A table loader is NOT yet present.** A `helm_table.dat` must be read once into
  a saved `HelmTableType` at startup (mirror how `nuc_eos` reads its table; call
  it from the same place under `#ifdef HAVE_BURN`). Flag this as a prerequisite.

## Transition / blending

`T_NSE = 5.8e9`, `T_interp = 5.0e9` already live in `GR1D_module.F90` under
`#ifdef HAVE_BURN`. Rules:
- T ≥ T_NSE → nuc_eos only. T ≤ T_interp → Helmholtz only.
- In between → call **both** and blend with a single smooth weight
  `w(T) ∈ [0,1]` (e.g. a smoothstep in `T` between `T_interp` and `T_NSE`).
  Blend **all** returned thermodynamic quantities with the **same** `w` so the
  state stays consistent; document the chosen weight function.
- When T is the unknown (`keytemp == 0`), decide the regime from the current/temp
  estimate carefully so the blend is single-valued and convergent — document the
  approach. Keep it simple; the EOS is cheap, so an extra evaluation is fine.

## Composition / state

Under `#ifdef HAVE_BURN`, `GR1D_module.F90` declares `nspec, aion(:), zion(:),
Yion(:,nzones)` and the thresholds above. `abar`/`zbar` for the Helmholtz call
are derived from `Yion` (`abar = 1/Σ Y_i`, `zbar = abar·Σ Z_i Y_i`). Reading the
initial composition into `Yion` is a *later* task — stub/placeholder is fine for
the EOS interface work.

## Burning coupling (already in place — for reference)

`Step.F90` (under `#ifdef HAVE_BURN`) calls `burn_newton(rho,temp,Yion(:,i),dts,
e_step,...)` per zone after the hydro update, then `eps += e_step`, **only where
T < T_NSE**. See `burn/CLAUDE.md` for the network integrator's own rules. Do not
edit the pynucastro-generated interface.

## Workflow expectations

- Work the **EOS interface first** (`eos.F90` + Helmholtz table loader), then
  composition I/O, then anything else.
- After any change, confirm `HAVE_BURN = 0` still produces the original code path
  (diff the non-guarded regions; build both with and without the flag).
- All reals `real*8`/`REAL(8)`, `implicit none`, match surrounding free-form style.
