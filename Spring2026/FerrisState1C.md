# A Computational Model of Mitochondrial One-Carbon Metabolism

*Xhesika Luarasi, Sky T. Pike, Judah De Paula — 2026*

## Abstract

A computational model of mitochondrial one-carbon metabolism, critical
in the production of the universal methylator S-adenosylmethionine, as
well as purines and thymidylate, has only a few times been
accomplished, and results have not agreed well with real-world
experiments. Additionally, new enzymes have been discovered and
further characterization of existing enzymes of mitochondrial
one-carbon metabolism has been accomplished since these computational
models were published, mostly using MatLab. This paper describes the
development of a Python-based computational model of mitochondrial
carbon metabolism that incorporates newly discovered enzymes and
kinetic parameters. Comparison with published bench-top experiment
results is performed and shows the potential of this model in
predicting one-carbon flux through the mitochondria, given the
concentration of one-carbon donors and the relative amount of each
mitochondrial one-carbon enzyme.

## 1. Introduction

- 1C metabolism
- Importance
- Compartmentalized
- Enzymes
- Reactions in cytoplasm
- Mitochondrial enzymes
- Reactions in mitochondria
- Many cytoplasmic models exist; few — maybe one, old — mitochondrial
  model
- Reasons to build this model:
  - Includes all mitochondrial 1C enzymes for the complete pathway at
    every stage
  - Few existing models incorporate the mitochondrial pathway
  - The pathway appears to direct 1C units from the methyl cycle to
    de novo purine synthesis
  - Contributes to understanding of aberrant 1C metabolism and
    dysregulation in diseases such as cancer

## 2. Methods

*(Python implementation — Judah De Paula)*

> **Revision note (2026-08-15):** this section was rewritten after a
> substantial architectural update to `Spring2026/` (commit `fae2aba`,
> "Same-state ODE architecture and mthfd2/2l equations encorporated").
> The rewrite replaced the previous fixed-step Euler integrator with a
> proper ODE solver, split the single MTHFD2 dehydrogenase step into
> four cofactor-specific fluxes across two independent enzymes, and
> replaced the fixed-flux GCS approximation with a full three-enzyme
> kinetic sub-model. It resolves the two numerical concerns raised in
> the previous draft of this section (timestep unit handling and
> per-reaction `dt` inconsistency — see §2.6). It also introduces two
> new bugs, flagged in §2.9 rather than fixed here. Both were fixed on
> 2026-10-09 in commit `d0e0ef7` (see §2.9).

### 2.1 Overview and code structure

The model is implemented in Python 3 as a small set of cooperating
modules under `FerrisState/Spring2026/`:

| Module | Role |
|---|---|
| `ModelParameters.py` | Compound name constants, the ordered `STATE_COMPOUNDS` tuple used by the solver, and all kinetic constants and the initial mitochondrial concentration pool (`MITO_POOL`) |
| `SerToFormateMmModel.py` | Reaction network — computes every reaction rate and every metabolite derivative from a single instantaneous state, and drives the ODE integration via `scipy.integrate.solve_ivp` |
| `GlycineCleavageSystem.py` | Glycine cleavage system (GCS) sub-model — three sequential enzyme rate laws (P-, T-, and L-protein), rates only, no state mutation |
| `ReportingTable.py` | Per-timepoint logging, CSV export, and per-compound time-course plotting (`matplotlib`) |
| `main.py` | Entry point — constructs the model from `ModelParameters`, runs it, and triggers reporting |

New dependencies introduced in this revision: `numpy` and
`scipy.integrate.solve_ivp`.

**Architectural change from the previous version:** the model
previously mutated a single `pool` dictionary in place, reaction by
reaction, within each manually-stepped loop iteration — so later
reactions in a step saw concentrations already updated by earlier ones
in the same step. The current version is a genuine same-state
architecture: `calculate_all_rates(state)` computes every reaction's
rate from one unmodified state, `calculate_derivatives(state, rates)`
turns those rates into `d[metabolite]/dt` for every tracked compound
without mutating anything, and `ode_system()` packages that as the
right-hand side of an ODE system that `scipy.integrate.solve_ivp`
integrates directly. No function anywhere in the reaction network
writes to a shared, partially-updated pool mid-step.

### 2.2 Kinetic rate laws

Five rate-law forms are used across the network:

**Single-substrate MM with substrate inhibition** (SHMT2 forward
direction only):

```
v = Vmax · [S] / (Km_app + [S] · (1 + [S] / Ki))
```

**Plain single-substrate MM** (SHMT2 reverse; also reused, with only
one substrate active, for Step 2b and Step 4, where water is treated
as constant and non-limiting):

```
v = Vmax · [S] / (Km + [S])
```

**Two-substrate MM** (irreversible general form; used for GCS
P-protein, and for Step 3):

```
v = Vmax · ([S1] / (Km1 + [S1])) · ([S2] / (Km2 + [S2]))
```

**Competitive-cofactor two-substrate MM** (MTHFD2/MTHFD2L, new this
revision — CH2THF saturation combined with two cofactors, NAD⁺ and
NADP⁺, competing for the same binding site):

```
v = Vmax · ([CH2THF] / (Km_CH2THF + [CH2THF]))
        · ([cofactor] / (Km_cofactor · (1 + [competing_cofactor] / Ki_competing_cofactor) + [cofactor]))
```

**GCS-specific multi-substrate forms** (new this revision, T- and
L-protein):

```
Ordered Bi-Bi (T-protein):
v = Vmax · A · B / (Ki_A·Km_B + Km_B·A + Km_A·B + A·B)      [Ki_A ≈ Km_A, stated approximation]

Ping-Pong Bi-Bi (L-protein):
v = Vmax · A · B / (Km_A·B + Km_B·A + A·B)
```

The previous version's `reversible_two_substrate_mm()` helper (a
single forward-minus-reverse rate shared across four different
reactions) has been removed. Only SHMT2 remains genuinely reversible
(forward and reverse rates computed and reported separately); every
other reaction in the network is now modeled as strictly
unidirectional.

### 2.3 Reaction network

Because every rate is now computed from one shared instantaneous
state rather than a sequence of in-place pool mutations, there is no
longer a meaningful "execution order" among reactions the way there
was in the previous Euler-stepped version — all of the following are
evaluated independently, from the same state, every time the solver
requests a derivative:

1. **SHMT2 (serine hydroxymethyltransferase 2), reversible:**
   `SER + THF ⇌ CH2THF + GLY`
2. **MTHFD2, NAD⁺-dependent branch:**
   `CH2THF + NAD+ → CH⁺THF + NADH`
3. **MTHFD2, NADP⁺-dependent branch:**
   `CH2THF + NADP+ → CH⁺THF + NADPH`
4. **MTHFD2L, NAD⁺-dependent branch:**
   `CH2THF + NAD+ → CH⁺THF + NADH`
5. **MTHFD2L, NADP⁺-dependent branch:**
   `CH2THF + NADP+ → CH⁺THF + NADPH`
   *(MTHFD2 and MTHFD2L are modeled as two independent enzyme
   populations that do not compete with or inhibit one another; within
   each enzyme, NAD⁺ and NADP⁺ compete for the same cofactor site —
   four separate fluxes in total, reported individually rather than
   collapsed into one, so NADH- vs. NADPH-producing flux and
   enzyme-of-origin stay distinguishable.)*
6. **MTHFD2/MTHFD2L cyclohydrolase (Step 2b), irreversible:**
   `CH⁺THF → CHOTHF + H⁺`
   (water constant/non-limiting, so the rate depends on `CH⁺THF`
   alone — previously modeled as reversible two-substrate MM including
   H2O and H⁺ explicitly; now single-substrate.)
7. **MTHFD1L (Step 3), irreversible:**
   `CHOTHF + ADP → Formate + THF + ATP`
   Modeled as substrate-level phosphorylation — ADP is consumed and
   ATP produced, the reverse pairing from the previous version (which
   consumed ATP and produced ADP). Inorganic phosphate (Pi) is not
   explicitly represented.
8. **ALDH1L2 (Step 4), irreversible:**
   `CHOTHF + H2O → Formate + THF`
   Reworked from the previous version's terminal step. Previously this
   step consumed formate and NADP⁺ to produce CO2 and NADPH
   (`Formate + NADP+ → CO2 + NADPH`) and was explicitly described as
   *not* representing ALDH1L2. It is now a direct hydrolysis of
   CHOTHF to formate + THF with no NADP⁺/CO2 involvement, and is
   explicitly labeled as ALDH1L2 hydrolase activity — a genuine change
   of assumption, not just a rename (see §2.7).
9. **Glycine cleavage system (GCS), three sequential enzyme steps,
   NAD⁺-dependent, replacing the previous fixed-flux approximation:**
   - P-protein: `GLY + H_OX → H_AM + CO2`
   - T-protein: `H_AM + THF → H_RED + CH2THF + NH3`
   - L-protein: `H_RED + NAD+ → H_OX + NADH + H+`
   - Net: `GLY + THF + NAD+ → CH2THF + CO2 + NH3 + NADH + H+`
   - `H_OX`, `H_AM`, and `H_RED` are the oxidized, aminomethylated, and
     reduced states of the H-protein's lipoyl-lysine swinging arm —
     new tracked state compounds, not present in the previous version.
     Note the cofactor also changed from NADP⁺ (previous fixed-flux
     version) to NAD⁺ (current kinetic version).

### 2.4 Kinetic parameters

| Reaction | Parameter | Value | Units |
|---|---|---|---|
| SHMT2 fwd | `Step1_Vmax_f` | 9756 | µM/h |
| SHMT2 fwd | `Step1_aKm_thf` | 20 (±8) | µM |
| SHMT2 fwd | `Step1_Ki_thf` | 99 (±17) | µM |
| SHMT2 fwd | `Step1_Ki_ser` | 602 (±370) | µM |
| SHMT2 rev | `Step1_Vmax_r` | 3608 | µM/h |
| SHMT2 rev | `Step1_aKm_gly` | 469.8 | µM |
| SHMT2 rev | `Step1_Km_ch2thf` | 980 | µM |
| Km_ser_app | `Step1_KM_SER_MAX` | 430.0 | µM |
| Km_ser_app | `Step1_ALPHA_KM_SER` | 158.0 | µM |
| Km_ser_app | `Step1_K05_THF` | 250.0 | µM |
| GCS P-protein | `GCS_P_Vmax` | 13200.0 | µM/h |
| GCS P-protein | `GCS_P_KM_GLY` | 6600.0 | µM |
| GCS P-protein | `GCS_P_KM_HOX` | 7.4 | µM |
| GCS T-protein | `GCS_T_Vmax` | 44600.0 | µM/h |
| GCS T-protein | `GCS_T_KM_HAM` | 584 | µM |
| GCS T-protein | `GCS_T_KM_THF` | 182100 | µM |
| GCS L-protein | `GCS_L_Vmax` | 13041 | µM/h |
| GCS L-protein | `GCS_L_KM_HRED` | 490 | µM |
| GCS L-protein | `GCS_L_KM_NAD` | 490 | µM |
| MTHFD2 + NAD⁺ | `MTHFD2_NAD_VMAX` | 477200.0 | µM/h |
| MTHFD2 + NAD⁺ | `MTHFD2_NAD_KM_CH2THF` | 354.0 | µM |
| MTHFD2 + NAD⁺ | `MTHFD2_NAD_KM_NAD` | 202.0 | µM |
| MTHFD2 + NAD⁺ | `MTHFD2_NAD_KI_NADP` | 352.0 (assumed, see §2.7) | µM |
| MTHFD2 + NADP⁺ | `MTHFD2_NADP_VMAX` | 1152000.0 | µM/h |
| MTHFD2 + NADP⁺ | `MTHFD2_NADP_KM_CH2THF` | 302.0 | µM |
| MTHFD2 + NADP⁺ | `MTHFD2_NADP_KM_NADP` | 352.0 | µM |
| MTHFD2 + NADP⁺ | `MTHFD2_NADP_KI_NAD` | 202.0 (assumed, see §2.7) | µM |
| MTHFD2L + NAD⁺ | `MTHFD2L_NAD_VMAX` | 1584000.0 | µM/h |
| MTHFD2L + NAD⁺ | `MTHFD2L_NAD_KM_CH2THF` | 130.0 | µM |
| MTHFD2L + NAD⁺ | `MTHFD2L_NAD_KM_NAD` | 147.0 | µM |
| MTHFD2L + NAD⁺ | `MTHFD2L_NAD_KI_NADP` | 537.0 (assumed, see §2.7) | µM |
| MTHFD2L + NADP⁺ | `MTHFD2L_NADP_VMAX` | 1296000.0 | µM/h |
| MTHFD2L + NADP⁺ | `MTHFD2L_NADP_KM_CH2THF` | 153.0 | µM |
| MTHFD2L + NADP⁺ | `MTHFD2L_NADP_KM_NADP` | 537.0 | µM |
| MTHFD2L + NADP⁺ | `MTHFD2L_NADP_KI_NAD` | 147.0 (assumed, see §2.7) | µM |
| Step 2b (cyclohydrolase) | `Step2b_Vmax_f` | 8000 | µM/h |
| Step 2b (cyclohydrolase) | `Step2b_Km_chpthf` | 60 | µM |
| Step 3 (MTHFD1L) | `Step3_Vmax_f` | 12000 | µM/h |
| Step 3 (MTHFD1L) | `Step3_Km_chothf` | 100 | µM |
| Step 3 (MTHFD1L) | `Step3_Km_adp` | 300 (in-code comment: "Temporary ?") | µM |
| Step 4 (ALDH1L2) | `Step4_Vmax` | 8640.0 (in-code comment: to be revised to 29.5 nmol/min/mg) | µM/h |
| Step 4 (ALDH1L2) | `Step4_Km_chothf` | 11.4 | µM |

**Two constants defined in `ModelParameters.py` are no longer
referenced by any rate function** in the current architecture:
`Step2b_Km_h2o` (40000 µM — vestigial from when Step 2b was modeled as
reversible two-substrate MM including water explicitly) and
`Step3_Km_atp` (300 µM — vestigial from when Step 3 used ATP as a
substrate rather than a product). Left in place, harmless, but worth a
cleanup pass.

Several constants (SHMT2 `Vmax_f`/`Vmax_r`, `aKm_gly`, `Km_ch2thf`) are
still annotated in-code as "Calculated. See word document" — that
worksheet is not currently part of the repository, so those
derivations aren't independently reproducible from the code alone.
Unchanged from the previous draft of this section.

### 2.5 Initial conditions

Simulation begins from a single mitochondrial matrix pool
(`MITO_POOL`), with concentrations in µM. Five new state compounds
were added this revision (GCS carrier states and the NAD⁺/NADH pool);
`CHOTHF`'s initial value changed from 1963 µM (previous version) to
50 µM, flagged in-code as `#doublecheck literature`:

| Compound | Initial concentration (µM) | Note |
|---|---|---|
| SER | 1963 | |
| THF | 20.8 | |
| NADP⁺ | 20 | |
| ATP | 300 | |
| H⁺ | 10000 | |
| H2O | 40000 | |
| CH2THF | 50 | |
| GLY | 1858 | |
| Formate | 0 | |
| ADP | 100 | |
| CHOTHF | 50 | changed from 1963 in the previous version; in-code: "doublecheck literature" |
| CH⁺THF | 1.407 | |
| NADPH | 200 | |
| CO2 | 0 | in-code: "Need a literature value" |
| NH3 | 0 | new; in-code: "Need a literature value" |
| H_OX | 10.0 | new; in-code: "Need a literature value" |
| H_AM | 0 | new |
| H_RED | 0 | new |
| NAD⁺ | 100.0 | new; in-code: "find ratio in literature" |
| NADH | 0 | new; in-code: "Find ration [sic] in literature" |

### 2.6 Simulation settings and numerical integration

The model now integrates the reaction network as a proper ODE system
using `scipy.integrate.solve_ivp` (explicit adaptive-step Runge-Kutta
4(5), `method='RK45'`, `rtol=1e-6`, `atol=1e-9`), replacing the
previous version's manually-stepped explicit Euler scheme entirely.
Configuration:

- `EXPERIMENT_RUNTIME_IN_HOURS = 1`
- `EXPERIMENT_TIMESTEPS = 1800`
- `INITIAL_SER_CONCENTRATION = 1963` µM

**`EXPERIMENT_TIMESTEPS` now means something different than before.**
It no longer controls the solver's internal step size — `solve_ivp`
chooses its own adaptive steps to meet the stated tolerances,
independent of this value. It now only determines how many evenly
spaced time points (`t_eval`, via `np.linspace`) the solution is
*reported* at (1801 points across the run, including t=0). This is a
reporting-resolution knob now, not an integration-accuracy knob.

**Both numerical concerns flagged in the previous draft of this
section are resolved:**

- The hours→seconds unit-conversion bug is fixed. `main.py` now
  computes `experiment_runtime_in_sec = EXPERIMENT_RUNTIME_IN_HOURS *
  3600` explicitly before constructing the model, and
  `ser_to_formate_mm_equation()` converts back to hours
  (`experiment_runtime_hr = experiment_runtime_in_sec / 3600.0`) to
  hand `solve_ivp` a `t_span` in the same units as the µM/h rate laws.
  Verified directly (see §2.9): a full simulated hour is what the
  solver actually receives and integrates over.
- The per-reaction `dt`-inconsistency (one reaction previously skipping
  the seconds→hours conversion its neighbors applied) is now
  structurally impossible rather than merely fixed — there is no
  per-reaction `dt` anywhere in the current architecture.
  `calculate_derivatives()` returns one shared `d[metabolite]/dt`
  vector for the whole network, and the solver integrates all of it
  together.

### 2.7 Assumptions

Assumptions carried over unchanged from the previous version:

- No CO2-consuming/venting reactions are modeled beyond the CO2
  produced by the GCS P-protein step.
- Serine and glycine transport rates across the mitochondrial membrane
  are not modeled; both are treated as available within the single
  well-mixed matrix pool.

Assumptions that changed this revision:

- **MTHFD2(L) and MTHFD1L are no longer modeled as reversible.**
  Previously reversible two-substrate MM with a stated "mainly forward
  in practice" justification; now genuinely unidirectional (no reverse
  rate term exists) for every step except SHMT2.
- **GCS is still forward-only, but is no longer a fixed literature
  flux.** It's now a concentration-dependent three-enzyme kinetic
  sub-model, and its cofactor changed from NADP⁺ to NAD⁺.
- **ALDH1L2 is now explicitly modeled** — reversing the previous
  version's stated assumption ("No dehydrogenase reaction for
  ALDH1L2"). It's modeled as hydrolase activity (`CHOTHF + H2O →
  Formate + THF`), not the NADP⁺-dependent dehydrogenase direction.

New assumptions introduced this revision:

- MTHFD2 and MTHFD2L are treated as independent, non-competing enzyme
  populations; NAD⁺ and NADP⁺ compete with each other only within a
  given enzyme's active site.
- Experimental mixed-cofactor inhibition constants (true `Ki` values
  for NAD⁺/NADP⁺ competition) are not currently available. As a stated
  baseline approximation, `Ki_competing_cofactor ≈ Km_competing_cofactor`
  for all four MTHFD2/MTHFD2L branches. `Km` and `Ki` are kept as
  separate parameters in code specifically so real `Ki` values can
  replace this approximation later without changing the rate-law
  structure.
- Inorganic phosphate (Pi) is not explicitly represented, despite Step
  3 now being modeled as substrate-level phosphorylation (ADP → ATP).
- Several new initial concentrations (CO2, NH3, H_OX, NAD⁺, NADH) are
  explicitly flagged in-code as placeholders pending a literature
  source, not yet sourced values.

### 2.8 Outputs

Each reported time point's full state, every individual reaction rate
(SHMT2 forward/reverse, all four MTHFD2/MTHFD2L branches, all three
GCS sub-rates, Step 2b/3/4), and both `time_hr` and `time_sec` are
logged to an in-memory table (`ReportingTable`), then:

- exported to `Spring2026_report.csv` (one row per reported time
  point), and
- rendered as one time-course plot per tracked quantity via
  `matplotlib`, saved to `plots/`.

The CSV and plots aren't committed to the repository (`plots/` and
`Spring2026_report.csv` are gitignored). They are generated by running
`main.py`. They were last regenerated on 2026-11-16 from the current
code (the `fae2aba` model plus the `d0e0ef7` fixes): 35 time-course
plots and an 1,801-row CSV covering the simulated hour. They reflect
the same-state ODE architecture, the new state compounds and the
per-branch rates, and can serve as the basis for the figures in §3.
Three older plots (`thf_to_ch2thf_*`, from the fixed-flux/Euler
version) may still sit in a local `plots/` folder. The current code no
longer produces them, since `fae2aba` renamed them `shmt2_*`, and they
shouldn't be used.

### 2.9 Known issues in `fae2aba` (resolved in `d0e0ef7`)

**Status (2026-11-16):** both bugs below were fixed on 2026-10-09 in
commit `d0e0ef7`, which restored a thin `run()` wrapper on
`SerToFormateMmModel` and updated `ReportingTable.py` to the
`Step1_`-prefixed constant names. `main.py` now runs to completion,
reproducing the values reported below (`Formate` ≈171.4 µM, `SER`
≈1774.1 µM at one hour). The original write-up is kept below for the
record.

**Open item:** two of Xhesika's branches have not been merged into
`main`. `my-edits` (2026-02-04) contains `Graphs.py`,
`OneCarbonMetabolism_Inhibition.py`, and versions of
`SMT2 Rate Bar Graph.py` and `THF conditioning.py` that differ from
those in `main`. `Xhesika-Luarasi-patch-1` dates from 2025-06-03. The
authors still need to decide whether either branch belongs in the
model.

Two bugs, found while reading the code for this section and confirmed
by direct execution. They weren't fixed at the time, since this is
Xhesika's commit and the fix was a call for the authors to make, not
something to patch silently in the middle of an active edit:

- **`main.py` cannot currently run to completion.** It calls
  `model.run()`, but `SerToFormateMmModel` no longer defines a `run()`
  method — it was removed as part of this revision, and nothing was
  added to replace it as the public entry point. Confirmed by running
  `main.py` directly: `AttributeError: 'SerToFormateMmModel' object
  has no attribute 'run'. Did you mean: 'ran'?`. The underlying method
  (`ser_to_formate_mm_equation`) is intact and callable directly; the
  simplest fix is restoring a thin `run()` wrapper (as in the previous
  version) or updating `main.py` to call
  `ser_to_formate_mm_equation()` directly.
- **`ReportingTable.pprint()` references stale constant names.** The
  Step 1 constants were renamed (`KM_SER_MAX` → `Step1_KM_SER_MAX`,
  and similarly for `ALPHA_KM_SER` and `K05_THF`) as part of this
  revision, but `ReportingTable.py`'s `pprint()` method — called
  immediately after `model.run()` in `main.py` — still references the
  old names. Confirmed directly: `NameError: name 'KM_SER_MAX' is not
  defined`. This would trigger regardless of the `run()` issue above,
  once that's fixed.

**With both bypassed directly** (calling
`ser_to_formate_mm_equation()` explicitly and skipping `report.pprint()`),
the ODE solver itself runs successfully to completion over the full
simulated hour and reports success (`solution.success == True`),
producing a defined final state — e.g. `Formate` rising from 0 to
≈171.4 µM and `SER` falling from 1963 to ≈1774.1 µM over the hour, in
one such run. Several fast-depleting pools (`CH⁺THF`, `CHOTHF`,
`H_OX`, `NAD⁺`, `NADP⁺`) settle to values that are negative only at
floating-point-noise scale (roughly 1e-11 to 1e-31 µM) rather than
meaningfully negative — consistent with near-total depletion, not
solver divergence, and each rate function already guards against
negative-concentration inputs (`if substrate <= 0: return 0.0`), which
keeps this from compounding. This confirms the ODE architecture itself
is sound; the two bugs above are entry-point/reporting issues, not
integration issues.

## 3. Results

- Primarily figures and tables of results looking at formate
  production and intermediate production
  - Varying concentrations of serine to glycine on CH2THF production
    and formate production
- Forward and reverse reactions until equilibrium
- Time to equilibrium

## 4. Conclusion / Discussion

- Partial mitochondrial and cytoplasmic 1C units
