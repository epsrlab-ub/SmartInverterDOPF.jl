# Three-phase Volt-VAr droop in a distribution OPF

The three-phase counterpart of the single-phase 33-bus case: **LinDist3Flow** as the
network model, on a real unbalanced LV feeder, with **all three droop encodings**.

Run from this directory:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"

julia --project=. LinDist3Flow_Lambda.jl        # MILP, Gurobi
julia --project=. LinDist3Flow_BigM.jl          # MILP, Gurobi
julia --project=. LinDist3Flow_Heaviside.jl     # NLP,  Ipopt
julia --project=. plot_network.jl               # feeder schematic
julia --project=. generate_results.jl           # refresh the documentation's data
```

`generate_results.jl` runs all three scripts and writes
`docs/src/assets/results/threephase/`, which is what the
[three-phase section of the tutorial](https://ra-emami.github.io/SmartInverterDOPF.jl/dev/tutorial_voltvar/#Three-phases)
is built from.

The three scripts share their skeleton **verbatim** — data, network model, inverter
constraints, objective, results and figures are byte-identical. They differ only in the
fenced `DROOP BLOCK`, and Heaviside additionally swaps Gurobi for Ipopt. Diff any two to
see exactly what an encoding costs.

Gurobi is [free for academic use](https://www.gurobi.com/academia/academic-program-and-licenses/);
Ipopt is open source, so the Heaviside script needs no licence at all.

## The three encodings agree

Same feeder, same inverters, same objective — only the droop encoding changes:

| encoding | class | solver | variables | binaries | constraints | solve | curtailed | max \|q − q_curve\| |
|:--|:--|:--|--:|--:|--:|--:|--:|--:|
| Lambda / SOS2 | MILP | Gurobi | 183,744 | 5,760 | 218,304 | 9.3 s | 46.32 kWh | 8.4e-07 |
| Big-M | MILP | Gurobi | 179,136 | 5,760 | 225,216 | 7.5 s | 46.32 kWh | 6.1e-15 |
| Heaviside | NLP | Ipopt | 171,072 | **0** | 207,936 | 30.5 s | 46.32 kWh | 1.1e-15 |

Identical curtailment to the digit, identical voltage range (0.9923 – 1.0207 p.u.), and
every dispatch point on the curve of its own class. Heaviside carries no integers at all
and pays for it in solve time; it converged cleanly to `LOCALLY_SOLVED` here, helped by
the linear host giving it no outer iteration to fight.

## What it does

| | |
|:--|:--|
| feeder | `network_5_Feeder_2` — real ENWL LV feeder, Kron-reduced to three wires |
| size | 194 buses, 193 lines, 489 m, 415/240 V |
| loads | 18 single-phase, split **4 / 5 / 9** across phases (3.7 / 6.0 / 7.7 kW) — genuinely unbalanced |
| PV | **12 smart inverters in 4 size classes**, 84 kW total, at the electrically farthest load buses, 4 per phase |
| horizon | 24 h at 15-minute resolution (96 steps), residential load shape, clear-sky irradiance |
| host | LinDist3Flow — phase-coupled linear branch flow |
| droop | Lambda / SOS2, Big-M, or Heaviside on the IEEE 1547 curve |
| objective | minimise total PV curtailment |

The point of the feeder choice is the unbalance. Loads are single-phase and unevenly
distributed, so the three phases genuinely diverge — which is the only reason to model
three phases at all.

## The twelve smart inverters

Four classes. Since `q̄ = S_max`, each class follows a **different droop curve** — same
breakpoint voltages, four saturation levels — which is what makes the dispatch figure
worth looking at.

| class | P rated | S max | q̄ | q̄ (p.u.) | sites |
|:--|--:|--:|--:|--:|--:|
| A | 3 kW | 3.30 kVA | 3.30 kvar | 0.0330 | 3 |
| B | 5 kW | 5.50 kVA | 5.50 kvar | 0.0550 | 3 |
| C | 8 kW | 8.80 kVA | 8.80 kvar | 0.0880 | 3 |
| D | 12 kW | 13.20 kVA | 13.20 kvar | 0.1320 | 3 |

Each phase gets one inverter of each class, and the class order is rotated by phase, so
size is confounded with neither phase nor distance from the substation:

| phase | bus 1 | bus 2 | bus 3 | bus 4 |
|:-:|:--|:--|:--|:--|
| 1 | 184 (A) | 74 (B) | 73 (C) | 45 (D) |
| 2 | 188 (B) | 157 (C) | 153 (D) | 145 (A) |
| 3 | 193 (C) | 179 (D) | 142 (A) | 149 (B) |

## Results

```
PV energy available : 476.9 kWh
PV energy delivered : 430.6 kWh
PV curtailment      :  46.32 kWh  (9.712 %)
voltage range       : 0.9923 – 1.0207 p.u.
  phase 1           : 0.9989 – 1.0127
  phase 2           : 0.9923 – 1.0207
  phase 3           : 0.9947 – 1.0150
```

Two checks run automatically at the end:

- **`max |q_dispatch − q_curve| = 8.4e-07 p.u.`** — every one of the 12 × 96 optimised
  operating points lies on the droop curve **of its own class**. The Lambda encoding is
  exact, and it does not care that the host is three-phase or that the fleet is mixed.
- **`LinDist3Flow vs exact AC: max |Δv| = 1.2e-03 p.u.`** — the linear host is checked
  against a full backward/forward sweep on the same dispatch, at the busiest PV step.

Each script writes three figures, suffixed by encoding — for example
`droop_dispatch_3ph.png` (Lambda), `droop_dispatch_3ph_bigm.png`,
`droop_dispatch_3ph_heaviside.png` — plus `network_schematic.png` from `plot_network.jl`.

## The formulation

For a line with 3×3 phase impedance `Z`, the phase-coupled voltage drop is

```
w_j = w_i − (aR·P_ij + aX·Q_ij),     aR[φ,ψ] =  2·Re(α^{ψ−φ}·Z[φ,ψ])
                                      aX[φ,ψ] =  2·Im(α^{ψ−φ}·Z[φ,ψ]),   α = e^{−j2π/3}
```

derived from `|V_j|² = |V_i − Z·I|²` with the quadratic term dropped and voltages assumed
near-balanced. The ±√3 cross terms in the published form of these matrices are just that
rotation written out. The code works in magnitude rather than squared magnitude
(`w_j − w_i ≈ 2·(v_j − v_i)` near nominal), matching the single-phase version, so the
droop breakpoints stay in ordinary p.u. voltage.

Two sanity checks are worth keeping in mind: for a single phase `α⁰ = 1` gives `aR = 2r`,
`aX = 2x`, recovering `w_j = w_i − 2(rP + xQ)`; and for diagonal `Z` the phases decouple
into three independent LinDistFlows.

Reference: Gan & Low, *Convex relaxations and linear approximation for optimal power flow
in multiphase radial networks*, PSCC 2014, [doi:10.1109/PSCC.2014.7038399](https://doi.org/10.1109/PSCC.2014.7038399).

## The droop block is unchanged

Every droop block here is the same as in the single-phase code, because each only ever
needs one thing from the host: a voltage variable at the inverter's own bus and phase.

```julia
vpv(i, t) = v[PV[i].bus, PV[i].phase, t]   # the whole three-phase interface
```

Nothing in any of the three encodings knows how many phases the network has, or that the
fleet is a mix of four sizes. That is why the three scripts can share their skeleton
verbatim, and why the results agree to the digit.

## Things you may want to change

At the top of `LinDist3Flow_Lambda.jl`:

- `N_PV_PER_PHASE` — sites per phase. With four classes, four per phase gives each phase
  one of each; other counts still cycle through the classes.
- `PV_CLASSES` — the class names and ratings. Add or remove entries freely; the placement,
  the droop curves and the figure all follow. Setting every class to the same rating
  collapses the four curves back to one.
- `VBP`, `QSHAPE` — the droop curve. `QSHAPE` is a free vector, so asymmetric curves (for
  example AS/NZS 4777.2's −0.6 / +0.44) work without code changes.
- `VLIM`, `SBASE_KVA`, `MIP_GAP`.

## A trap worth flagging

The validation sweep must start from a properly rotated three-phase set (1∠0°, 1∠−120°,
1∠+120°). Seeding all three phases at 1∠0° is silent and expensive: the mutual impedance
terms then add instead of largely cancelling, and the sweep reports roughly **twice** the
true voltage deviation — which looks exactly like the linear host being badly wrong. It
cost me a detour; the comment in `sweep()` marks the spot.

## Data and licence

| file | what |
|:--|:--|
| `data/network_5_Feeder_2.bmopf.json` | the feeder, unmodified, in BMOPF JSON format |
| `data/network_5_Feeder_2_report.md` | its network summary, as published |
| `data/load_profile_residential_15min.json` | residential load shape, 96 steps |
| `data/solar_profile.json` | clear-sky irradiance, 96 steps |

The feeder comes from
[BMOPFDraftData](https://github.com/frederikgeth/BMOPFDraftData) (Frederik Geth),
`output/ENWLvariants/Three-wire-Kron-reduced/`, derived from the CSIRO four-wire LV
dataset ([10.25919/jaae-vc35](https://doi.org/10.25919/jaae-vc35)), itself derived from
Electricity North West's LV Network Solutions data.

**Licence: CC BY 4.0**, commercial use permitted, attribution required. Note that the
per-case `meta.license` stamp present on the curated `benchmarks/` cases is *absent* on
these `output/` files — the licence comes from the repository's README table.
