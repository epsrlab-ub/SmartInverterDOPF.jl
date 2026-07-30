# SmartInverterDOPF.jl

*Embedding IEEE 1547 Volt-VAr droop curves into distribution optimal power flow.*

A smart inverter does not take a reactive-power set-point — it follows a Volt-VAr curve
based on its own terminal voltage. An OPF that ignores that curve returns a dispatch the
inverter will never deliver. This package puts the curve inside the optimisation, three
different ways, and shows that they agree.

```julia
using SmartInverterDOPF, Gurobi

case = load_case()
res  = solve_dopf(case, Gurobi.Optimizer; method = :lambda)

kWh(case, sum(res.PVC))     # PV energy curtailed over the day
extrema(res.V)              # voltage range across the feeder
```

## The three encodings

The Volt-VAr law is a five-segment piecewise-linear function — a definition by cases,
which is exactly what a solver cannot read. There are three standard ways to rewrite it
as constraints a solver accepts:

| `method` | class | idea | needs |
|:--|:--|:--|:--|
| `:bigm` | MILP | one binary per segment activates that segment's voltage window and affine law | MILP solver |
| `:lambda` | MILP | the operating point is a convex combination of breakpoints, with SOS2 forcing adjacency | MILP solver |
| `:heaviside` | NLP | segment masks built from unit steps, summed into one closed-form expression | NLP solver |

All three are exact — they reproduce the curve rather than approximating it — and on the
bundled case study they return the same dispatch to within solver tolerance. They differ
in the solver technology they demand and in how they scale with inverters × time steps.

The [Tutorial](@ref "Embedding the Volt-VAr droop into a distribution OPF") derives all
three, verifies that every optimised operating point lands on the curve, and compares
them side by side.

## What's in the box

- A **current-voltage AC-OPF** (IVACOPF) host model, solved by successive linearisation,
  with convergence checked against the exact nonlinear power-flow identity.
- The **IEEE 33-bus** feeder over 24 h at 15-minute resolution, with per-class load
  shapes and a clear-sky PV profile.
- Three smart inverters, an inverter capability polygon, and a curtailment-minimising
  objective.
- A **backward/forward sweep** power flow for the no-inverter reference case.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/ra-emami/SmartInverterDOPF.jl")
```

You will also need a solver. `:bigm` and `:lambda` need an MILP solver; `:heaviside`
needs an NLP solver such as [Ipopt](https://github.com/jump-dev/Ipopt.jl).

!!! warning "MILP solver choice"
    The results throughout this documentation were produced with **Gurobi**. HiGHS does
    not complete the successive-linearisation loop on this case — see the note in the
    [Tutorial](@ref "Embedding the Volt-VAr droop into a distribution OPF").

## Citing

If this material is useful in your work, please cite the papers it builds on, listed in
the tutorial's [References](@ref) section.
