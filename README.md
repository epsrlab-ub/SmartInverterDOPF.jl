# SmartInverterDOPF.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://ra-emami.github.io/SmartInverterDOPF.jl/dev/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Embedding IEEE 1547 Volt-VAr droop curves into distribution optimal power flow.

A smart inverter does not accept a reactive-power set-point. It follows a Volt-VAr
curve, deciding from its own terminal voltage how much reactive power to inject or
absorb. An OPF that ignores that curve returns a dispatch the inverter is never going to
deliver. This package puts the curve inside the optimisation — three different ways —
and shows that all three agree.

**[Read the tutorial →](https://ra-emami.github.io/SmartInverterDOPF.jl/dev/tutorial_voltvar/)**

## Quick start

```julia
using Pkg
Pkg.add(url = "https://github.com/ra-emami/SmartInverterDOPF.jl")
Pkg.add(["Gurobi", "Ipopt"])
```

```julia
using SmartInverterDOPF, Gurobi

case = load_case()
res  = solve_dopf(case, Gurobi.Optimizer; method = :lambda)

kWh(case, sum(res.PVC))     # PV energy curtailed over the day, kWh
extrema(res.V)              # voltage range across the feeder, p.u.
```

> **Solver note.** Of the three MILP solvers tried, only Gurobi completes the
> successive-linearisation loop on this case. Each open-source solver was given 20×
> Gurobi's total time for the same encoding and withdrawn if it could not finish
> ([`scripts/solver_benchmark.jl`](scripts/solver_benchmark.jl)):
>
> | solver | `:bigm` | `:lambda` |
> |:--|:--|:--|
> | Gurobi | converged, 4 passes, 19.2 s | converged, 6 passes, 15.9 s |
> | HiGHS | `INFEASIBLE` on pass 2 | `INFEASIBLE` on pass 2 |
> | GLPK | withdrawn — pass 1 unfinished at 20× | withdrawn — pass 1 unfinished at 20× |
>
> The model is not infeasible; Gurobi solves the same sequence to proven optimality. The
> solvers return different optimal solutions to the first subproblem, which sends the
> linearisation down different trajectories. Committed results are generated with Gurobi,
> which is [free for academic use](https://www.gurobi.com/academia/academic-program-and-licenses/).
> The `:heaviside` encoding needs no MILP solver at all — only Ipopt, which is open
> source — and reaches the same answer.

## The three encodings

The Volt-VAr law is a five-segment piecewise-linear function — a definition by cases,
which is exactly what a solver cannot read. Three standard rewrites make it tractable:

| `method` | class | idea | needs |
|:--|:--|:--|:--|
| `:bigm` | MILP | one binary per segment activates that segment's voltage window and affine law | MILP solver |
| `:lambda` | MILP | operating point as a convex combination of breakpoints, SOS2 forcing adjacency | MILP solver |
| `:heaviside` | NLP | segment masks from unit steps, summed into one closed-form expression | NLP solver |

All three are exact and, on the bundled case study, return the same dispatch to within
solver tolerance. They differ in the solver technology they demand and in how they scale
with inverters × time steps.

## Case study

The IEEE 33-bus radial feeder over 24 h at 15-minute resolution (96 steps), with three
PV smart inverters at buses 7, 18 and 33, per-class load shapes, a clear-sky irradiance
profile, and a curtailment-minimising objective. The host model is a current-voltage
AC-OPF solved by successive linearisation, with convergence checked against the exact
nonlinear power-flow identity.

## Repository layout

```
src/          the package: case data, droop encodings, DOPF
data/         IEEE 33-bus feeder, load profiles, solar profile (JSON)
scripts/      generate_results.jl — regenerates the precomputed documentation results
              solver_benchmark.jl — measures each MILP solver under the 20x rule
docs/         Documenter site; builds without any solver
test/         test suite
```

## Regenerating the documentation results

The documentation is built from precomputed results committed under
`docs/src/assets/results/`, so building the docs needs no optimisation solver. To
recompute them:

```bash
julia --project=scripts scripts/generate_results.jl          # Gurobi + Ipopt
julia --project=scripts scripts/generate_results.jl highs    # HiGHS + Ipopt (see solver note)
julia --project=scripts scripts/generate_results.jl glpk     # GLPK + Ipopt (see solver note)
```

To reproduce the solver comparison itself — Gurobi, HiGHS and GLPK on both MILP
encodings, with the 20× withdrawal rule:

```bash
julia --project=scripts scripts/solver_benchmark.jl
```

## Building the documentation locally

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## License

MIT — see [LICENSE](LICENSE).
