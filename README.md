# SmartInverterDOPF.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://ra-emami.github.io/SmartInverterDOPF.jl/dev/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Embedding IEEE 1547 Volt-VAr droop curves into distribution optimal power flow.

A smart inverter does not accept a reactive-power set-point. It follows a Volt-VAr
curve, deciding from its own terminal voltage how much reactive power to inject or
absorb. An OPF that ignores that curve returns a dispatch the inverter is never going to
deliver. This package puts the curve inside the optimisation — three different ways, on
either of two distribution OPF host models — and shows that the encodings agree.

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
res  = solve_dopf(case, Gurobi.Optimizer; method = :lambda)   # host = :ivacopf (default)

kWh(case, sum(res.PVC))     # PV energy curtailed over the day, kWh
extrema(res.V)              # voltage range across the feeder, p.u.
```

## The two host models

`method` picks how the droop curve is encoded; `host` picks the network model it sits
inside. They are independent — the droop constraints are identical in both hosts.

| `host` | model | accuracy | solve |
|:--|:--|:--|:--|
| `:ivacopf` *(default)* | current-voltage AC-OPF ([doi:10.1109/OJIA.2024.3367547](https://doi.org/10.1109/OJIA.2024.3367547), extended in [doi:10.1016/j.epsr.2026.113613](https://doi.org/10.1016/j.epsr.2026.113613)) | very accurate, near-exact AC | **iterative** — re-linearised until the residual of the exact power-flow identity clears a tolerance |
| `:lindistflow` | linearised branch flow ([doi:10.1109/SMARTGRID.2010.5622021](https://doi.org/10.1109/SMARTGRID.2010.5622021)) | approximate — losses dropped, small voltage deviations assumed | run **once**, much faster and far lower computational effort |

```julia
solve_dopf(case, Gurobi.Optimizer; method = :lambda, host = :lindistflow)
```

> **Solver note.** The committed results are generated with **Gurobi**, which is
> [free for academic users](https://www.gurobi.com/academia/academic-program-and-licenses/).
> We also tried the open-source MILP solvers HiGHS and GLPK; neither worked out — one
> returned an infeasible status inside the successive-linearisation loop, the other was
> too slow to finish. If no MILP licence is available, the `:heaviside` encoding needs
> only Ipopt, which is open source, and reaches the same answer.

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
docs/         Documenter site; builds without any solver
test/         test suite
```

## Regenerating the documentation results

The documentation is built from precomputed results committed under
`docs/src/assets/results/`, so building the docs needs no optimisation solver. To
recompute them:

```bash
julia --project=scripts scripts/generate_results.jl
```

## Building the documentation locally

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## License

MIT — see [LICENSE](LICENSE).
