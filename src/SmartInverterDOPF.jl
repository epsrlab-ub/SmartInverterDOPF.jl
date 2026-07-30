"""
    SmartInverterDOPF

Embedding the IEEE 1547 Volt-VAr (Q-V) droop curve of smart inverters into a
distribution optimal power flow.

The package provides one distribution OPF — a current-voltage AC-OPF solved by
successive linearisation — and three interchangeable encodings of the droop law:

| method | class | requires |
|:--|:--|:--|
| `:bigm` | mixed-integer linear | an MILP solver |
| `:lambda` | mixed-integer linear (SOS2) | an MILP solver |
| `:heaviside` | non-smooth nonlinear | an NLP solver |

All three describe the same curve and return the same dispatch; they differ in the
solver technology they need and in how they scale.

```julia
using SmartInverterDOPF, Gurobi
case = load_case()
res  = solve_dopf(case, Gurobi.Optimizer; method = :lambda)
```

!!! warning
    HiGHS does not complete the successive-linearisation loop on the bundled case: it
    solves the first pass, then reports `INFEASIBLE` on the second, for both MILP
    encodings, at every MIP gap tried. Gurobi solves the same sequence to proven
    optimality.
"""
module SmartInverterDOPF

using JuMP
using JSON3

export Case, DroopCurve, DOPFResult,
       load_case, solve_dopf, base_case_voltages,
       ieee1547_curve, droop_q, kWh, nbus, ndg

include("case.jl")
include("droop.jl")
include("dopf.jl")

end # module
