# Embedding the Volt-VAr droop into a distribution OPF

*Rahmat Emami Mirak*

A smart inverter does not accept a reactive-power set-point. It follows a Volt-VAr
curve: it measures its own terminal voltage and decides, on its own, how much reactive
power to inject or absorb. An optimal power flow that ignores that curve will happily
return a reactive dispatch the inverter is never going to deliver.

This tutorial shows three ways to put the curve *inside* the optimisation, so that every
dispatch point the solver returns is one the inverter would actually produce. All three
are exact — none of them approximates the curve — and on the case study below all three
return the same answer. What separates them is the solver technology they demand and how
they scale.

```@setup tut
using JSON3, Plots, Printf, Markdown
gr(fmt = :svg, size = (760, 420), legendfontsize = 8, titlefontsize = 10,
   guidefontsize = 9, tickfontsize = 8, framestyle = :box, grid = true,
   gridalpha = 0.15, dpi = 150)

RES  = joinpath("assets", "results")
case = JSON3.read(read(joinpath(RES, "case.json"), String))
runs = Dict(m => JSON3.read(read(joinpath(RES, "$m.json"), String))
            for m in ("bigm", "lambda", "heaviside"))

const NAMES = Dict("bigm" => "Big-M", "lambda" => "Lambda / SOS2", "heaviside" => "Heaviside")
const ORDER = ["bigm", "lambda", "heaviside"]
hours = range(0, 24 - 24/96, length = 96)          # 96 quarter-hourly steps
fmt(x, n) = @sprintf("%.*f", n, x)
sci(x)    = @sprintf("%.2e", x)
md(rows...) = Markdown.parse(join(rows, "\n"))

# Table builders. These live here rather than in the visible blocks so that the page
# shows tables, not the string-mangling that produces them — while still deriving every
# number from the committed results rather than hard-coding it.
breakpoint_table() = md(
    "| | " * join(["``V^{\\text{bp}}_$i``" for i in 1:6], " | ") * " |",
    "|:--|" * repeat("--:|", 6),
    "| voltage (p.u.) | " * join(fmt.(Vbp, 2), " | ") * " |",
    "| ``q/\\bar q`` | " * join(fmt.(qshape, 0), " | ") * " |")

inverter_table() = md(
    "| bus | PV rating (kW) | inverter rating (kVA) | reactive capability ``\\bar q`` (kVAr) |",
    "|--:|--:|--:|--:|",
    join(["| $(case.DG_SET[i]) | $(fmt(case.Pdg_max_kW[i], 0)) | " *
          "$(fmt(case.Sdg_max_kVA[i], 0)) | $(fmt(case.qbar_kVAr[i], 0)) |"
          for i in eachindex(case.DG_SET)], "\n"))

iteration_table(m) = md(
    "| iteration | solve (s) | objective | linearisation residual | solver status |",
    "|--:|--:|--:|--:|:--|",
    join(["| $(r.iter) | $(fmt(r.seconds, 2)) | $(fmt(r.objective, 6)) | " *
          "$(sci(r.residual)) | `$(r.status)` |" for r in runs[m].iterations], "\n"))

deviation_table() = md(
    "| method | max ``\\lvert q^G_i - q_i(v_i)\\rvert`` (p.u.) |",
    "|:--|--:|",
    join(["| $(NAMES[m]) | $(sci(runs[m].max_droop_deviation)) |" for m in ORDER], "\n"))

comparison_table() = md(
    "| method | class | solver | variables | binaries | constraints | iters | " *
    "solve (s) | curtailed (kWh) | curtailed (%) | losses (kWh) | voltage range (p.u.) |",
    "|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|",
    join(map(ORDER) do m
        r = runs[m]
        "| $(NAMES[m]) | $(r.nbin == 0 ? "NLP" : "MILP") | `$(r.solver)` | $(r.nvar) | " *
        "$(r.nbin) | $(r.ncon) | $(r.n_iterations) | $(fmt(r.solve_seconds, 1)) | " *
        "$(fmt(r.E_curt_kWh, 1)) | $(fmt(r.curt_percent, 3)) | $(fmt(r.loss_kWh, 1)) | " *
        "$(fmt(r.Vmin, 4)) – $(fmt(r.Vmax, 4)) |"
    end, "\n"))

base_case_sentence() = md(
    "**$(fmt(case.V_base_lo, 4)) p.u.** at the ends of the feeder, against a limit of " *
    "$(fmt(case.Vmin_limit, 2)) p.u.")

Vbp, qshape = collect(Float64, case.Vbp), collect(Float64, case.qshape)
```

## Prerequisites

Julia with [JuMP](https://jump.dev), a mathematical-programming solver, and this
package. Two of the three methods produce a mixed-integer linear program and need an
MILP solver; the third produces a nonlinear program and needs an NLP solver.

```julia
using Pkg
Pkg.add(["JuMP", "Gurobi", "Ipopt"])
Pkg.add(url = "https://github.com/ra-emami/SmartInverterDOPF.jl")
```

!!! warning "Solver choice is not free here"
    The results on this page were produced with **Gurobi** (MILP) and **Ipopt** (NLP).
    The open-source MILP solver HiGHS does *not* complete the successive-linearisation
    loop on this case: it solves the first pass, then reports `INFEASIBLE` on the second,
    for both `:bigm` and `:lambda`, at every MIP gap from `1e-3` down to `0`.

    The model is not infeasible — Gurobi solves the same sequence to proven optimality.
    The likely cause is that the two solvers return *different* optimal solutions to the
    first subproblem, which sends the linearisation down different trajectories; one of
    them lands on a subproblem that is genuinely infeasible. That is a fragility of this
    successive-linearisation scheme, not a defect in either solver, and it is worth
    knowing about before building on it.

## Why the curve has to live inside the OPF

IEEE 1547-2018 requires every interconnecting distributed energy resource to be capable
of Volt-VAr control. The utility enables the function and sets the curve; the inverter
then runs it autonomously as a local feedback law. An advanced distribution management
system can coordinate hundreds of these inverters through a distribution OPF — but only
if the OPF knows the law each one is following.

Leave the curve out and the OPF treats ``q_i^G`` as a free decision variable inside the
inverter's apparent-power circle. It will pick whatever value minimises the objective.
The inverter, meanwhile, is looking at its own terminal voltage and producing something
else entirely. The dispatch is not merely suboptimal; it is not physically realisable.

Put the curve in, and the feasible set shrinks to exactly the points the fleet can
actually reach. As a bonus, once the curve is an algebraic object inside the model, its
breakpoints can themselves become decision variables — which is how droop curves get
optimised rather than merely respected.

## The IEEE 1547 Volt-VAr law

The curve is piecewise linear in five segments, defined by six breakpoint voltages
``V^{\text{bp}}_1 \le \dots \le V^{\text{bp}}_6`` and the reactive set-point at each. Writing
``\bar q_i`` for the reactive capability of inverter ``i``:

```math
q_i(v_i) \;=\;
\begin{cases}
\bar q_i, & V^{\text{bp}}_1 \le v_i \le V^{\text{bp}}_2 \quad\text{(full injection)}\\[4pt]
\bar q_i \dfrac{V^{\text{bp}}_3 - v_i}{V^{\text{bp}}_3 - V^{\text{bp}}_2},
        & V^{\text{bp}}_2 \le v_i \le V^{\text{bp}}_3 \quad\text{(sloped)}\\[6pt]
0, & V^{\text{bp}}_3 \le v_i \le V^{\text{bp}}_4 \quad\text{(dead-band)}\\[4pt]
-\bar q_i \dfrac{v_i - V^{\text{bp}}_4}{V^{\text{bp}}_5 - V^{\text{bp}}_4},
        & V^{\text{bp}}_4 \le v_i \le V^{\text{bp}}_5 \quad\text{(sloped)}\\[6pt]
-\bar q_i, & V^{\text{bp}}_5 \le v_i \le V^{\text{bp}}_6 \quad\text{(full absorption)}
\end{cases}
```

Low voltage means inject reactive power to hold the voltage up; high voltage means
absorb it. Between the two sits a dead-band in which the inverter does nothing, so that
small fluctuations do not provoke needless reactive flow.

```@example tut
Vbp, qshape = collect(Float64, case.Vbp), collect(Float64, case.qshape)
plot(Vbp, qshape, lw = 3, color = :steelblue, label = "IEEE 1547 Q-V droop",
     xlabel = "terminal voltage  v  (p.u.)", ylabel = "q / q̄",
     title = "The Volt-VAr characteristic", ylims = (-1.35, 1.35), legend = :topright)
scatter!(Vbp, qshape, ms = 5, color = :steelblue, label = "breakpoints")
vspan!([Vbp[3], Vbp[4]], color = :grey, alpha = 0.12, label = "dead-band")
hline!([0], color = :black, lw = 0.6, alpha = 0.5, label = false)
annotate!([(Vbp[1] + 0.007, 1.13, text("inject", 8, :left, :steelblue)),
           (Vbp[6] - 0.007, -1.13, text("absorb", 8, :right, :steelblue))])
```

The curve used throughout this tutorial is the package default:

```@example tut
breakpoint_table()   # hide
```

## Why an "if-else" cannot go straight into a solver

The law above is a definition by cases, and that is exactly what a solver cannot read.
Two distinct obstacles follow.

**Conditional logic.** Which of the five expressions applies depends on ``v_i``, which is
itself a decision variable. Branching on the value of an unknown is not an algebraic
constraint, and an algebraic constraint is the only thing a solver accepts.

**Non-differentiability.** Even setting the branching aside, the slope jumps at every
breakpoint. Newton and interior-point methods build their steps from derivatives, and at
a kink the derivative does not exist.

There are two ways out, and they define the rest of this tutorial:

- **Introduce integer variables** to encode the logic exactly. The model becomes an
  MILP. This is the Big-M and Lambda/SOS2 route.
- **Write the logic in closed algebraic form** using step functions. The model stays
  integer-free but becomes non-smooth, so it needs an NLP solver. This is the Heaviside
  route.

## The host model

The droop is a self-contained module. Whatever distribution OPF you use, it exposes a
voltage magnitude at each inverter bus; the droop module adds the relationship tying
that inverter's reactive output to that voltage:

```
        ┌────────────────────────────┐
        │            DOPF            │
        │  network model + limits    │
        └───────┬────────────▲───────┘
         exposes│ vᵢ      qᵢᴳ│ sets
        ┌───────▼────────────┴───────┐
        │      Q-V droop module      │
        │   the IEEE 1547 curve      │
        └────────────────────────────┘
```

Nothing in the three encodings below depends on the host. They work equally with
LinDistFlow, a current-voltage AC-OPF, an unbalanced three-phase AC-OPF, or a full
nonlinear AC-OPF — they only require that the host expose ``v_i`` at each inverter bus.

The host used here is a **current-voltage AC optimal power flow** (IVACOPF), in which the
network is written in rectangular current and voltage coordinates. Ohm's law and KCL are
then exactly linear. Two things remain bilinear — the ``v \cdot I`` power balance and the
``|I|^2`` branch loss — and these are handled by successive linearisation: each is
expanded about the previous iterate and the model is re-solved until the residual of the
*exact* loss identity falls below a tolerance. Convergence is checked against the true
nonlinear relation, not the linearised one, so the converged point satisfies the real
power flow.

## The setup

The IEEE 33-bus radial feeder, over a full day at 15-minute resolution — 96 time steps.
Three PV systems with smart inverters sit at buses 7, 18 and 33. Loads follow separate
industrial, commercial and residential shapes; PV follows a clear-sky irradiance profile.

```@example tut
inverter_table()   # hide
```

Each inverter is rated 10 % above its PV array, so there is headroom for reactive support
even at full irradiance. Bus voltages are limited to
``[`` `` `` ``0.95, 1.05`` ``]`` p.u., and the objective is to **minimise total PV
curtailment** over the day.

```@example tut
avail = [collect(Float64, a) for a in case.avail_kW]
plot(hours, sum(avail), lw = 2, color = :darkorange, fillrange = 0, fillalpha = 0.18,
     label = "total PV available", xlabel = "hour of day", ylabel = "kW",
     title = "Available PV generation across the three inverters",
     xticks = 0:3:24, xlims = (0, 24), legend = :topleft)
for (i, b) in enumerate(case.DG_SET)
    plot!(hours, avail[i], lw = 1.4, ls = :dash, label = "bus $b")
end
plot!()
```

Why does a curtailment objective have anything to do with voltage at all? Because the
droop ties the two together. Active power injection raises the local voltage; the droop
reads that voltage and sets reactive power accordingly; and reactive flow moves voltages
across the whole feeder. The optimiser wants every kilowatt it can get — the droop
decides what taking it costs everywhere else.

## Method A — Big-M

Introduce a binary ``\delta_{b}`` for each of the five segments, exactly one of which is
active:

```math
\sum_{b=1}^{5}\delta_{b}=1 .
```

Each binary must switch on its own voltage window. For a segment ``b`` spanning
``[V^{\text{bp}}_{b}, V^{\text{bp}}_{b+1}]``, the pair

```math
v_i \;\ge\; V^{\text{bp}}_{b} - M\,(1-\delta_{b}), \qquad
v_i \;\le\; V^{\text{bp}}_{b+1} + M\,(1-\delta_{b})
```

is vacuous when ``\delta_b = 0`` and binding when ``\delta_b = 1``.

The sloped segments need one more step. Their law involves the product ``\delta_b v_i``
of a binary and a continuous variable, which is bilinear. Because ``\delta_b`` is binary
and ``v_i`` is bounded, that product is linearised **exactly** by an auxiliary variable
``W_b := \delta_b v_i``:

```math
-M(1-\delta_b) \;\le\; v_i - W_b \;\le\; M(1-\delta_b), \qquad
V^{\text{bp}}_{b}\,\delta_b \;\le\; W_b \;\le\; V^{\text{bp}}_{b+1}\,\delta_b .
```

The second pair doubles as the voltage window, so no separate window constraint is
needed for sloped segments. The droop law is then a single linear equation in which
every coefficient is a constant:

```math
q_i^G = \delta_1\bar q_i
      + \alpha_1 W_2 + \delta_2\frac{\bar q_i V^{\text{bp}}_3}{V^{\text{bp}}_3 - V^{\text{bp}}_2}
      + \alpha_2 W_4 + \delta_4\frac{\bar q_i V^{\text{bp}}_4}{V^{\text{bp}}_5 - V^{\text{bp}}_4}
      - \delta_5\bar q_i ,
```

with slopes ``\alpha_1 = -\bar q_i/(V^{\text{bp}}_3-V^{\text{bp}}_2)`` and
``\alpha_2 = -\bar q_i/(V^{\text{bp}}_5-V^{\text{bp}}_4)``.

In JuMP:

```julia
@variable(model, δ[1:5, DG_SET, HOUR_SET, QUARTER_SET], Bin)
@variable(model, W2[DG_SET, HOUR_SET, QUARTER_SET])
@variable(model, W4[DG_SET, HOUR_SET, QUARTER_SET])

@constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
    sum(δ[i, d, h, m] for i in 1:5) == 1)

# flat segments 1, 3, 5 — the binary only switches on a voltage window
for (i, lo, hi) in ((1, 1, 2), (3, 3, 4), (5, 5, 6))
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        v[d, h, m] >= Vbp[lo] - Mbig * (1 - δ[i, d, h, m]))
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        v[d, h, m] <= Vbp[hi] + Mbig * (1 - δ[i, d, h, m]))
end

# sloped segments 2 and 4 — W = δ·v, whose bounds double as the window
for (i, W, lo, hi) in ((2, W2, 2, 3), (4, W4, 4, 5))
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        v[d, h, m] - W[d, h, m] >= -Mbig * (1 - δ[i, d, h, m]))
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        v[d, h, m] - W[d, h, m] <=  Mbig * (1 - δ[i, d, h, m]))
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        W[d, h, m] >= Vbp[lo] * δ[i, d, h, m])
    @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
        W[d, h, m] <= Vbp[hi] * δ[i, d, h, m])
end
```

!!! tip "Choose M as tightly as you can justify"
    ``M`` only has to dominate the largest possible violation of a deactivated
    constraint, which here is set by the voltage bounds. A needlessly large ``M`` leaves
    the LP relaxation loose, the branch-and-bound tree deep, and the solve slow. The
    value used here is `1.1`.

The cost of exactness is bookkeeping: five binaries per inverter per time step, plus two
auxiliary continuous variables. And if the breakpoints themselves become decision
variables, ``W_b`` turns into a product of two continuous unknowns and the reformulation
is no longer exact — the reason breakpoint optimisation is usually done in the Lambda
formulation instead.

## Method B — Lambda / SOS2

Any point on a line segment is a weighted average of its two endpoints. So give every
breakpoint a weight ``\lambda_b`` and write *both* coordinates with the same weights:

```math
v_i = \sum_{b=1}^{6}\lambda_b V^{\text{bp}}_b, \qquad
q_i^G = \sum_{b=1}^{6}\lambda_b q^{\text{bp}}_b, \qquad
\sum_{b=1}^{6}\lambda_b = 1, \qquad \lambda_b \ge 0 .
```

Because one set of weights builds both coordinates, ``v_i`` and ``q_i^G`` move together
along the curve. The breakpoint ordinates ``q^{\text{bp}}_b`` are constants, so these are
linear constraints.

There is a catch. Nothing above stops the solver from blending *non-adjacent*
breakpoints — mixing ``\lambda_1`` and ``\lambda_5``, say — which produces points in the
interior of the convex hull rather than on the curve. The fix is the classical
**SOS2** condition: at most two weights may be nonzero, and they must be adjacent.
In mixed-integer form, with one binary ``z_b`` per segment:

```math
\sum_{b=1}^{5} z_b = 1, \qquad
\lambda_1 \le z_1, \qquad
\lambda_b \le z_{b-1} + z_b \;\; (b = 2,\dots,5), \qquad
\lambda_6 \le z_5 .
```

Exactly one segment is active, and a weight may be nonzero only if it touches that
segment — so precisely the two ends of the active segment blend, and nothing else.

```julia
@variable(model, λ[1:6, DG_SET, HOUR_SET, QUARTER_SET] >= 0)
@variable(model, z[1:5, DG_SET, HOUR_SET, QUARTER_SET], Bin)

@constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
    sum(λ[i, d, h, m] for i in 1:6) == 1)
@constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
    sum(z[i, d, h, m] for i in 1:5) == 1)
@constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
    λ[1, d, h, m] <= z[1, d, h, m])
@constraint(model, [i in 2:5, d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
    λ[i, d, h, m] <= z[i-1, d, h, m] + z[i, d, h, m])
@constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
    λ[6, d, h, m] <= z[5, d, h, m])

@constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
    v[d, h, m] == sum(λ[i, d, h, m] * Vbp[i] for i in 1:6))
@constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
    Qdg[d, h, m] == sum(λ[i, d, h, m] * qpts[d][i] for i in 1:6))
```

!!! note "SOS2 without the binaries"
    Most MILP solvers support SOS2 natively via `MOI.SOS2`, which lets the solver
    branch on the set directly instead of on explicit binaries. The formulation above is
    written out longhand because it is portable and because it makes the logic visible —
    which is the point of a tutorial.

The Lambda form has a decisive practical advantage over Big-M: the breakpoint voltages
``V^{\text{bp}}_b`` appear *linearly*. Make them decision variables and the constraint
``v_i = \sum_b \lambda_b V^{\text{bp}}_b`` becomes bilinear in ``\lambda`` and
``V^{\text{bp}}`` — still a single well-understood bilinear term, rather than the tangle
Big-M produces. This is why adaptive-droop work is normally built on Lambda.

## Method C — Heaviside

The third route abandons integers entirely. An "if" is just an on/off switch, and the
unit step is exactly that:

```math
H(x) = \begin{cases} 1, & x \ge 0\\ 0, & x < 0\end{cases}
```

Shift it to flip at a breakpoint and subtract two of them, and you get a **window** that
equals 1 on one segment and 0 everywhere else:

```math
\mathcal{W}_b(v_i) \;=\; H\!\left(v_i - V^{\text{bp}}_{b}\right) - H\!\left(v_i - V^{\text{bp}}_{b+1}\right)
```

which is precisely the condition ``V^{\text{bp}}_b \le v_i \le V^{\text{bp}}_{b+1}`` — the
if-else of segment ``b``, written without logic and without binaries. Multiply each
segment's law by its own window and add them up. The windows are disjoint, so at any
voltage all but one vanish and the sum collapses to the single active law:

```math
q_i^G \;=\; \bar q_i\,\mathcal{W}_1
      \;+\; \alpha_1\!\left(v_i - V^{\text{bp}}_3\right)\mathcal{W}_2
      \;+\; \alpha_2\!\left(v_i - V^{\text{bp}}_4\right)\mathcal{W}_4
      \;-\; \bar q_i\,\mathcal{W}_5
```

with the sloped terms anchored at their zero crossings, so that segment 2 gives
``\bar q_i`` at ``V^{\text{bp}}_2`` and ``0`` at ``V^{\text{bp}}_3`` — exactly the sloped law
of the curve. The dead-band contributes nothing and needs no term at all.

```julia
Hstep(x) = op_ifelse(op_greater_than_or_equal_to(x, 0), 1.0, 0.0)
α1 = Dict(d => -qbar[d] / (Vbp[3] - Vbp[2]) for d in DG_SET)
α2 = Dict(d => -qbar[d] / (Vbp[5] - Vbp[4]) for d in DG_SET)

@constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
    Qdg[d, h, m] ==
        qbar[d] * (Hstep(v[d,h,m] - Vbp[1]) - Hstep(v[d,h,m] - Vbp[2]))
      + α1[d] * (v[d,h,m] - Vbp[3]) * (Hstep(v[d,h,m] - Vbp[2]) - Hstep(v[d,h,m] - Vbp[3]))
      + α2[d] * (v[d,h,m] - Vbp[4]) * (Hstep(v[d,h,m] - Vbp[4]) - Hstep(v[d,h,m] - Vbp[5]))
      - qbar[d] * (Hstep(v[d,h,m] - Vbp[5]) - Hstep(v[d,h,m] - Vbp[6])))
```

`op_ifelse` and `op_greater_than_or_equal_to` are JuMP's nonlinear operators
(JuMP ≥ 1.15); they build the expression correctly outside a macro.

No extra variables at all — one algebraic expression per inverter per time step. The
price is paid in solver behaviour, and it is not merely theoretical. ``H(\cdot)`` is
discontinuous, so the derivative is undefined at every breakpoint and the problem is
non-convex. On this case study Ipopt does reach the same answer as the MILP encodings,
but the first pass — started from a flat voltage profile, far from the solution — stops
at `ALMOST_LOCALLY_SOLVED`, short of its own convergence tolerance:

```@example tut
iteration_table("heaviside")   # hide
```

Once the linearisation point is close enough that most inverters sit comfortably inside a
single segment, the later passes converge cleanly. But that first pass is the
non-smoothness showing up in practice: an interior-point method cannot get reliable
derivative information at the kinks, and a harder case — more inverters, voltages sitting
nearer the breakpoints — is where this bites. In production this encoding is normally
*smoothed*, the step replaced by a sigmoid ``H(x) \approx (1 + e^{-kx})^{-1}``, which
restores differentiability at the cost of no longer representing the curve exactly.

## Verification: does the dispatch actually lie on the curve?

This is the check that matters. Each method is exact only if every one of the
``3 \times 96 = 288`` optimised operating points lands on the droop.

```@example tut
plts = map(ORDER) do m
    r = runs[m]
    p = plot(Vbp, qshape, lw = 2.5, color = :steelblue, label = "Q-V droop",
             title = NAMES[m], xlabel = "v at inverter bus (p.u.)",
             ylabel = m == "bigm" ? "q / q̄" : "", ylims = (-1.35, 1.35),
             xlims = (0.87, 1.11), legend = :topright)
    scatter!(p, collect(Float64, r.Vdg_flat), collect(Float64, r.Qdg_norm),
             m = :+, ms = 4, msw = 1.6, color = :crimson, label = "dispatch")
    p
end
plot(plts..., layout = (1, 3), size = (1000, 350), left_margin = 4Plots.mm,
     bottom_margin = 4Plots.mm)
```

Every point sits on the line. Numerically:

```@example tut
deviation_table()   # hide
```

All three are at the level of floating-point round-off, which is what "exact" means here:
the encodings do not approximate the curve, they reproduce it.

Note that the operating points cluster in the sloped region below nominal and in the
dead-band. The saturated tails are never reached — on this feeder the voltage never drops
to ``V^{\text{bp}}_2`` nor rises to ``V^{\text{bp}}_5``. The flat segments still have to be
in the model, because the solver must be free to consider them, but they do no work here.

## Side by side

```@example tut
comparison_table()   # hide
```

The three rows agree on every physical quantity. Curtailed energy matches to four
significant figures, losses to five, and the voltage range is identical to four
decimals.
The residual differences are the NLP solver's convergence tolerance, not a modelling
difference — which is the empirical statement of the claim that these are three
encodings of one curve.

What differs is the machinery. Heaviside adds no variables at all; Big-M and Lambda each
add 1440 binaries — five per inverter per time step — and that count grows with
inverters × time steps, which is the scaling wall for the integer methods. Big-M reached
its answer in fewer successive-linearisation passes here, but iteration counts of this
kind are case-specific and should not be read as a general ranking.

!!! warning "Solve times are indicative only"
    These timings come from a single run on one machine with one solver configuration,
    and the successive-linearisation loop rebuilds the model from scratch each pass. Use
    them to compare orders of magnitude, not to rank solvers.

### What the curtailment figure is, and what it is not

The curtailed-energy column is reported here as the output of this host model on this
case — it is the quantity the three encodings are being compared *on*, and they agree on
it. It should not be read as a physical curtailment requirement of the feeder, because
the usual suspects turn out not to be binding:

| probe | outcome |
|:--|:--|
| voltage limits widened from ``[0.95, 1.05]`` to ``[0.90, 1.10]`` | curtailment unchanged |
| voltage limits effectively removed, ``[0.50, 1.50]`` | curtailment unchanged |
| MIP gap tightened from ``10^{-3}`` to ``0`` | identical to two decimals |
| inverter capability at the most-curtailed step | ``\lvert S\rvert / S_{\max} \approx 0.38`` |
| feeder voltage at the most-curtailed step | minimum ``0.983`` p.u., far inside limits |

At the worst-curtailed quarter-hour the bus-7 inverter holds about 915 kW of 2200 kW
available, with zero reactive output and no constraint active at its own bus or anywhere
else on the feeder. Whatever sets that level, it is a property of the host formulation
rather than of the droop encoding — which is precisely why it does not disturb the
comparison above: all three encodings sit inside the same host and inherit it equally.

Because the encodings are the subject here, the host is taken as given. Anyone building
on this model for curtailment studies should chase that down first.

## What the droop actually buys

The comparison so far has been between encodings. The more useful comparison is against
not having the inverters at all. The base case is the same feeder and the same demand
with no smart inverters, solved by a backward/forward sweep power flow.

```@example tut
buses = 1:case.n_bus
p = plot(xlabel = "bus", ylabel = "voltage (p.u.)", legend = :bottomleft,
         title = "Daily voltage envelope, with and without smart inverters",
         xticks = [1, 5, 10, 15, 20, 25, 30, 33], xlims = (1, 33))
plot!(p, buses, collect(Float64, case.V_base_max), lw = 2, ls = :dash,
      color = :grey65, label = "no inverters — max")
plot!(p, buses, collect(Float64, case.V_base_min), lw = 2, ls = :dash,
      color = :grey35, label = "no inverters — min")
plot!(p, buses, collect(Float64, runs["lambda"].V_with_max), lw = 2.2,
      color = :darkorange2, label = "with inverters — max")
plot!(p, buses, collect(Float64, runs["lambda"].V_with_min), lw = 2.2,
      color = :dodgerblue4, label = "with inverters — min")
hline!(p, [case.Vmin_limit, case.Vmax_limit], ls = :dot, lw = 1.5, color = :red,
       label = "limits")
```

Without the inverters the feeder violates its lower voltage limit: the minimum across the
day reaches

```@example tut
base_case_sentence()   # hide
```

With the droop-controlled inverters the whole envelope sits inside the band: the
worst-case voltage rises to 0.9501 p.u. and the feeder complies.

The reactive support is not a marginal improvement here — it is what makes the feeder
operable at all. Rerun the same case with the curve flattened to ``q \equiv 0``, so the
inverters still deliver active power but provide no reactive support, and the OPF is
**infeasible** at these voltage limits. There is no active-power dispatch that keeps this
feeder inside ``[0.95, 1.05]`` without Volt-VAr control.

Following a single inverter through the day shows the mechanism:

```@example tut
r  = runs["lambda"]
i  = 1                                  # the inverter at bus 7
b  = case.DG_SET[i]
V  = collect(Float64, r.Vdg_series[i])
Q  = collect(Float64, r.Qdg_series[i])
P  = collect(Float64, r.Pdg_series[i])
A  = collect(Float64, case.avail_kW[i])

p1 = plot(hours, V, lw = 2, color = :dodgerblue4, label = "v at bus $b",
          ylabel = "voltage (p.u.)", xticks = 0:3:24, xlims = (0, 24), legend = :bottomright)
hline!(p1, [case.Vmin_limit], ls = :dot, color = :red, label = "lower limit")
hline!(p1, [Vbp[3], Vbp[4]], ls = :dash, color = :grey55, alpha = 0.8,
       label = "dead-band edges")

p2 = plot(hours, Q, lw = 2, color = :seagreen, label = "reactive output",
          ylabel = "kVAr", xticks = 0:3:24, xlims = (0, 24), legend = :topright)
hline!(p2, [0], color = :black, lw = 0.6, alpha = 0.5, label = false)

p3 = plot(hours, A, lw = 2, ls = :dash, color = :grey45, label = "available",
          xlabel = "hour of day", ylabel = "kW", xticks = 0:3:24, xlims = (0, 24),
          legend = :topleft)
plot!(p3, hours, P, lw = 2, color = :darkorange2, fillrange = 0, fillalpha = 0.15,
      label = "delivered")

plot(p1, p2, p3, layout = (3, 1), size = (780, 700), link = :x,
     left_margin = 5Plots.mm)
```

The inverter injects reactive power whenever its voltage sits below the dead-band, backs
off to zero inside it, and the active power it delivers tracks the available irradiance
except where the voltage limit forces a cut.

## Reproducing these results

The figures and tables on this page are drawn from results committed to the repository,
so building the documentation needs no solver. To regenerate them:

```bash
julia --project=scripts scripts/generate_results.jl
```

That runs all three methods with Gurobi and Ipopt and rewrites
`docs/src/assets/results/`. To run a single method yourself:

```julia
using SmartInverterDOPF, Gurobi

case = load_case()
res  = solve_dopf(case, Gurobi.Optimizer; method = :lambda)

println("curtailed: ", round(kWh(case, sum(res.PVC)), digits = 2), " kWh")
println("voltage:   ", round(minimum(res.V), digits = 4), " – ",
                       round(maximum(res.V), digits = 4), " p.u.")
```

Swap `:lambda` for `:bigm` or `:heaviside`; the latter needs an NLP solver such as
`Ipopt.Optimizer`.

## Takeaways

**Embedding is a correctness requirement, not a refinement.** Smart inverters follow
their curve, not a set-point. Only a droop-aware OPF returns a dispatch the fleet will
actually deliver.

**Two exact families, one curve.** Integer encodings (Big-M, Lambda/SOS2) give an MILP;
non-smooth algebra (Heaviside) gives an NLP. All three reproduce the curve to
round-off and return the same dispatch. The choice is which solver world you want to
work in.

**Scale picks the method.** Binaries multiply with inverters × time steps, which is what
eventually breaks the MILP route on large fleets. The integer-free encoding avoids that
but hands the difficulty to the NLP solver, where non-smoothness shows up as degraded
convergence — visible here in a first pass that stops at `ALMOST_LOCALLY_SOLVED`.

**Optimising the breakpoints is the open frontier.** Treat ``V^{\text{bp}}`` as decision
variables and the Lambda formulation acquires a single bilinear term, while Big-M's exact
product linearisation collapses. That is why adaptive-droop formulations are normally
built on Lambda.

## References

**Standard and host models**

1. IEEE Std 1547-2018, *IEEE Standard for Interconnection and Interoperability of
   Distributed Energy Resources with Associated Electric Power Systems Interfaces*.
2. M. E. Baran and F. F. Wu, "Network reconfiguration in distribution systems for loss
   reduction and load balancing," *IEEE Transactions on Power Delivery*, 1989.
   — LinDistFlow
3. Z. Soltani, M. Khorsand, and S. Ma, "Current–Voltage Unbalanced Distribution AC
   Optimal Power Flow for Advanced Distribution Management System Applications,"
   *IEEE Open Journal of Industry Applications*, 2024. — the IVACOPF host model used here

**Embedding the droop in a distribution OPF**

4. A. Savasci, A. Inaolaji, and S. Paudyal, "Distribution Grid Optimal Power Flow
   Integrating Volt-VAr Droop of Smart Inverters," *IEEE Green Technologies Conference*, 2021. — Big-M
5. A. Inaolaji, A. Savasci, and S. Paudyal, "Distribution Grid Optimal Power Flow with
   Volt-VAr and Volt-Watt Settings of Smart Inverters," *IEEE IAS Annual Meeting*, 2021.
   — Lambda / SOS2
6. A. Inaolaji, A. Savasci, and S. Paudyal, "Distribution Grid OPF in Unbalanced
   Multiphase Networks with Volt-VAr and Volt-Watt Droop Settings of Smart Inverters,"
   *IEEE Transactions on Industry Applications*, 2022. — Lambda, three-phase
7. A. Inaolaji, A. Savasci, and S. Paudyal, "Optimal Droop Settings of Smart Inverters,"
   *IEEE Photovoltaic Specialists Conference (PVSC)*, 2021. — Heaviside
8. R. Emami Mirak and A. Inaolaji, "Adaptive and Fair Optimization of Smart Inverter
   Droop Curves in Distribution Grids," *Electric Power Systems Research*, vol. 262,
   2027, Art. no. 113613. — Lambda / SOS2 with adaptive breakpoints
