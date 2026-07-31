# Measure what each MILP solver actually does on this case.
#
#   julia --project=scripts scripts/solver_benchmark.jl
#
# Runs :bigm and :lambda under Gurobi, HiGHS and GLPK and writes one JSON per
# (solver, method) into docs/src/assets/results/benchmark/. Gurobi goes first and sets
# the reference time; every open-source solver then gets a budget of WITHDRAW_FACTOR
# times Gurobi's total solve time for the same encoding. A solver that cannot finish
# inside that budget is recorded as `withdrawn` — the documentation reports the budget
# it was given and the iteration it was on when the clock ran out, rather than leaving
# it running or quietly dropping it.
#
# Every outcome is recorded, including failures: `solve_dopf` raises on any terminal
# status outside OPTIMAL / LOCALLY_SOLVED / ALMOST_LOCALLY_SOLVED, and the message it
# raises names the iteration and the status, which is what gets written here.

using SmartInverterDOPF
using JSON3

const ROOT   = dirname(@__DIR__)
const OUTDIR = joinpath(ROOT, "docs", "src", "assets", "results", "benchmark")
mkpath(OUTDIR)

# a solver is withdrawn once it has burned this many times Gurobi's total solve time
const WITHDRAW_FACTOR = 20

using Gurobi, HiGHS, GLPK

solver_setup = Dict(
    "Gurobi" => (Gurobi.Optimizer, Dict{String,Any}("MIPGap" => 0.001, "OutputFlag" => 0)),
    "HiGHS"  => (HiGHS.Optimizer,  Dict{String,Any}("mip_rel_gap" => 0.001, "output_flag" => false)),
    "GLPK"   => (GLPK.Optimizer,   Dict{String,Any}("mip_gap" => 0.001, "msg_lev" => 0)),
)

case = load_case(joinpath(ROOT, "data"))

"""
Run one (solver, method) pair under a wall-clock budget and return a plain Dict
describing what happened. Never throws: a solver failing is a result, not an error.
"""
function bench(solver::String, method::Symbol, budget::Float64)
    opt, attrs = solver_setup[solver]
    @info "benchmarking" solver method budget_s=round(budget, digits=1)

    t0 = time()
    payload = try
        res = solve_dopf(case, opt; method = method, attributes = attrs,
                         time_limit_sec = budget)
        Dict{String,Any}(
            "outcome"       => res.converged ? "completed" : "max_iter",
            "converged"     => res.converged,
            "solve_seconds" => res.solve_seconds,
            "n_iterations"  => length(res.iterations),
            "iterations"    => [Dict("iter" => r.iter, "seconds" => r.seconds,
                                     "objective" => r.objective, "residual" => r.residual,
                                     "status" => r.status) for r in res.iterations],
            "E_curt_kWh"    => kWh(case, sum(res.PVC)),
            "Vmin" => minimum(res.V), "Vmax" => maximum(res.V),
        )
    catch e
        msg = sprint(showerror, e)
        # solve_dopf's own message is "iteration N terminated with status S"
        it     = match(r"iteration (\d+)", msg)
        status = match(r"status (\w+)", msg)
        Dict{String,Any}(
            "outcome"        => occursin("TIME_LIMIT", msg) ? "withdrawn" : "failed",
            "converged"      => false,
            "failed_at_iter" => it     === nothing ? nothing : parse(Int, it.captures[1]),
            "status"         => status === nothing ? nothing : status.captures[1],
            "error"          => first(msg, 400),
        )
    end
    payload["wall_seconds"] = time() - t0
    payload["solver"]       = solver
    payload["method"]       = String(method)
    payload["budget_seconds"] = budget

    open(joinpath(OUTDIR, "$(lowercase(solver))_$(method).json"), "w") do io
        JSON3.pretty(io, payload)      # flush immediately: partial runs stay useful
    end
    return payload
end

results  = Dict{Tuple{String,Symbol},Any}()
baseline = Dict{Symbol,Float64}()

# Gurobi first — it sets the reference time for both encodings.
for method in (:bigm, :lambda)
    r = bench("Gurobi", method, 3600.0)
    results[("Gurobi", method)] = r
    baseline[method] = get(r, "solve_seconds", NaN)
    @info "  Gurobi reference" method seconds=round(baseline[method], digits = 1)
end

for solver in ("HiGHS", "GLPK"), method in (:bigm, :lambda)
    budget = WITHDRAW_FACTOR * baseline[method]
    results[(solver, method)] = bench(solver, method, budget)
end

# a compact index so the docs do not have to know the file naming
index = Dict{String,Any}(
    "withdraw_factor" => WITHDRAW_FACTOR,
    "reference"       => "Gurobi",
    "baseline_seconds" => Dict(String(m) => baseline[m] for m in (:bigm, :lambda)),
    "julia_version"   => string(VERSION),
    "runs" => [Dict("solver" => s, "method" => String(m),
                    "outcome" => results[(s,m)]["outcome"],
                    "wall_seconds" => results[(s,m)]["wall_seconds"],
                    "budget_seconds" => results[(s,m)]["budget_seconds"],
                    "n_iterations" => get(results[(s,m)], "n_iterations", 0),
                    "failed_at_iter" => get(results[(s,m)], "failed_at_iter", nothing),
                    "status" => get(results[(s,m)], "status", nothing),
                    "solve_seconds" => get(results[(s,m)], "solve_seconds", nothing))
               for s in ("Gurobi", "HiGHS", "GLPK"), m in (:bigm, :lambda)],
)
open(joinpath(OUTDIR, "index.json"), "w") do io
    JSON3.pretty(io, index)
end

println("\n", rpad("solver", 8), rpad("method", 9), rpad("outcome", 12),
        rpad("wall (s)", 11), rpad("budget (s)", 12), "detail")
for s in ("Gurobi", "HiGHS", "GLPK"), m in (:bigm, :lambda)
    r = results[(s,m)]
    detail = r["outcome"] in ("completed", "max_iter") ?
        "$(r["n_iterations"]) iters, curtail $(round(r["E_curt_kWh"], digits=1)) kWh" :
        "iter $(r["failed_at_iter"]) -> $(r["status"])"
    println(rpad(s, 8), rpad(String(m), 9), rpad(r["outcome"], 12),
            rpad(round(r["wall_seconds"], digits = 1), 11),
            rpad(round(r["budget_seconds"], digits = 1), 12), detail)
end
println("\nwrote ", OUTDIR)
