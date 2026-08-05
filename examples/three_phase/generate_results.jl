# Regenerate the precomputed three-phase results the documentation is built from.
#
#   julia --project=. examples/three_phase/generate_results.jl
#
# Each of the three shipped scripts is executed in its own module and its results are
# harvested, so the numbers in the documentation come from exactly the code a reader
# runs — there is no second copy of the model here to drift out of step.
#
# Writes docs/src/assets/results/threephase/{case,lambda,bigm,heaviside}.json.

using JSON3, Printf

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, "..", ".."))
const OUT  = joinpath(ROOT, "docs", "src", "assets", "results", "threephase")
mkpath(OUT)

"Run one of the standalone scripts in a private module and hand back its bindings."
function run_script(file)
    @info "running" file
    m = Module(Symbol("TP_", file))
    Core.eval(m, :(using Base))
    Base.include(m, joinpath(HERE, file))
    return m
end

# ---- the exact AC audit, applied to a solved module ------------------------------------
function audit(m)
    dq(v, qb) = v <= m.VBP[2] ? qb :
                v <= m.VBP[3] ? qb * (m.VBP[3] - v) / (m.VBP[3] - m.VBP[2]) :
                v <= m.VBP[4] ? 0.0 :
                v <= m.VBP[5] ? -qb * (v - m.VBP[4]) / (m.VBP[5] - m.VBP[4]) : -qb
    gap = 0.0; dev_true = 0.0; dev_model = 0.0; nviol = 0
    tmin = Inf; tmax = -Inf
    for t in 1:m.T
        Vt = m.sweep(t)
        gap  = max(gap, maximum(abs.(m.V[:, :, t] .- Vt)))
        tmin = min(tmin, minimum(Vt)); tmax = max(tmax, maximum(Vt))
        nviol += count(<(m.VLIM[1] - 1e-6), Vt) + count(>(m.VLIM[2] + 1e-6), Vt)
        for i in 1:m.npv
            b, φ, qb = m.PV[i].bus, m.PV[i].phase, m.PV[i].Smax
            dev_true  = max(dev_true,  abs(m.Qdg_v[i, t] - dq(Vt[b, φ], qb)))
            dev_model = max(dev_model, abs(m.Qdg_v[i, t] - dq(m.V[b, φ, t], qb)))
        end
    end
    return (; gap, dev_true, dev_model, nviol, true_lo = tmin, true_hi = tmax)
end

# ---- run all three ----------------------------------------------------------------------
mods = Dict("lambda"    => run_script("LinDist3Flow_Lambda.jl"),
            "bigm"      => run_script("LinDist3Flow_BigM.jl"),
            "heaviside" => run_script("LinDist3Flow_Heaviside.jl"))

# ---- case description, taken from the Lambda run ---------------------------------------
L = mods["lambda"]
kW(x) = x * L.SBASE / 1e3

loads_ph = [count(b -> L.Pload_pk[b, φ] > 0, 1:L.nb) for φ in 1:3]
loadkw_ph = [kW(sum(L.Pload_pk[:, φ])) for φ in 1:3]

case = Dict(
    "name" => L.CASE,
    "n_bus" => L.nb, "n_line" => L.nbr, "n_load" => length(L.net.load),
    "Vbase_V" => L.VBASE, "Sbase_kVA" => L.SBASE_KVA, "Zbase_ohm" => L.ZBASE,
    "length_m" => sum(Float64(l.length) for (_, l) in pairs(L.net.line)),
    "loads_per_phase" => loads_ph, "load_kW_per_phase" => loadkw_ph,
    "load_kW_total" => kW(sum(L.Pload_pk)),
    "Vmin_limit" => L.VLIM[1], "Vmax_limit" => L.VLIM[2],
    "Vbp" => L.VBP, "qshape" => L.QSHAPE,
    "n_steps" => L.T,
    "classes" => [Dict("name" => c[1], "P_kW" => c[2],
                       "S_kVA" => L.S_OVER_P * c[2],
                       "qbar_pu" => L.S_OVER_P * c[2] * 1e3 / L.SBASE)
                  for c in L.PV_CLASSES],
    "sites" => [Dict("bus" => L.BUSES[g.bus], "phase" => g.phase,
                     "class" => L.PV_CLASSES[g.cls][1],
                     "class_idx" => g.cls,
                     "Z_ohm" => L.DIST[L.BUSES[g.bus]] * L.ZBASE,
                     "P_kW" => kW(g.Pmax), "S_kVA" => kW(g.Smax))
                for g in L.PV],
    "PV_kW_total" => sum(kW(g.Pmax) for g in L.PV),
)
open(joinpath(OUT, "case.json"), "w") do io; JSON3.pretty(io, case); end

# ---- per-method results -----------------------------------------------------------------
for (tag, m) in mods
    a = audit(m)
    E_avail = m.kWh(sum(m.Pavail)); E_curt = m.kWh(sum(m.PVC_v))
    payload = Dict(
        "method" => m.METHOD,
        "solver" => tag == "heaviside" ? "Ipopt" : "Gurobi",
        "model_class" => tag == "heaviside" ? "NLP" : "MILP",
        "nvar" => m.num_variables(m.model),
        "nbin" => count(m.is_binary, m.all_variables(m.model)),
        "ncon" => m.num_constraints(m.model; count_variable_in_set_constraints = false),
        "solve_seconds" => m.secs, "status" => string(m.status),
        "E_avail_kWh" => E_avail, "E_curt_kWh" => E_curt,
        "curt_percent" => 100 * E_curt / E_avail,
        "Vmin" => minimum(m.V), "Vmax" => maximum(m.V),
        "Vmin_phase" => [minimum(m.V[:, φ, :]) for φ in 1:3],
        "Vmax_phase" => [maximum(m.V[:, φ, :]) for φ in 1:3],
        # per-time-step envelopes, for the daily figure
        "Vmax_t" => [[maximum(m.V[:, φ, t]) for t in 1:m.T] for φ in 1:3],
        "Vmin_t" => [[minimum(m.V[:, φ, t]) for t in 1:m.T] for φ in 1:3],
        # per-inverter series, for the droop-verification figure
        "Vdg_series" => [[m.V[m.PV[i].bus, m.PV[i].phase, t] for t in 1:m.T] for i in 1:m.npv],
        "Qdg_series" => [[m.Qdg_v[i, t] for t in 1:m.T] for i in 1:m.npv],
        "P_avail_kW" => [kW(sum(m.Pavail[:, t])) for t in 1:m.T],
        "P_disp_kW"  => [kW(sum(m.Pdg_v[:, t]))  for t in 1:m.T],
        # verification
        "max_droop_deviation" => a.dev_model,
        "audit" => Dict("v_gap" => a.gap, "droop_residual_true_v" => a.dev_true,
                        "n_limit_violations" => a.nviol,
                        "true_Vmin" => a.true_lo, "true_Vmax" => a.true_hi),
    )
    open(joinpath(OUT, "$tag.json"), "w") do io; JSON3.pretty(io, payload); end
    @printf("%-10s %-4s %6.1f s  curtail %6.2f kWh (%.3f %%)  V %.4f–%.4f  dev %.2e\n",
            tag, payload["model_class"], m.secs, E_curt, payload["curt_percent"],
            payload["Vmin"], payload["Vmax"], a.dev_model)
end

println("\nwrote ", OUT)
