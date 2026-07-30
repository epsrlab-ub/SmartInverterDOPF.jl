# Current-voltage AC optimal power flow (IVACOPF), solved by successive linearisation,
# with the Volt-VAr droop plugged in as a self-contained module.

"""
    DOPFResult

Outcome of [`solve_dopf`](@ref).

| field | meaning |
|:--|:--|
| `V` | bus voltage magnitude, `[bus, hour, quarter]`, p.u. |
| `Vdg`, `Qdg`, `Pdg`, `PVC` | per-inverter voltage, reactive output, active output and curtailment, `[dg, hour, quarter]`, p.u. |
| `Ploss` | total network active loss over the day, p.u. |
| `iterations` | per-iteration solve time, objective and linearisation residual |
| `nvar`, `nbin`, `ncon` | size of the model actually handed to the solver |
| `converged` | whether the residual fell below `tol` |
| `solve_seconds` | summed solver time over all iterations |
"""
struct DOPFResult
    V::Array{Float64,3}
    Vdg::Array{Float64,3}
    Qdg::Array{Float64,3}
    Pdg::Array{Float64,3}
    PVC::Array{Float64,3}
    Ploss::Float64
    iterations::Vector{NamedTuple{(:iter, :seconds, :objective, :residual, :status),
                                  Tuple{Int,Float64,Float64,Float64,String}}}
    nvar::Int
    nbin::Int
    ncon::Int
    converged::Bool
    solve_seconds::Float64
end

"""
    solve_dopf(case, optimizer; method = :lambda, curve = ieee1547_curve(),
               max_iter = 15, tol = 1e-6, silent = true)

Minimise total PV curtailment over the day subject to the IVACOPF network model and
the Volt-VAr droop encoded by `method` (`:bigm`, `:lambda` or `:heaviside`).

The AC power flow enters through two bilinear identities — the `v·I` power balance and
the `|I|²` loss — which are linearised about the previous iterate and refreshed until
the residual of the exact loss identity falls below `tol`. `optimizer` must therefore
match the droop encoding: an MILP solver for `:bigm` and `:lambda`, an NLP solver for
`:heaviside`.

```julia
using SmartInverterDOPF, Gurobi
res = solve_dopf(load_case(), Gurobi.Optimizer; method = :lambda)
```
"""
function solve_dopf(c::Case, optimizer;
                    method::Symbol = :lambda,
                    curve::DroopCurve = ieee1547_curve(),
                    max_iter::Int = 15,
                    tol::Float64 = 1e-6,
                    silent::Bool = true,
                    time_limit_sec = 3600,
                    attributes = Dict{String,Any}())

    HOUR_SET, QUARTER_SET = c.HOUR_SET, c.QUARTER_SET
    BUS_SET, BRANCH_SET   = c.BUS_SET, c.BRANCH_SET
    Bi_BRANCH_SET         = c.Bi_BRANCH_SET
    DG_SET, NON_DG_SET    = c.DG_SET, c.NON_DG_SET
    SLACK_SET             = [c.slack]
    R, X                  = c.R, c.X
    nb                    = length(BUS_SET)
    qbar                  = Dict(d => c.Sdg_max[d] for d in DG_SET)   # reactive capability

    # linearisation point, initialised at a flat start
    v_r_pr    = ones(nb, 24, 4)
    v_im_pr   = zeros(nb, 24, 4)
    Ibs_r_pr  = zeros(nb, 24, 4)
    Ibs_im_pr = zeros(nb, 24, 4)
    Ibr_r_pr  = Dict((br, h, q) => 0.0 for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET)
    Ibr_im_pr = Dict((br, h, q) => 0.0 for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET)

    log = NamedTuple{(:iter, :seconds, :objective, :residual, :status),
                     Tuple{Int,Float64,Float64,Float64,String}}[]
    result = nothing
    nvar = nbin = ncon = 0
    converged = false
    total_solve = 0.0

    for it in 1:max_iter
        model = Model(optimizer)
        silent && set_silent(model)
        for (k, val) in attributes
            set_optimizer_attribute(model, k, val)
        end
        time_limit_sec === nothing || set_time_limit_sec(model, time_limit_sec)

        # ---- host variables ---------------------------------------------------------
        @variable(model, 0 <= Pgen[SLACK_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Qgen[SLACK_SET, HOUR_SET, QUARTER_SET])
        @variable(model, c.Vmin <= v[BUS_SET, HOUR_SET, QUARTER_SET] <= c.Vmax, start = c.Vnom)
        @variable(model, v_r[BUS_SET, HOUR_SET, QUARTER_SET], start = c.Vnom)
        @variable(model, v_im[BUS_SET, HOUR_SET, QUARTER_SET], start = 0)
        @variable(model, Ibr_r[BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Ibr_im[BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Psnd[Bi_BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Qsnd[Bi_BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, 0 <= Ploss[BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, 0 <= Qloss[BRANCH_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Ibs_r[BUS_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Ibs_im[BUS_SET, HOUR_SET, QUARTER_SET])
        @variable(model, 0 <= Pdg[DG_SET, HOUR_SET, QUARTER_SET])
        @variable(model, Qdg[DG_SET, HOUR_SET, QUARTER_SET], start = 0)
        @variable(model, 0 <= PVC[DG_SET, HOUR_SET, QUARTER_SET])

        # ---- the droop module -------------------------------------------------------
        add_droop!(model, method, curve, v, Qdg, DG_SET, HOUR_SET, QUARTER_SET, qbar)

        # ---- inverter capability: a 16-segment polygon inscribing the S-circle ------
        k = 16
        for l in 1:k
            θ = l * π / k
            @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                cos(θ) * Pdg[d,h,m] + sin(θ) * Qdg[d,h,m] <=  c.Sdg_max[d])
            @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                cos(θ) * Pdg[d,h,m] + sin(θ) * Qdg[d,h,m] >= -c.Sdg_max[d])
        end
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Pdg[d,h,m] <= c.Pdg_max_vary[d][h,m])
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Pdg[d,h,m] <= c.Pdg_max[d])
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Qdg[d,h,m] <= c.Sdg_max[d])

        # ---- slack reference --------------------------------------------------------
        @constraint(model, [i in SLACK_SET, h in HOUR_SET, m in QUARTER_SET], v_r[i,h,m]  == c.Vnom)
        @constraint(model, [i in SLACK_SET, h in HOUR_SET, m in QUARTER_SET], v_im[i,h,m] == 0)

        # ---- Ohm's law along each branch (exact, linear) ----------------------------
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            v_r[i,h,m] - v_r[j,h,m] == R[(i,j)]*Ibr_r[(i,j),h,m] - X[(i,j)]*Ibr_im[(i,j),h,m])
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            v_im[i,h,m] - v_im[j,h,m] == R[(i,j)]*Ibr_im[(i,j),h,m] + X[(i,j)]*Ibr_r[(i,j),h,m])

        # ---- KCL at every bus (exact, linear) ---------------------------------------
        @constraint(model, [bus in BUS_SET, h in HOUR_SET, m in QUARTER_SET],
            Ibs_r[bus,h,m] == sum(Ibr_r[(bus,j),h,m] for (i,j) in BRANCH_SET if i == bus)
                            - sum(Ibr_r[(i,bus),h,m] for (i,j) in BRANCH_SET if j == bus))
        @constraint(model, [bus in BUS_SET, h in HOUR_SET, m in QUARTER_SET],
            Ibs_im[bus,h,m] == sum(Ibr_im[(bus,j),h,m] for (i,j) in BRANCH_SET if i == bus)
                             - sum(Ibr_im[(i,bus),h,m] for (i,j) in BRANCH_SET if j == bus))

        # ---- power balance: v·I linearised about the previous iterate ---------------
        Plin(i,h,m) = v_r_pr[i,h,m]*Ibs_r[i,h,m]  + v_im_pr[i,h,m]*Ibs_im[i,h,m] +
                      Ibs_r_pr[i,h,m]*v_r[i,h,m]  + Ibs_im_pr[i,h,m]*v_im[i,h,m] -
                      v_r_pr[i,h,m]*Ibs_r_pr[i,h,m] - v_im_pr[i,h,m]*Ibs_im_pr[i,h,m]
        Qlin(i,h,m) = v_im_pr[i,h,m]*Ibs_r[i,h,m] - v_r_pr[i,h,m]*Ibs_im[i,h,m] +
                      Ibs_r_pr[i,h,m]*v_im[i,h,m] - Ibs_im_pr[i,h,m]*v_r[i,h,m] -
                      v_im_pr[i,h,m]*Ibs_r_pr[i,h,m] + v_r_pr[i,h,m]*Ibs_im_pr[i,h,m]

        @constraint(model, [i in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
            Pgen[i,h,m] == Plin(i,h,m))
        @constraint(model, [i in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
            Qgen[i,h,m] == Qlin(i,h,m))
        @constraint(model, [i in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
            -c.Pload[i,h,m] == Plin(i,h,m))
        @constraint(model, [i in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
            -c.Qload[i,h,m] == Qlin(i,h,m))
        @constraint(model, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            -c.Pload[i,h,m] + Pdg[i,h,m] == Plin(i,h,m))
        @constraint(model, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            -c.Qload[i,h,m] + Qdg[i,h,m] == Qlin(i,h,m))

        # ---- branch loss: |I|² linearised about the previous iterate ----------------
        Isq(i,j,h,m) = 2*Ibr_r_pr[((i,j),h,m)]*Ibr_r[(i,j),h,m]   - Ibr_r_pr[((i,j),h,m)]^2 +
                       2*Ibr_im_pr[((i,j),h,m)]*Ibr_im[(i,j),h,m] - Ibr_im_pr[((i,j),h,m)]^2
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            Ploss[(i,j),h,m] == R[(i,j)] * Isq(i,j,h,m))
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            Qloss[(i,j),h,m] == X[(i,j)] * Isq(i,j,h,m))

        # ---- sending-end power bookkeeping ------------------------------------------
        @constraint(model, [bus in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
            Pgen[bus,h,m] == sum(Psnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Pload[bus,h,m])
        @constraint(model, [bus in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Pdg[bus,h,m] == sum(Psnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Pload[bus,h,m])
        @constraint(model, [bus in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
            0 == sum(Psnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Pload[bus,h,m])
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            Psnd[(i,j),h,m] == Ploss[(i,j),h,m] - Psnd[(j,i),h,m])
        @constraint(model, [bus in SLACK_SET, h in HOUR_SET, m in QUARTER_SET],
            Qgen[bus,h,m] == sum(Qsnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Qload[bus,h,m])
        @constraint(model, [bus in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Qdg[bus,h,m] == sum(Qsnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Qload[bus,h,m])
        @constraint(model, [bus in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET],
            0 == sum(Qsnd[(i,j),h,m] for (i,j) in Bi_BRANCH_SET if i == bus) + c.Qload[bus,h,m])
        @constraint(model, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
            Qsnd[(i,j),h,m] == Qloss[(i,j),h,m] - Qsnd[(j,i),h,m])

        # ---- voltage magnitude, linearised about the previous iterate ---------------
        @constraint(model, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET],
            v[i,h,m] == (v_r_pr[i,h,m]/sqrt(v_r_pr[i,h,m]^2 + v_im_pr[i,h,m]^2))*v_r[i,h,m]
                      + (v_im_pr[i,h,m]/sqrt(v_r_pr[i,h,m]^2 + v_im_pr[i,h,m]^2))*v_im[i,h,m])

        # ---- curtailment and objective ----------------------------------------------
        @constraint(model, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            PVC[d,h,m] == c.Pdg_max_vary[d][h,m] - Pdg[d,h,m])
        @objective(model, Min, sum(PVC[d,h,m] for d in DG_SET, h in HOUR_SET, m in QUARTER_SET))

        nvar = num_variables(model)
        nbin = count(is_binary, all_variables(model))
        ncon = num_constraints(model; count_variable_in_set_constraints = false)

        secs = @elapsed optimize!(model)
        total_solve += secs
        status = termination_status(model)
        # ALMOST_LOCALLY_SOLVED is tolerated, and recorded: the Heaviside encoding is
        # non-smooth at every breakpoint, so an interior-point solver routinely stops
        # just short of its convergence tolerance. That is the price of the encoding,
        # not a failure of the model, and the tutorial reports it rather than hiding it.
        status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ||
            error("iteration $it terminated with status $status")

        # residual of the exact (bilinear) loss identity: how far the linearisation is
        # from the true power flow
        residual = maximum(abs(value((v_r[i,h,m] - v_r[j,h,m]) * Ibr_r[(i,j),h,m]
                                   + (v_im[i,h,m] - v_im[j,h,m]) * Ibr_im[(i,j),h,m]
                                   - Ploss[(i,j),h,m]))
                           for (i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET)

        push!(log, (iter = it, seconds = secs, objective = objective_value(model),
                    residual = residual, status = string(status)))

        # refresh the linearisation point
        for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET
            v_r_pr[b,h,q]    = value(v_r[b,h,q])
            v_im_pr[b,h,q]   = value(v_im[b,h,q])
            Ibs_r_pr[b,h,q]  = value(Ibs_r[b,h,q])
            Ibs_im_pr[b,h,q] = value(Ibs_im[b,h,q])
        end
        for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET
            Ibr_r_pr[(br,h,q)]  = value(Ibr_r[br,h,q])
            Ibr_im_pr[(br,h,q)] = value(Ibr_im[br,h,q])
        end

        result = (V   = [value(v[i,h,m])    for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET],
                  Vdg = [value(v[d,h,m])    for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
                  Qdg = [value(Qdg[d,h,m])  for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
                  Pdg = [value(Pdg[d,h,m])  for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
                  PVC = [value(PVC[d,h,m])  for d in DG_SET,  h in HOUR_SET, m in QUARTER_SET],
                  Ploss = sum(value(Ploss[(i,j),h,m])
                              for (i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET))

        if residual < tol
            converged = true
            break
        end
    end

    return DOPFResult(result.V, result.Vdg, result.Qdg, result.Pdg, result.PVC,
                      result.Ploss, log, nvar, nbin, ncon, converged, total_solve)
end

"""
    base_case_voltages(case)

Voltage magnitudes for the same feeder and demand with **no** smart inverters, obtained
by a backward/forward sweep power flow. This is the reference the droop-aware dispatch
is compared against. The branch list of this feeder is already in topological order.
"""
function base_case_voltages(c::Case)
    nb, nbr = length(c.BUS_SET), length(c.BRANCH_SET)
    V = zeros(nb, 24, 4)
    for h in c.HOUR_SET, q in c.QUARTER_SET
        Vc  = ones(ComplexF64, nb)
        Ibr = zeros(ComplexF64, nbr)
        for _ in 1:30
            Ibus = [conj((c.Pload[b,h,q] + im*c.Qload[b,h,q]) / Vc[b]) for b in c.BUS_SET]
            for k in nbr:-1:1                                    # backward: currents
                (_, j) = c.BRANCH_SET[k]
                Ibr[k] = Ibus[j] + sum(Ibr[t] for t in 1:nbr
                                       if c.BRANCH_SET[t][1] == j; init = 0.0 + 0.0im)
            end
            for k in 1:nbr                                       # forward: voltages
                (i, j) = c.BRANCH_SET[k]
                Vc[j] = Vc[i] - (c.R[(i,j)] + im*c.X[(i,j)]) * Ibr[k]
            end
        end
        V[:,h,q] = abs.(Vc)
    end
    return V
end
