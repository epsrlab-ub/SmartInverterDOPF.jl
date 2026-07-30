# Network, load and PV data for the tutorial case study.

"""
    Case

A radial distribution feeder together with a full day of load and irradiance data,
all in per unit on `Sbase`/`Vbase`, plus the smart-inverter fleet.

Time is indexed as `(h, m)` with `h ∈ HOUR_SET` (1:24) and `m ∈ QUARTER_SET` (1:4),
i.e. 96 quarter-hourly steps.
"""
struct Case
    # bases
    Sbase::Float64                     # VA
    Vbase::Float64                     # V
    Zbase::Float64                     # Ω
    # topology
    BUS_SET::UnitRange{Int}
    BRANCH_SET::Vector{Tuple{Int,Int}}
    Bi_BRANCH_SET::Vector{Tuple{Int,Int}}
    slack::Int
    R::Dict{Tuple{Int,Int},Float64}    # p.u.
    X::Dict{Tuple{Int,Int},Float64}    # p.u.
    # time
    HOUR_SET::UnitRange{Int}
    QUARTER_SET::UnitRange{Int}
    # demand, p.u., [bus, hour, quarter]
    Pload::Array{Float64,3}
    Qload::Array{Float64,3}
    # smart inverters
    DG_SET::Vector{Int}
    NON_DG_SET::Vector{Int}
    Pdg_max::Dict{Int,Float64}                 # nameplate active rating, p.u.
    Sdg_max::Dict{Int,Float64}                 # apparent rating, p.u.
    Pdg_max_vary::Dict{Int,Matrix{Float64}}    # irradiance-scaled ceiling, [hour, quarter]
    # operating limits
    Vmin::Float64
    Vmax::Float64
    Vnom::Float64
end

nbus(c::Case) = length(c.BUS_SET)
ndg(c::Case)  = length(c.DG_SET)

"""
    load_case(datadir = joinpath(pkgdir(SmartInverterDOPF), "data"); kwargs...)

Read `system_33bus.json`, `load_profiles_15min.json` and `solar_profile.json` from
`datadir` and assemble a [`Case`](@ref).

# Keyword arguments
- `dg_buses = [7, 18, 33]`: buses hosting a smart inverter.
- `dg_psize = [2.2e6, 1.0e6, 1.5e6]`: PV nameplate active power in W.
- `s_over_p = 1.1`: inverter apparent rating as a multiple of the active rating,
  i.e. the headroom available for reactive support.
- `vlimits = (0.95, 1.05)`: bus voltage limits in p.u.
"""
function load_case(datadir::AbstractString = joinpath(pkgdir(@__MODULE__), "data");
                   dg_buses = [7, 18, 33],
                   dg_psize = [2.2e6, 1.0e6, 1.5e6],
                   s_over_p = 1.1,
                   vlimits  = (0.95, 1.05))

    sys      = JSON3.read(read(joinpath(datadir, "system_33bus.json"), String))
    profiles = JSON3.read(read(joinpath(datadir, "load_profiles_15min.json"), String))
    solar    = JSON3.read(read(joinpath(datadir, "solar_profile.json"), String))

    Sbase = sys.base.S_base_kVA * 1e3
    Vbase = sys.base.V_base_kV  * 1e3
    Zbase = Vbase^2 / Sbase

    HOUR_SET, QUARTER_SET = 1:24, 1:4
    BUS_SET = 1:length(sys.buses)
    slack   = 1

    # peak demand per bus and its load class (1 industrial, 2 commercial, 3 residential;
    # class 0 marks the substation, which carries no load of its own)
    Pload_std = [b.P_load_kW   * 1e3 for b in sys.buses]
    Qload_std = [b.Q_load_kVAr * 1e3 for b in sys.buses]
    class     = Dict(b.bus => b.load_class for b in sys.buses)

    Pmult = hcat(profiles.load_percent.P.industrial,
                 profiles.load_percent.P.commercial,
                 profiles.load_percent.P.residential) ./ 100
    Qmult = hcat(profiles.load_percent.Q.industrial,
                 profiles.load_percent.Q.commercial,
                 profiles.load_percent.Q.residential) ./ 100

    Pload = zeros(length(BUS_SET), 24, 4)
    Qload = zeros(length(BUS_SET), 24, 4)
    for bus in BUS_SET[2:end], h in HOUR_SET, q in QUARTER_SET
        t = q + (h - 1) * length(QUARTER_SET)
        Pload[bus, h, q] = Pload_std[bus] * Pmult[t, class[bus]] / Sbase
        Qload[bus, h, q] = Qload_std[bus] * Qmult[t, class[bus]] / Sbase
    end

    R = Dict{Tuple{Int,Int},Float64}()
    X = Dict{Tuple{Int,Int},Float64}()
    for br in sys.branches
        R[(br.from, br.to)] = br.R_ohm / Zbase
        X[(br.from, br.to)] = br.X_ohm / Zbase
    end
    BRANCH_SET    = [(br.from, br.to) for br in sys.branches]
    Bi_BRANCH_SET = vcat(BRANCH_SET, [(b, a) for (a, b) in BRANCH_SET])

    # irradiance as a fraction of nameplate, reshaped to [hour, quarter]
    solar_profile = permutedims(reshape(collect(Float64, solar.G_percent),
                                        length(QUARTER_SET), length(HOUR_SET))) ./ 100

    DG_SET  = collect(Int, dg_buses)
    Pdg_max = Dict(DG_SET[i] => dg_psize[i] / Sbase for i in eachindex(DG_SET))
    Sdg_max = Dict(DG_SET[i] => s_over_p * dg_psize[i] / Sbase for i in eachindex(DG_SET))
    Pdg_max_vary = Dict(d => Pdg_max[d] .* solar_profile for d in DG_SET)
    NON_DG_SET   = setdiff(BUS_SET, DG_SET, [slack])

    return Case(Sbase, Vbase, Zbase,
                BUS_SET, BRANCH_SET, Bi_BRANCH_SET, slack, R, X,
                HOUR_SET, QUARTER_SET, Pload, Qload,
                DG_SET, NON_DG_SET, Pdg_max, Sdg_max, Pdg_max_vary,
                vlimits[1], vlimits[2], 1.0)
end

"""
    kWh(case, x_pu)

Convert a per-unit power summed over quarter-hourly steps into kWh.
"""
kWh(c::Case, x) = x * c.Sbase / 1e3 / 4
