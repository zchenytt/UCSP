
"""
Timeless quantities in this module
    ReserveCurve comes from 3% load + 5% wind (load depends on t, while wind don't depends on t, for simplicity, uses its AR(1):μ)
"""
module Case2383
import PowerModels, PGLib
const CaseName = "pglib_opf_case2383wp_k.m"

const WindZone = [1,2,3]
const WindN = [
    [494, 10, 16, 17, 18],
    [42, 43, 44, 30, 31],
    [105, 64, 1426, 63, 67]
]
const WindP = [ # PMax
    [2.34, 4.0, 7.2, 10.8, 25.2],
    [4.1, 4.1, 4.1, 8.0, 10.0],
    [4.3, 4.5, 4.95, 6.5, 7.5]
]

get_Case_Dict()=(PowerModels.silence();PowerModels.make_basic_network(PGLib.pglib(CaseName)))
_0!(F; a=9e-7)=for (i,e)=enumerate(F) if -a<e<a setindex!(F, 0., i) end end
get_PTDF(CaD)=(F=PowerModels.calc_basic_ptdf_matrix(CaD);_0!(F);F)

struct Bbyz
    v::Vector{Vector{Int}} # zone relations
    m::Matrix{Vector{Int}} # branch indicies, organized into zone pair
    P::Dict{Int, Float64} # explicitly extract PfMax for type stability
    Bbyz(N)=new([Int[] for _=1:N-1], Matrix{Vector{Int}}(undef,N,N), Dict{Int, Float64}())
end
function branch_by_zone(CaD; PfK=1.0)
    o = Bbyz(6)
    CaDˈbus, CaDˈbr = CaD["bus"], CaD["branch"]
    for i=1:length(CaDˈbr)
        CaDˈbrˈi = CaDˈbr[string(i)]
        f,t = CaDˈbrˈi["f_bus"],CaDˈbrˈi["t_bus"]
        f,t = CaDˈbus[string(f)]["zone"],CaDˈbus[string(t)]["zone"]
        f == t && continue
        # Now here are the lines that are inter-zone
        o.P[i] = PfK * CaDˈbrˈi["rate_a"]
        f,t = min(f,t), max(f,t)
        if t in o.v[f]
            push!(o.m[f,t], i)
        else
            o.m[f,t], _ = [i], push!(o.v[f], t)
        end
    end
    foreach(sort!, o.v)
    o
end
function get_load(rng, CaD)
    N, G, I = [Int[] for _=1:6], [Float64[] for _=1:6], [Int[] for _=1:6] # `I` is index of Curve
    CaDˈbus,CaDˈload = CaD["bus"],CaD["load"]
    Threshold = 0.001
    for g=1:length(CaDˈload)
        CaDˈloadˈg = CaDˈload[string(g)]
        p = CaDˈloadˈg["pd"]
        p < Threshold && continue
        n = CaDˈloadˈg["load_bus"]
        z = CaDˈbus[string(n)]["zone"]
        curve_index = rand(rng, 1:6) # It happens that #Zone == #LoadCurve == 6
        _,_,_ = push!(N[z], n),push!(G[z], p),push!(I[z], curve_index)
    end
    (n=N, g=G, i=I)
end

get_type(PMa) = if PMa < 0.5 # used by `case2383`
        i = 6 # gasCT (Smallest)
    elseif PMa < 0.505
        i = 5 # oil
    elseif PMa < 0.9
        i = 4 # hydro
    elseif PMa < 1.4
        i = 3 # gasCC
    elseif PMa < 10.
        i = 2 # coal
    else
        i = 1 # nuclear (Biggest)
end

"""
Hold-time generator, same for Dn/Up Times
"""
_HT(rng, Ty) = if Ty == 1 # nuclear
        23
    elseif Ty == 2 # coal
        rand(rng, 12:23)
    elseif Ty == 3 # gasCC
        rand(rng, 6:23)
    elseif Ty == 4 # hydro
        rand(rng, 1:6)
    elseif Ty == 5 # oil
        rand(rng, 3:12)
    elseif Ty == 6 # gasCT
        rand(rng, 1:4)
end

"""
(Dn, UP)
This is Int, not readily (RD, RU). They depends on GMax, GMin
"""
_RampPeriods(rng, Ty) = if Ty == 1 # nuclear
        23,23
    elseif Ty == 2 # coal
        rand(rng, 6:17),rand(rng, 8:23)
    elseif Ty == 3 # gasCC
        rand(rng, 2:8),rand(rng, 3:10)
    elseif Ty == 4 # hydro
        2,2
    elseif Ty == 5 # oil
        rand(rng, 2:6), rand(rng, 2:6)
    elseif Ty == 6 # gasCT
        2,2
end

"""
relative_cost::Int
"""
_SUcostOprcostRatio(rng, Ty) = if Ty == 1 # nuclear
        rand(rng, 50:150)
    elseif Ty == 2 # coal
        rand(rng, 20:80)
    elseif Ty == 3 # gasCC
        rand(rng, 10:40)
    elseif Ty == 4 # hydro
        rand(rng, 1:5)
    elseif Ty == 5 # oil
        rand(rng, 10:50)
    elseif Ty == 6 # gasCT
        rand(rng, 2:15)
end
_SDcostOprcostRatio(rng, Ty) = if Ty == 1 # nuclear
        rand(rng, 5:20)      # procedural + safety-driven shutdowns, rare but costly
    elseif Ty == 2 # coal
        rand(rng, 2:10)      # thermal stress, ramp-down inefficiency
    elseif Ty == 3 # gasCC
        rand(rng, 1:6)       # relatively flexible, moderate shutdown cost
    elseif Ty == 4 # hydro
        rand(rng, 0:2)       # essentially negligible operational shutdown cost
    elseif Ty == 5 # oil
        rand(rng, 2:8)       # similar to gas but slightly higher handling cost
    elseif Ty == 6 # gasCT
        rand(rng, 1:5)       # fast-start peaker, low shutdown penalty
end

"""
SD/SU rate is calculated from RD/RU as clamp(R2S * Rrate, Gmin, GMax)
"""
_r2s(rng) = rand(rng, .8:.1:1.5)
_Ci() = [Int[] for _=1:6]
_Cf() = [Float64[] for _=1:6]
function get_gen(rng, CaD)
    N, RUI, RDI, HUT, HDT, SCI, DCI = _Ci(), _Ci(), _Ci(), _Ci(), _Ci(), _Ci(), _Ci()
    G, MmRatio, Gmin, Clin, R2S = _Cf(), _Cf(), _Cf(), _Cf(), _Cf()
    CaDˈbus,CaDˈgen = CaD["bus"],CaD["gen"]
    Threshold = 0.001 * sum(v["pmax"] for (_,v)=CaDˈgen)
    for g=1:length(CaDˈgen)
        CaDˈgenˈg = CaDˈgen[string(g)]
        p = CaDˈgenˈg["pmax"]
        p < Threshold && continue
        n = CaDˈgenˈg["gen_bus"]
        z = CaDˈbus[string(n)]["zone"]
        push!(N[z], n) # node
        push!(G[z], p) # PMax
        Ty = get_type(p) # Type TODO need to export this info?
        pm = CaDˈgenˈg["pmin"]; push!(Gmin[z], pm) # Pmin
        ra = p / pm; push!(MmRatio[z], ra) # PMax/Pmin-Ratio 
        push!(HUT[z], _HT(rng, Ty)); push!(HDT[z], _HT(rng, Ty)) # minimun Dn/Up time
        C = CaDˈgenˈg["cost"][2]/5000; C < 0.2 && (C=rand(rng, 0.3:1e-4:0.5)); push!(Clin[z], C) # LinCost
        push!(SCI[z], _SUcostOprcostRatio(rng, Ty)) # Start-Up Cost (indirect)
        push!(DCI[z], _SDcostOprcostRatio(rng, Ty)) # Shut-Dw Cost (indirect)
        ø1,ø2 = _RampPeriods(rng, Ty); push!(RDI[z], ø1); push!(RUI[z], ø2) # normal Ramp (indirect)
        push!(R2S[z], _r2s(rng)) # Start/Shut Ramp (indirect)
    end
    (
        n = N,
        pmin = Gmin,
        pmax = G,
        MmRatio = MmRatio, # if <1.1, don't add ramp constrs
        UT = HUT, # ::Int
        DT = HDT,
        C = Clin,
        sci = SCI, # start-up cost::Int
        dci = DCI, 
        rdi = RDI, # ramp-down::Int
        rui = RUI,
        r2s = R2S # a ratio::Float64
    )
end

const CaD = get_Case_Dict()
const F = get_PTDF(CaD)

end
