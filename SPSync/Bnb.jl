module Bnb

const N_Dumb, N_Root = 0, 1
mutable struct Stats
    tm::Float64 # time spent on solving LP masters
    tn::Float64 # time spent on solveing LP subproblems
    tsb::Float64 # time spent on MIP subproblems to generate SB cuts
    tev::Float64 # time spent on MIP subproblems to eval ub during training (to facilitate pruning)
    cb::Int # Number of times that a Ben cut cut off a Θs_check point
    csb::Int # NUmber of times that an SBen cut cut off a Θs_check point
    np::Int # Number of pruning events
end
function log_stats(NDˈS)
    (; tm, tn, tsb, tev, cb, csb, np) = NDˈS
    f = x -> round(x; digits = 1)
    tn, tsb, tev = f(tn), f(tsb), f(tev)
    tm = round(tm; digits = 5)
    println("lp_time=($tm, $tn), mip_time=($tsb, $tev); Ncut_off=($cb, $csb); Npruned=$np")
end
struct BnBNodeRef
    f::Int # From which node; Root node (node 1) is assumed to be from node 0
    # here we intentionally use Cdouble to facilitate GRB API call
    b::Tuple{Vector{Cdouble},Vector{Cdouble}} # (Xlb, Xub) [index 1 and 2]
    c::Vector{Vector{Float64}} # MyCut
end
nXfixed(r)=count(l==u for (l,u)=zip(r...))
struct NDT # a general info bundle
    # node classification
    nc::Dict{Int,Vector{Int}} # e.g. 0 => [1] ; 1 => [1, 2]; Node 1 and 2 has 1 fixed branching decision
    # Close: You shouldn't explore that node any more.
    o::Dict{Int, Bool} # Open: You can still explore that node.
    e::Dict{Int, Bool} # isleaf? for each node
    s::Dict{Int, Bool} # Do you intend to use sBen cut at this node?
    l::Dict{Int, Float64} # best lb at each node
    r::Dict{Int, BnBNodeRef} # reference for each node
    S::Stats # stats
    NDT(N) = new(
        Dict(0 => [N_Root]),
        Dict(N_Root => true),
        Dict(N_Root => true),
        Dict(N_Root => true),
        Dict(N_Root => Inf), # lb=Inf is deemed invalid
        Dict(N_Root => BnBNodeRef(N_Dumb, (fill(0.,N), fill(1.,N)), Vector{Float64}[])),
        Stats(0.,0.,0.,0.,0,0,0)
    )
end
# This function is used to label new nodes
currentMaxN(NDˈnc) = maximum(p -> maximum(p.second; init=N_Root), NDˈnc)
function global_lb(ND)
    lb, lD = Inf, ND.l
    for (k,v)=ND.e  v && (lb = min(lD[k], lb)) end
    lb == Inf ? error("global_lb == Inf") : return lb 
end

function select_i_X(Xres, NDˈrˈb)
    (lbv, ubv), (ib, Xb, dmin), Nf = NDˈrˈb, (0, 0., Inf), 0
    for ((i,x), l, u) = zip(enumerate(Xres), lbv, ubv)
        (l != u && isfrac(x)) || continue
        d = abs(x - .5)
        d < dmin && ((ib, Xb, dmin) = (i, x, d))
        Nf += 1
    end
    ib, Xb, Nf # if ib==0, Xres should be Int-Feasible
end
function branch_and_closeˈnonleaf_Nnow!(ND, Nnow, iBr, vBr, use_sB)
    b = round(Bool, vBr)
    N1 = new_a_node!(ND, Nnow, iBr,  b, use_sB) # small label to more promising branch
    N2 = new_a_node!(ND, Nnow, iBr, !b, use_sB) # the other possibility
    ND.o[Nnow] = ND.e[Nnow] = false
    N1, N2
end
function new_a_node!(ND, Nnow, iBr, b, use_sB)
    N = 1+currentMaxN(ND.nc) # Label a new node
    ND.r[N] = tRef = BnBNodeRef(Nnow,
        map(copy, ND.r[Nnow].b),
        map(copy, ND.r[Nnow].c)
    )
    tRef.b[2-b][iBr] = b # perform a new branching step
    ND.o[N] = ND.e[N] = true
    ND.l[N], ND.s[N] = ND.l[Nnow], use_sB & ND.s[Nnow]
    k = nXfixed(tRef.b)
    haskey(ND.nc, k) ? push!(ND.nc[k], N) : (ND.nc[k] = [N])
    N
end
function pick_an_open_node(ND, srchMod)
    f, (; nc, o) = ifelse(srchMod, reverse, identity), ND
    for l=f(0:maximum(keys(nc))), n=sort!(nc[l])
        o[n] && return n
    end;  return N_Dumb
end
function log_gaps(lb, ub; color = 52)
    agap = ub - lb
    rgap = agap / ub
    printstyled("lb= $lb, ub= $ub, agap= $agap, rgap= $rgap\n"; color = color)
    agap, rgap
end

isfrac(x) = ≈(x, 0.5; atol = 0.49999)
isIntFeas(chkˈX) = !any(isfrac, chkˈX)
end
