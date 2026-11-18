"""
    provide solver agnostic APIs
    if change solver, modify low-level implementations in this file
"""
module Ops
import JuMP, Gurobi
const GRB_CONFIG = Dict{String,Any}("Threads"=>1,"OutputFlag"=>0)
const N_Dumb, N_Root = 0, 1

ni2gi4s(s) = #= GRB index in mst =# s-1
ni2gi4x(#= subscript of x_i =# i, S) = i+S-1
gi2ni4x(i, S) = i+1-S # inverse of the above
struct BnBNodeRef
    f::Int # From which node; Root node (node 1) is assumed to be from node 0
    # here we intentionally use Cdouble to facilitate GRB API call
    b::Tuple{Vector{Cdouble},Vector{Cdouble}} # (Xlb, Xub) [index 1 and 2]
end
nXfixed(r)=count(l==u for (l,u)=zip(r...))
struct NDT # a general info bundle
    # node classification
    nc::Dict{Int,Vector{Int}} # e.g. 0 => [1] ; 1 => [1, 2]; Node 1 and 2 has 1 fixed branching decision
    o::Dict{Int, Bool} # isopen? (to guide exploration) for each node
    e::Dict{Int, Bool} # isleaf? for each node
    l::Dict{Int, Float64} # best lb at each node
    r::Dict{Int, BnBNodeRef} # reference for each node
    NDT(N, MaxNxFixed) = new(
        Dict(i => i==0 ? [N_Root] : Int[] for i=0:MaxNxFixed),
        Dict(N_Root => true),
        Dict(N_Root => true),
        Dict(N_Root => Inf), # lb=Inf is deemed invalid
        Dict(N_Root => BnBNodeRef(N_Dumb, (fill(0.,N), fill(1.,N))))
    )
end

# This function is used to label new nodes
currentMaxN(NDˈnc) = maximum(p -> maximum(p.second; init=1), NDˈnc)
function global_lb(ND)
    lb, lD = Inf, ND.l
    for (k,v)=ND.e  v && (lb = min(lD[k], lb)) end
    lb == Inf ? error("lb == Inf") : return lb 
end

# high level
    function load_mst_sol!(m)
        (; N, S, chk) = m; (; X, Θ, ul), o = chk, m.m.moi_backend
        load_F64_vec!(o, 0, S, Θ)
        load_F64_vec!(o, S, N, X)
        ul[S+1] = JuMP.objective_bound(m.m)
    end
    recover_θ_obj_coeff!(m)=set_obj_coeff!(m.m.moi_backend, 0, m.S, fill(1/m.S, m.S))

# Utility
    isfrac(x) = ≈(x, 0.5; atol = 0.49999)
    isIntFeas(chkˈX) = !any(isfrac, chkˈX)
    function assert_optimal(model, model_name)
        t = JuMP.termination_status(model)
        t == JuMP.OPTIMAL && return
        s = string(model_name, ": ", t)
        @error(s); error(s)
    end

# Low Level
    load_F64_vec!(o, start, len, v)    = Gurobi.GRBgetdblattrarray(o, "X", start, len, v)
    load_F64_vec_RC!(o, start, len, v) = Gurobi.GRBgetdblattrarray(o, "RC", start, len, v)
    set_lb!(o, start, len, v)          = Gurobi.GRBsetdblattrarray(o, "LB", start, len, v)
    set_ub!(o, start, len, v)          = Gurobi.GRBsetdblattrarray(o, "UB", start, len, v)
    add_le_constr!(o, read_len, ci, cd, Rhs) = Gurobi.GRBaddconstr(o, read_len, ci, cd, Cchar('<'), Rhs, C_NULL)
    set_obj_coeff!(o, start, len, v)   = Gurobi.GRBsetdblattrarray(o, "Obj", start, len, v)

# GRB related
    Env() = Gurobi.Env(GRB_CONFIG)
    Model() = Model(Env())
    function Env(N::Int)
        v = Vector{Gurobi.Env}(undef, N)
        Threads.@threads for i=eachindex(v)
            v[i] = Env()
        end
        v
    end
    function Model(e::Gurobi.Env)
        m = JuMP.direct_model(Gurobi.Optimizer(e))
        JuMP.set_string_names_on_creation(m, false)
        m
    end
    function Model!(mv::Vector{JuMP.Model}, i, #=with existing ones=# ev::Vector{Gurobi.Env})
        m = JuMP.direct_model(Gurobi.Optimizer(ev[i]))
        JuMP.set_string_names_on_creation(m, false)
        mv[i] = m
    end
    Model!(mv::Vector{JuMP.Model}, #=with existing ones=# ev::Vector{Gurobi.Env}) = Threads.@threads for i=eachindex(mv)
        Model!(mv, i, ev)
    end

function printinfo()
    a,b,c,d = Threads.nthreads(:default),Threads.nthreads(:interactive),Threads.threadid(),Threads.threadpool()
    println("Settings> Threads=($a,$b), tid=$c, pool=$d")
end
end


# 21: deep blue
# 27: blue
# 36: green
# 208: brown
