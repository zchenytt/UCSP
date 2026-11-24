# 21: deep blue
# 27: blue
# 52: black red
# 36: green
# 93: purple
# 208: brown
"""
    provide solver agnostic APIs
    if change solver, modify low-level implementations in this file
"""
module Ops
import JuMP, Gurobi

const GRB_CONFIG = Dict{String,Any}("Threads"=>1,"OutputFlag"=>0)
# Low Level
    load_F64_vec!(o, start, len, v)    = Gurobi.GRBgetdblattrarray(o, "X", start, len, v)
    load_F64_vec_RC!(o, start, len, v) = Gurobi.GRBgetdblattrarray(o, "RC", start, len, v)
    load_F64_vec_Pi!(o, start, len, v) = Gurobi.GRBgetdblattrarray(o, "Pi", start, len, v)
    set_lb!(o, start, len, v)          = Gurobi.GRBsetdblattrarray(o, "LB", start, len, v)
    set_ub!(o, start, len, v)          = Gurobi.GRBsetdblattrarray(o, "UB", start, len, v)
    set_obj_coeff!(o, start, len, v)   = Gurobi.GRBsetdblattrarray(o, "Obj", start, len, v)
    add_le_constr!(o, read_len, ci, cd, Rhs) = Gurobi.GRBaddconstr(o, read_len, ci, cd, Cchar('<'), Rhs, C_NULL)

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

function assert_optimal(model, model_name)
    t = JuMP.termination_status(model)
    t == JuMP.OPTIMAL && return
    s = string(model_name, ": ", t)
    @error(s); error(s)
end
function assert_nodelimit(model, model_name)
    t = JuMP.termination_status(model)
    (t == JuMP.OPTIMAL || t == JuMP.NODE_LIMIT) && return
    s = string(model_name, ": ", t)
    @error(s); error(s)
end
function printinfo()
    a,b,c,d = Threads.nthreads(:default),Threads.nthreads(:interactive),Threads.threadid(),Threads.threadpool()
    println("Settings> Threads=($a,$b), tid=$c, pool=$d")
end

end
