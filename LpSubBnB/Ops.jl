"""
    provide solver agnostic APIs
    if change solver, modify low-level implementations in this file
"""
module Ops
import Gurobi, JuMP

struct BnBNodeRef
    f::Int # From which node; Root node (node 1) is assumed to be from node 0
    l::Vector{Cdouble} # here we intentionally use Cdouble to facilitate GRB API call
    u::Vector{Cdouble}
end
ni2gi4s(s) = #= GRB index in mst =# s-1
ni2gi4x(#= subscript of x_i =# i, S) = i+S-1
gi2ni4x(i, S) = i+1-S # inverse of the above

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

end

