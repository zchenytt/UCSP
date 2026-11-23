"""
implement the sync algorithm
"""
module SyncAlgo
import ..Ops, JuMP, ..Sin_Tha_Models, ..Sin_Tha_Cgen

function explore_node!(ND, Nnow, ubs, Xres, nnt, new_m, li0, zi0, S, xH, pAH, GD, LD, liD, ato)
    (lbv, ubv), vCut, sTr, NDˈS = ND.r[Nnow].b, ND.r[Nnow].c, ND.s[Nnow], ND.S
    Needsˈinit, N = isempty(vCut), length(lbv)
    m = Sin_Tha_Models.mst!(new_m, S, xH, pAH, GD; θ_has_obj= !Needsˈinit)
    Ops.set_lb!(m.m.moi_backend, 1, N, lbv);Ops.set_ub!(m.m.moi_backend, 1, N, ubv)
    for c = vCut # pre-build CTPLN model
        xBˈCn, n = c[N+1], nnt.t[1][1]; (; Ci, Cd) = n
        for (Pi, l, u, i) = zip(c, lbv, ubv, eachindex(Cd))
            iF = l == u
            Cd[i] = (!iF)Pi # Fixed term has no slope
            xBˈCn += iF * u * Pi # Fixed terms Fold into Cn
        end; Cd[N+1] = -1
        Ops.add_le_constr!(m.m.moi_backend, N+1, Ci, Cd, -xBˈCn) # 🟩
    end
    if Needsˈinit
        JuMP.optimize!(m.m); Ops.assert_optimal(m.m, "mst_ini"); Sin_Tha_Cgen.load_mst_sol!(m)
        Sin_Tha_Cgen.ini_cut_gen!(vCut, lbv, ubv, m, nnt.t[1], ato)
        Ops.set_obj_coeff!(m.m.moi_backend, 0, 1, [1.]) # 🟦
    end
    NchkCutOff = 0
    for cut_gen = 1:typemax(Int)
        JuMP.optimize!(m.m); Ops.assert_optimal(m.m, "cutgen_mst"); Sin_Tha_Cgen.load_mst_sol!(m)
        NDˈS.tm += JuMP.solve_time(m.m)
        m.chk.Θ[1] > ubs && (ND.o[Nnow] = false; NDˈS.np += 1; break)
        NchkCutOff += Added = Sin_Tha_Cgen.Ben_attempt!(vCut, lbv, ubv, m, nnt.t[1], ato, NDˈS)
        Added && (NDˈS.cb += 1; continue)
        if sTr
            NchkCutOff += Added = Sin_Tha_Cgen.sB_attempt!(vCut, lbv, ubv, m, nnt.t, ato, NDˈS)
            Added && (NDˈS.csb += 1; continue)
        end
        Xres .= m.chk.X # load solution info
        break # no more need to explore this node
    end
    ND.l[Nnow] = m.chk.Θ[1] # load lb info
    ND.o[Nnow] || return NchkCutOff
    Mu = similar(vCut, Cdouble)
    Ops.load_F64_vec_Pi!(m.m.moi_backend, 0, length(Mu), Mu)
    deleteat!(vCut, Mu .== 0)
    return NchkCutOff
end
function get_valid_ub(X, nˈe, S, ato, NDˈS; rgap = 0.01)
    ato[] = 0.
    aTo = Threads.Atomic{Float64}(0)
    Threads.@threads for n = nˈe
        Ops.set_lb!(n.m.moi_backend, 0, length(X), X)
        Ops.set_ub!(n.m.moi_backend, 0, length(X), X)
        JuMP.set_attribute(n.m, "MIPGap", rgap)
        JuMP.optimize!(n.m); Ops.assert_optimal(n.m, "eval_mip")
        Threads.atomic_max!(aTo, JuMP.solve_time(n.m))
        obv = JuMP.objective_value(n.m)
        Threads.atomic_add!(ato, obv/S)
    end
    NDˈS.tev += aTo[]
    ato[]
end

end

