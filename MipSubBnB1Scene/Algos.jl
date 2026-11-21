module Algos
import ..QsMst, ..Ops, JuMP

function get_valid_ub(NDˈS, X, n; rgap = 0.01)
    Ops.set_lb!(n.m.moi_backend, 0, length(X), X)
    Ops.set_ub!(n.m.moi_backend, 0, length(X), X)
    JuMP.set_attribute(n.m, "MIPGap", rgap)
    print("Eval: "); JuMP.optimize!(n.m)
    NDˈS.tev += t = JuMP.solve_time(n.m)
    println("Finished in $(round(t))s, rgap=$rgap.")
    Ops.assert_optimal(n.m, "eval_mip")
    JuMP.objective_value(n.m)
end
function explore_node!(ND, Nnow, ubs, Xres, ng, new_m, li0, zi0, S, xH, pAH, GD, LD, liD)
    (lbv, ubv), vCut, sTr, NDˈS = ND.r[Nnow].b, ND.r[Nnow].c, ND.s[Nnow], ND.S
    Needsˈinit, N = isempty(vCut), length(lbv)
    m = QsMst.mst!(new_m, li0, zi0, S, xH, pAH, GD, LD, liD; θ_has_obj= !Needsˈinit)
    Ops.set_lb!(m.m.moi_backend, S, N, lbv); Ops.set_ub!(m.m.moi_backend, S, N, ubv)
    for c = vCut # pre-build CTPLN model
        xBˈCn, n = c[N+2], ng.t[1]; Cd = n.Cd
        for (Pi, l, u, i) = zip(c, lbv, ubv, eachindex(Cd))
            iF = l == u
            Cd[i] = (!iF)Pi
            xBˈCn += iF * u * Pi
        end; Cd[N+1] = -1
        Ops.add_le_constr!(m.m.moi_backend, N+1, n.Ci, Cd, -xBˈCn) # 🟩
    end
    if Needsˈinit
        JuMP.optimize!(m.m); Ops.assert_optimal(m.m, "mst_init"); Ops.load_mst_sol!(m)
        n = ng.t[1]; nˈXl, nˈrC = n.Xl, n.rC # the LP
        nˈXl .= m.chk.X
        Ops.set_lb!(n.m.moi_backend, 0, N, nˈXl);Ops.set_ub!(n.m.moi_backend, 0, N, nˈXl)
        JuMP.optimize!(n.m); Ops.assert_optimal(n.m, "sub_lp_very_1st")
        Ops.load_F64_vec_RC!(n.m.moi_backend, 0, N, nˈrC)
        ObjVal, newCut = JuMP.objective_value(n.m), [(lbv[i] != ubv[i])RC for (i,RC)=enumerate(nˈrC)]
        BenˈCn = ObjVal - newCut'nˈXl
        n.Cd[1:N] .= newCut; n.Cd[N+1] = -1
        Ops.add_le_constr!(m.m.moi_backend, N+1, n.Ci, n.Cd, -BenˈCn) # 🟩
        push!(vCut, push!(newCut, n.s, BenˈCn)) # 🟩🟩
        Ops.recover_θ_obj_coeff!(m) # 🟦
    end
    NchkCutOff = 0
    for cut_gen = 1:typemax(Int)
        JuMP.optimize!(m.m); Ops.assert_optimal(m.m, "mst"); Ops.load_mst_sol!(m)
        NDˈS.tm += JuMP.solve_time(m.m)
        m.chk.ul[S+1] > ubs && (ND.o[Nnow] = false; NDˈS.np += 1; break) # prune this node
        n = ng.t[1]; nˈXl, nˈrC = n.Xl, n.rC # the LP
        nˈXl .= m.chk.X; Θs_chk = m.chk.Θ[n.s] # local static ✂️
        Ops.set_lb!(n.m.moi_backend, 0, N, nˈXl);Ops.set_ub!(n.m.moi_backend, 0, N, nˈXl)
        JuMP.optimize!(n.m); Ops.assert_optimal(n.m, "sub_lp")
        NDˈS.tn += JuMP.solve_time(n.m)
        ObjVal,_ = JuMP.objective_value(n.m), Ops.load_F64_vec_RC!(n.m.moi_backend, 0, N, nˈrC)
        Vio_s = ObjVal-Θs_chk; rgapCut = abs(Vio_s / ObjVal)
        # println("ObV= $ObjVal, Θs= $Θs_chk, Vio_s= $Vio_s, rgapBen= $rgapCut")
        if rgapCut > Ops.CutRTol.b && Vio_s > 0
            newCut = [(lbv[i] != ubv[i])RC for (i,RC)=enumerate(nˈrC)]
            BenˈCn = ObjVal - newCut'nˈXl
            n.Cd[1:N] .= newCut; n.Cd[N+1] = -1
            Ops.add_le_constr!(m.m.moi_backend, N+1, n.Ci, n.Cd, -BenˈCn) # 🟩
            NchkCutOff += 1; push!(vCut, push!(newCut, n.s, BenˈCn)) # 🟩🟩
            NDˈS.cb += 1; continue # add a Ben cut
        end
        if sTr # attempt to solve MIP subproblem to add sB cut
            n = ng.t[2]; nPi = n.Xl; # the MIP model to derive the cn term
            for (i,RC)=enumerate(nˈrC)  nPi[i]=-(lbv[i] != ubv[i])RC end
            Ops.set_obj_coeff!(n.m.moi_backend, 0, N, nPi)
            Ops.set_lb!(n.m.moi_backend, 0, N, lbv);Ops.set_ub!(n.m.moi_backend, 0, N, ubv)
            JuMP.set_attribute(n.m, "NodeLimit", 2); JuMP.optimize!(n.m); Ops.assert_nodelimit(n.m, "mip_cn")
            NDˈS.tsb += JuMP.solve_time(n.m)
            sBˈCn = JuMP.objective_bound(n.m)
            sBˈub = sum((lbv[i] != ubv[i])nˈXl[i]RC for (i,RC)=enumerate(nˈrC); init=sBˈCn) 
            Vio_s = sBˈub-Θs_chk; rgapCut = Vio_s / sBˈub
            # printstyled("sBˈub= $sBˈub, Θs= $Θs_chk, Vio_s= $Vio_s, rgapSB= $rgapCut\n"; color = 21)
            if rgapCut > Ops.CutRTol.s && Vio_s > 0
                newCut = [(lbv[i] != ubv[i])RC for (i,RC)=enumerate(nˈrC)]
                n.Cd[1:N] .= newCut; n.Cd[N+1] = -1
                Ops.add_le_constr!(m.m.moi_backend, N+1, n.Ci, n.Cd, -sBˈCn) # 🟩
                NchkCutOff += 1; push!(vCut, push!(newCut, n.s, sBˈCn)) # 🟩🟩
                NDˈS.csb += 1; continue # add a sB cut
            end
        end
        Xres .= m.chk.X # load solution info
        break # no more need to explore this node
    end
    ND.l[Nnow] = m.chk.ul[S+1] # load lb info
    ND.o[Nnow] || return NchkCutOff
    Mu = similar(vCut, Cdouble)
    Ops.load_F64_vec_Pi!(m.m.moi_backend, 0, length(Mu), Mu)
    deleteat!(vCut, Mu .== 0)
    return NchkCutOff
end
function select_i_X(Xres, NDˈrˈb)
    (lbv, ubv), (ib, Xb, dmin), Nf = NDˈrˈb, (0, 0., Inf), 0
    for (i, x)=enumerate(Xres)
        (lbv[i] != ubv[i] && Ops.isfrac(x)) || continue
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
    N = 1+Ops.currentMaxN(ND.nc) # Label a new node
    ND.r[N] = tRef = Ops.BnBNodeRef(Nnow,
        map(copy, ND.r[Nnow].b),
        map(copy, ND.r[Nnow].c)
    )
    tRef.b[2-b][iBr] = b # perform a new branching step
    ND.o[N] = ND.e[N] = true
    ND.l[N], ND.s[N] = ND.l[Nnow], use_sB & ND.s[Nnow]
    k = Ops.nXfixed(tRef.b)
    haskey(ND.nc, k) ? push!(ND.nc[k], N) : (ND.nc[k] = [N])
    N
end
function pick_an_open_node(ND, srchMod)
    f, (; nc, o) = ifelse(srchMod, reverse, identity), ND
    for l=f(0:maximum(keys(nc))), n=sort!(nc[l])
        o[n] && return n
    end;  return Ops.N_Dumb
end
function log_gaps(lb, ub; color = 52)
    agap = ub - lb
    rgap = agap / ub
    printstyled("lb= $lb, ub= $ub, agap= $agap, rgap= $rgap\n"; color = color)
    agap, rgap
end

end
