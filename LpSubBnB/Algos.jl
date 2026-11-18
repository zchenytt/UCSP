module Algos
import ..QsMst, ..Ops, JuMP

function explore_node!(ND, Nnow, Xres, n, new_m, li0, zi0, S, xH, pAH, GD, LD, liD)
    m = QsMst.mst!(new_m, li0, zi0, S, xH, pAH, GD, LD, liD; θ_has_obj = false);
    lbv, ubv = ND.r[Nnow].b
    Ops.set_lb!(m.m.moi_backend, S, m.N, lbv)
    Ops.set_ub!(m.m.moi_backend, S, m.N, ubv)
    JuMP.optimize!(m.m); Ops.assert_optimal(m.m, "mst"); Ops.load_mst_sol!(m)
    n.Xl .= m.chk.X;
    Ops.set_lb!(n.m.moi_backend, 0, n.N, n.Xl)
    Ops.set_ub!(n.m.moi_backend, 0, n.N, n.Xl)
    JuMP.optimize!(n.m); Ops.assert_optimal(n.m, "sub")
    Ops.load_F64_vec_RC!(n.m.moi_backend, 0, n.N, n.Pi)
    BenˈCn, iCi = JuMP.objective_value(n.m), 0
    for (X, RC, giX, ni) = zip(n.Xl, n.Pi, range(S; length=m.N), 1:m.N)
        lbv[ni] != ubv[ni] || continue
        BenˈCn -= (X)RC
        abs(RC) > 1e-13 || continue
        iCi += 1; n.Ci[iCi], n.Cd[iCi] = giX, RC
    end 
    n.Ci[iCi+1],n.Cd[iCi+1]=Ops.ni2gi4s(n.s),-1 # send to GRB
    Ops.add_le_constr!(m.m.moi_backend, iCi+1, n.Ci, n.Cd, -BenˈCn) # 🟩🟩
    Ops.recover_θ_obj_coeff!(m) # 🟦
    for k = 1:typemax(Int)
        JuMP.optimize!(m.m); Ops.assert_optimal(m.m, "mst"); Ops.load_mst_sol!(m)
        # print("Nnow=$Nnow, k=$k, lb=$(m.chk.ul[S+1]), ")
        n.Xl .= m.chk.X;
        Ops.set_lb!(n.m.moi_backend, 0, n.N, n.Xl)
        Ops.set_ub!(n.m.moi_backend, 0, n.N, n.Xl)
        JuMP.optimize!(n.m); Ops.assert_optimal(n.m, "sub")
        BenˈCn, iCi = JuMP.objective_value(n.m), 0
        Vio_s = BenˈCn-m.chk.Θ[n.s];
        # println("Vio_s=$Vio_s")
        Vio_s > 1e-6 || break
        Ops.load_F64_vec_RC!(n.m.moi_backend, 0, n.N, n.Pi)
        for (X, RC, giX, ni) = zip(n.Xl, n.Pi, range(S; length=m.N), 1:m.N)
            lbv[ni] != ubv[ni] || continue
            BenˈCn -= (X)RC
            abs(RC) > 1e-13 || continue
            iCi += 1; n.Ci[iCi], n.Cd[iCi] = giX, RC
        end 
        n.Ci[iCi+1],n.Cd[iCi+1]=Ops.ni2gi4s(n.s),-1 # send to GRB
        Ops.add_le_constr!(m.m.moi_backend, iCi+1, n.Ci, n.Cd, -BenˈCn) # 🟩🟩
    end
    ND.l[Nnow] = m.chk.ul[S+1] # load lb info
    Xres .= m.chk.X # load solution info
    JuMP.objective_value(n.m) # TODO load some other info
end
function make_br_decision(Xres)
    ib, Xb, dmin = 0, 0., Inf
    for (i, x)=enumerate(Xres)
        0<x<1 || continue
        d = abs(x - .5)
        d < dmin && (ib = i; Xb = x; dmin = d)
    end
    ib, round(Bool, Xb) # if ib==0, Xres should be Int-Feasible
end
function branch_and_closeˈnonleaf_Nnow!(ND, Nnow, iBr, vBr)
    new_a_node!(ND, Nnow, iBr,  vBr) # small label to more promising branch
    new_a_node!(ND, Nnow, iBr, !vBr) # the other possibility
    ND.o[Nnow] = ND.e[Nnow] = false
end
function new_a_node!(ND, Nnow, iBr, vBr)
    Nint = 1+Ops.currentMaxN(ND.nc); # Label a new node
    ND.r[Nint] = tRef = Ops.BnBNodeRef(Nnow, map(copy, ND.r[Nnow].b));
    tRef.b[2-vBr][iBr] = vBr; # perform a new branching step
    ND.o[Nint] = ND.e[Nint] = true
    push!(ND.nc[Ops.nXfixed(tRef.b)], Nint)
end
function pick_an_open_node(ND)
    (; nc, o) = ND
    for l=0:maximum(keys(nc)), n=nc[l]
        o[n] && return n
    end;        return 0
end

end
