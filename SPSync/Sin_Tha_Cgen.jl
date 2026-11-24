"""
functions to generate all kinds of cuts
"""
module Sin_Tha_Cgen
import ..Ops, JuMP
const CutRTol = (b = .03, s = .001) # Ben 3%, sB .5%

function load_mst_sol!(m)
    m.chk.Θ[1] = JuMP.objective_value(m.m)
    Ops.load_F64_vec!(m.m.moi_backend, 1, m.N, m.chk.X)
end
function get_lp_RCˈCn!(n, lbv, ubv, ato)
    (; m, N, rC, Xl, Cd) = n; o = m.moi_backend
    Ops.set_lb!(o, 0, N, Xl); Ops.set_ub!(o, 0, N, Xl)
    JuMP.optimize!(m); Ops.assert_optimal(m, "sub_lp")
    Threads.atomic_max!(ato, JuMP.solve_time(m))
    _, Obv = Ops.load_F64_vec_RC!(o, 0, N, rC), JuMP.objective_value(m)
    Cd[N+2] = Obv-sum(((l != u)p)x for (l,u,x,p)=zip(lbv,ubv,Xl,rC))
end
function get_sBˈCn!(n1, n2, lbv, ubv, ato)
    (; N, rC) = n1; m = n2.m; o = m.moi_backend
    for ((i,p),l,u)=zip(enumerate(rC), lbv, ubv)
        rC[i] = (l != u) * -1 * p
    end
    Ops.set_obj_coeff!(o,0,N,rC); Ops.set_lb!(o,0,N,lbv); Ops.set_ub!(o,0,N,ubv)
    JuMP.set_attribute(m, "NodeLimit", 2); JuMP.optimize!(m); Ops.assert_nodelimit(m, "mip_cn")
    Threads.atomic_max!(ato, JuMP.solve_time(m))
    n1.Cd[N+2] = sBˈCn = JuMP.objective_bound(m)
end
function ini_cut_gen!(vCut, lbv, ubv, m, nˈtˈ1, ato)
    Pr, N = 1 / m.S, m.N
    for n=nˈtˈ1 n.Xl .= m.chk.X end
    Threads.@threads(for n=nˈtˈ1 get_lp_RCˈCn!(n, lbv, ubv, ato) end)
    n1=nˈtˈ1[1]; Cd=n1.Cd; Pi=view(Cd, 1:N); Pi .= Cn = 0.
    for n=nˈtˈ1
           Cn += Pr * n.Cd[N+2]
        @. Pi += Pr * n.rC
    end
    for (i,l,u)=zip(eachindex(Pi), lbv, ubv)
        l == u && (Pi[i] = 0.)
    end
    push!(vCut, push!(collect(Pi), Cn))
    Cd[N+1]=-1; Ops.add_le_constr!(m.m.moi_backend, N+1, n1.Ci, Cd, -Cn)
end
function Ben_attempt!(vCut, lbv, ubv, m, nˈtˈ1, ato, NDˈS)
    Pr, N = 1 / m.S, m.N
    (; Θ, X) = m.chk; for n=nˈtˈ1 n.Xl .= X end
    ato[] = 0.
    Threads.@threads(for n=nˈtˈ1 get_lp_RCˈCn!(n, lbv, ubv, ato) end)
    NDˈS.tn += ato[]
    n1=nˈtˈ1[1]; Cd=n1.Cd; Pi=view(Cd, 1:N); Pi .= Cn = 0.
    for n=nˈtˈ1
           Cn += Pr * n.Cd[N+2]
        @. Pi += Pr * n.rC
    end
    for (i,l,u)=zip(eachindex(Pi), lbv, ubv)
        l == u && (Pi[i] = 0.)
    end # Here, Pi is derived. And it also applies to sB slope
    ub = Pi'X + Cn; Vio_s = ub - Θ[1]; rgapCut = abs(Vio_s / ub)
    # println("Ben> Vio_s = $Vio_s, rgapCut = $rgapCut")
    if (addBcut = rgapCut > CutRTol.b && Vio_s > 0)
        push!(vCut, push!(collect(Pi), Cn))
        Cd[N+1]=-1; Ops.add_le_constr!(m.m.moi_backend, N+1, n1.Ci, Cd, -Cn)
    end
    addBcut
end
function sB_attempt!(vCut, lbv, ubv, m, nˈt, ato, NDˈS)
    Pr, N = 1 / m.S, m.N; nt1, nt2 = nˈt
    ato[] = 0.
    Threads.@threads(for i=eachindex(nt2) get_sBˈCn!(nt1[i], nt2[i], lbv, ubv, ato) end)
    NDˈS.tsb += ato[]
    Cn = sum(n.Cd[N+2]Pr for n = nt1)
    n1=nt1[1]; Cd=n1.Cd; Pi=view(Cd, 1:N)
    (; Θ, X) = m.chk
    ub = Pi'X + Cn; Vio_s = ub - Θ[1]; rgapCut = abs(Vio_s / ub)
    # println("sB> Vio_s = $Vio_s, rgapCut = $rgapCut")
    if (add_sB = rgapCut > CutRTol.s && Vio_s > 0)
        push!(vCut, push!(collect(Pi), Cn))
        Cd[N+1]=-1; Ops.add_le_constr!(m.m.moi_backend, N+1, n1.Ci, Cd, -Cn)
    end
    add_sB
end

end
