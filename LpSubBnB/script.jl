include("src/Settings.jl");
include("src/Case2383.jl");
include("src/TS.jl");
include("src/Static.jl");
include("src/Ops.jl");
include("src/QsMst.jl");
import JuMP, Random
_gModel() = Settings.Model(Settings.Env());
begin
    SEEDi, T, S = 3, 23, 1
    rng = Random.Xoshiro(SEEDi);
    LD = Case2383.get_load(rng, Case2383.CaD);
    GD = Case2383.get_gen(rng, Case2383.CaD);
    liD = Case2383.branch_by_zone(Case2383.CaD).P;
    li0, zi0 = TS.def_t0_indices(rng); # absolute location indices
    li1 = li0 + 2(T+1)-1; # absolute terminal of planning: The last program is at li1-static
    xH,pAH = Static.get_History!(_gModel(), rng, li0, zi0, Case2383.F, liD, LD, GD, Case2383.WindN);
end;
n = QsMst.proto!(1, _gModel(), li0, li1, zi0, T, S, xH, pAH, GD, LD, liD);

# Global containers
    NInt = 0; # no node at the beginning
    NDO = Dict{Int, Bool}(); # true = open (still exploring), false = closed (e.g. becoming non-leaf)
    NDL = Dict{Int, Float64}(); # lower bounds
    NDR = Dict{Int, Ops.BnBNodeRef}();
NInt += 1; # create a new node
Nnow = NInt; # currently exploring Root node
NDO[Nnow] = true; # node Nnow is at OPEN status
NDR[Nnow] = Ops.BnBNodeRef(0,fill(Cdouble(0),n.N),fill(Cdouble(1),n.N)); # node dict for reference purpose (thus should never delete item)

# Root node
    m = QsMst.mst!(_gModel(), li0, zi0, S, xH, pAH, GD, LD, liD; θ_has_obj = false);
    Ops.set_lb!(m.m.moi_backend, S, m.N, NDR[Nnow].l)
    Ops.set_ub!(m.m.moi_backend, S, m.N, NDR[Nnow].u)
    JuMP.optimize!(m.m); Ops.assert_optimal(m.m, "mst")
    Ops.load_mst_sol!(m)
    n.Xl .= m.chk.X;
    Ops.set_lb!(n.m.moi_backend, 0, n.N, n.Xl)
    Ops.set_ub!(n.m.moi_backend, 0, n.N, n.Xl)
    JuMP.optimize!(n.m); Ops.assert_optimal(n.m, "sub")
    Ops.load_F64_vec_RC!(n.m.moi_backend, 0, n.N, n.Pi)
    BenˈCn=ObjVal=JuMP.objective_value(n.m)
    iCi = 0
    for (X, RC, giX, ni) = zip(n.Xl, n.Pi, range(S; length=m.N), 1:m.N)
        if NDR[Nnow].l[ni] != NDR[Nnow].u[ni]
            BenˈCn -= (X)RC
            abs(RC) < 1e-13 && continue
            iCi += 1; n.Ci[iCi], n.Cd[iCi] = giX, RC
        end
    end 
    n.Ci[iCi+1],n.Cd[iCi+1]=Ops.ni2gi4s(n.s),-1 # send to GRB
    Ops.add_le_constr!(m.m.moi_backend, iCi+1, n.Ci, n.Cd, -BenˈCn) # 🟩🟩
    Ops.recover_θ_obj_coeff!(m) # 🟩
    for k = 1:typemax(Int)
        JuMP.optimize!(m.m); Ops.assert_optimal(m.m, "mst")
        Ops.load_mst_sol!(m)
        print("Nnow=$Nnow, k=$k, lb=$(m.chk.ul[S+1]), ")
        n.Xl .= m.chk.X;
        Ops.set_lb!(n.m.moi_backend, 0, n.N, n.Xl)
        Ops.set_ub!(n.m.moi_backend, 0, n.N, n.Xl)
        JuMP.optimize!(n.m); Ops.assert_optimal(n.m, "sub")
        BenˈCn=ObjVal=JuMP.objective_value(n.m)
        Vio_s = ObjVal-m.chk.Θ[n.s]; println("Vio_s=$Vio_s")
        if Vio_s > 1e-6
            Ops.load_F64_vec_RC!(n.m.moi_backend, 0, n.N, n.Pi)
            iCi = 0
            for (X, RC, giX, ni) = zip(n.Xl, n.Pi, range(S; length=m.N), 1:m.N)
                if NDR[Nnow].l[ni] != NDR[Nnow].u[ni]
                    BenˈCn -= (X)RC
                    abs(RC) < 1e-13 && continue
                    iCi += 1; n.Ci[iCi], n.Cd[iCi] = giX, RC
                end
            end 
            n.Ci[iCi+1],n.Cd[iCi+1]=Ops.ni2gi4s(n.s),-1 # send to GRB
            Ops.add_le_constr!(m.m.moi_backend, iCi+1, n.Ci, n.Cd, -BenˈCn) # 🟩🟩
        else
            NDL[Nnow] = m.chk.ul[S+1] # record the lower bound
            break
        end
    end

# Observe "m.chk.X" and then create new nodes by branching
    NInt += 1; # create a new node, branchAct = (3 => true)
    NDO[NInt] = true; # mark this new node as OPEN
    lbV, ubV = copy(NDR[Nnow].l), copy(NDR[Nnow].u); # inherit!
    ubV[3] = false; # The incremental branching act
    NDR[NInt] = Ops.BnBNodeRef(Nnow,lbV,ubV);
    NInt += 1; # create a new node, branchAct = (3 => false)
    NDO[NInt] = true; # mark this new node as OPEN
    lbV, ubV = copy(NDR[Nnow].l), copy(NDR[Nnow].u); # inherit!
    lbV[3] = true; # The incremental branching act
    NDR[NInt] = Ops.BnBNodeRef(Nnow,lbV,ubV);
    NDO[Nnow] = false; # as we've partitioned the father node, close it

Nnow = 3; # 🟧 select an open node randomly
