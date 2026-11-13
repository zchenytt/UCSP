include("src/Settings.jl");
include("src/Case2383.jl");
include("src/TS.jl");
include("src/Static.jl");
include("src/TSSP.jl"); # upper level UC ★ 
include("src/F2SP.jl"); # upper level UC ■
include("src/ED.jl"); # lower level ED
using WGLMakie
# http://localhost:9384/browser-display

import JuMP, Random
_gModel() = Settings.Model(Settings.Env());
receding_horizon!(ImplementCost, optTime, mod, li0, li1, zi0, T, S, xH, pAH, GD, LD, liD) = for i = 0:length(ImplementCost)-1
    #= Upper level =#
        m = mod.Model!(_gModel(), li0+i,li1,zi0+i, T, S, xH, pAH, GD, LD, liD; PfK=1.0);
        JuMP.set_attribute(m.m, "MIPGap", 0.5/100)
        JuMP.optimize!(m.m)
        ter = JuMP.termination_status(m.m)
        ter ∈ (JuMP.INTERRUPTED, JuMP.OPTIMAL) || error("upper level ter = $ter")
        optTime[li0+i] += JuMP.solve_time(m.m)
        uc_cost = mod.uv_immediate_cost(m)
        ImplementCost[li0+i] += uc_cost
        mod.update_xH_vec!(m, xH, GD)
    #= Lower level=#
        m = ED.static_lp!(_gModel(), li0+i,zi0+i, xH, pAH, GD, LD, liD; PfK=1.0)
        JuMP.optimize!(m.m)
        ter = JuMP.termination_status(m.m)
        ter == JuMP.OPTIMAL || error("ED ter = $ter")
        optTime[li0+i] += JuMP.solve_time(m.m)
        ed_cost = JuMP.objective_value(m.m)
        ImplementCost[li0+i] += ed_cost
        ED.update_pAH!(pAH, m, GD)
        # println("mod=$mod, No.SP=$i, uc_cost = $uc_cost, ed_cost = $ed_cost")
end;
addsum!(v) = for i=2:length(v) v[i] += v[i-1] end;
# function test(XXSP, SEEDi; NSP = 48, T = 23, S = 3)
function test(XXSP, SEEDi; T = 23, S = 3)
    NSP = 2*(T+1)
    rng = Random.Xoshiro(SEEDi);
    LD = Case2383.get_load(rng, Case2383.CaD)
    GD = Case2383.get_gen(rng, Case2383.CaD)
    liD = Case2383.branch_by_zone(Case2383.CaD).P
    li0, zi0 = TS.def_t0_indices(rng); # absolute location indices
    li1 = li0 + 2(T+1)-1 # absolute terminal of planning: The last program is at li1-static
    xH,pAH = Static.get_History!(_gModel(), rng, li0, zi0, Case2383.F, liD, LD, GD, Case2383.WindN);
    ImplementCost = Dict(i => 0. for i = range(li0; length=NSP));
    optTime = Dict(i => 0. for i = range(li0; length=NSP));
    receding_horizon!(ImplementCost, optTime, XXSP, li0, li1, zi0, T, S, xH, pAH, GD, LD, liD)
    totcost = sum(values(ImplementCost))
    println("total cost of $XXSP is $totcost")
    (; XXSP, SEEDi, NSP, T, S, LD, GD, liD, li0, li1, zi0, xH, pAH, ImplementCost, optTime)
end;

# 🟠🟧 This can only be run once
for SEEDi = 0:typemax(Int)
    println("SEEDi = $SEEDi")
    fint = test(F2SP, SEEDi)
    mynt = test(TSSP, SEEDi)
end

# 🟠🟧 This can only be run once

rang = range(mynt.li0;length=mynt.NSP);
fiv = [fint.ImplementCost[i] for i = rang];
myv = [mynt.ImplementCost[i] for i = rang];
addsum!(fiv);addsum!(myv);

# TODO figure out why is fin-binding more cost effecive in RH optimization.

WGLMakie.activate!(resize_to = :body);
f = Figure();
ax = Axis(f[1,1]);

scatterlines!(ax, rang, fiv; marker = :star4, color = :blue);
scatterlines!(ax, rang, myv; marker = :xcross, color = :green);

f
