include("src/Settings.jl");
include("src/Case2383.jl");
include("src/TS.jl");
include("src/Static.jl");
include("src/TSSP.jl"); # upper level UC ★ 
include("src/ED.jl"); # lower level ED

import JuMP, Random
_gModel() = Settings.Model(Settings.Env());

receding_horizon!(ImplementCost, optTime, li0, zi0, T, S, xH, pAH, GD, LD, liD) = for i = 0:length(ImplementCost)-1
    #= Upper level =#
        m = TSSP.x_now!(_gModel(), li0+i,zi0+i, T, S, xH, pAH, GD, LD, liD; PfK=1.0);
        JuMP.set_attribute(m.m, "MIPGap", 0.5/100)
        JuMP.optimize!(m.m); JuMP.termination_status(m.m) ∈ (JuMP.INTERRUPTED, JuMP.OPTIMAL) || error("2SSP_UC error")
        println("No.SP = $i, time = $(JuMP.solve_time(m.m))")
        optTime[li0+i] += JuMP.solve_time(m.m)
        ImplementCost[li0+i] += TSSP.uv_immediate_cost(m)
        TSSP.update_xH_vec!(m, xH, GD)
    #= Lower level=#
        m = ED.static_lp!(_gModel(), li0+i,zi0+i, xH, pAH, GD, LD, liD; PfK=1.0)
        JuMP.optimize!(m.m); JuMP.termination_status(m.m) == JuMP.OPTIMAL || error("ED error")
        optTime[li0+i] += JuMP.solve_time(m.m)
        ImplementCost[li0+i] += JuMP.objective_value(m.m)
        ED.update_pAH!(pAH, m, GD)
end;

SEEDi = 1; # Case SEEDi=3 is very hard
NSP = 24 * 5 # Num of SPs: one week
T = 23 # Num of look-ahead periods per SP
S = 3 # Num of scenarios used to construct each SP
rng = Random.Xoshiro(0);
LD = Case2383.get_load(rng, Case2383.CaD);
GD = Case2383.get_gen(rng, Case2383.CaD);
liD = Case2383.branch_by_zone(Case2383.CaD).P;
li0, zi0 = TS.def_t0_indices(rng); # absolute location indices
xH,pAH = Static.get_History!(_gModel(), rng, li0, zi0, Case2383.F, liD, LD, GD, Case2383.WindN);

ImplementCost = Dict(i => 0. for i = range(li0; length=NSP));
optTime = Dict(i => 0. for i = range(li0; length=NSP));
receding_horizon!(ImplementCost, optTime, li0, zi0, T, S, xH, pAH, GD, LD, liD)
