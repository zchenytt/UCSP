include("src/Case2383.jl");
include("src/TS.jl");
include("src/Static.jl");
include("src/Ops.jl");
include("src/QsMst.jl");
include("src/Algos.jl");

import JuMP, Random
SEEDi, T, S = 3, 23, 1
function test(SEEDi, MaxDepth)
    println("Test Case> SEED = $SEEDi")
    rng = Random.Xoshiro(SEEDi);
    LD = Case2383.get_load(rng, Case2383.CaD);
    GD = Case2383.get_gen(rng, Case2383.CaD);
    liD = Case2383.branch_by_zone(Case2383.CaD).P;
    li0, zi0 = TS.def_t0_indices(rng); # absolute location indices
    li1 = li0 + 2(T+1)-1; # absolute terminal of planning: The last program is at li1-static
    xH,pAH = Static.get_History!(Ops.Model(), rng, li0, zi0, Case2383.F, liD, LD, GD, Case2383.WindN);
    n = QsMst.proto!(1, Ops.Model(), li0, li1, zi0, T, S, xH, pAH, GD, LD, liD);
    for MaxNxFixed = 0:MaxDepth
        ND = Ops.NDT(n.N, MaxNxFixed);
        ubs,hashOfIntSols,Xres = Inf,UInt[],fill(NaN, n.N);
        for eXpl = 1:typemax(Int)
            Nnow = Algos.pick_an_open_node(ND)
            if Nnow == 0
                # println("Nnow = $Nnow, all nodes are closed.")
                break
            end
            ObjVal = Algos.explore_node!(ND, Nnow, Xres, n, Ops.Model(), li0, zi0, S, xH, pAH, GD, LD, liD)
            iBr, vBr = Algos.make_br_decision(Xres)
            if iBr == 0 # encounter an Int-Feas. Sol
                push!(hashOfIntSols, hash(round.(Bool, Xres)))
                if ObjVal < ubs
                    ubs = ObjVal
                    # println("Found a better Int-Feas Solution, ubs=$ubs")
                end
                ND.o[Nnow] = false; continue
            end
            if Ops.nXfixed(ND.r[Nnow].b)==maximum(keys(ND.nc))
                # println("Node is still Int-infeas., close forcibly due to Max Br Depth Reached")
                ND.o[Nnow] = false; continue
            end
            Algos.branch_and_closeˈnonleaf_Nnow!(ND, Nnow, iBr, vBr)
        end
        glb = Ops.global_lb(ND)
        agap = ubs-glb
        printstyled("MaxNxFixed=$MaxNxFixed, ubs=$ubs, glb=$glb, agap=$agap\n"; color = :magenta)
        abs(agap) < 0.05 && break
    end
end;
for i = 0:40 test(i, 8) end

# Results:
julia> for i = 0:40 test(i, 8) end
Test Case> SEED = 0
MaxNxFixed=0, ubs=Inf, glb=4487.533087216761, agap=Inf
MaxNxFixed=1, ubs=Inf, glb=4488.681465144553, agap=Inf
MaxNxFixed=2, ubs=4488.687859538824, glb=4488.687859538824, agap=0.0
Test Case> SEED = 1
MaxNxFixed=0, ubs=7224.562062437532, glb=7224.562062437533, agap=-9.094947017729282e-13
Test Case> SEED = 2
MaxNxFixed=0, ubs=Inf, glb=6172.612701333131, agap=Inf
MaxNxFixed=1, ubs=Inf, glb=6172.612976997613, agap=Inf
MaxNxFixed=2, ubs=6172.824283225421, glb=6172.824283225421, agap=0.0
Test Case> SEED = 3
MaxNxFixed=0, ubs=Inf, glb=1254.8864165807247, agap=Inf
MaxNxFixed=1, ubs=Inf, glb=1277.0052629046554, agap=Inf
MaxNxFixed=2, ubs=Inf, glb=1278.6744734347537, agap=Inf
MaxNxFixed=3, ubs=Inf, glb=1278.7039574424202, agap=Inf
MaxNxFixed=4, ubs=Inf, glb=1278.7485286035342, agap=Inf
MaxNxFixed=5, ubs=1278.9571147693384, glb=1278.9571147693384, agap=0.0
Test Case> SEED = 4
MaxNxFixed=0, ubs=Inf, glb=7207.105023416946, agap=Inf
MaxNxFixed=1, ubs=Inf, glb=7207.259377788028, agap=Inf
MaxNxFixed=2, ubs=Inf, glb=7207.260405495851, agap=Inf
MaxNxFixed=3, ubs=7207.272437187335, glb=7207.261402865845, agap=0.011034321489205468
Test Case> SEED = 5
MaxNxFixed=0, ubs=4270.198891433941, glb=4270.19889143394, agap=9.094947017729282e-13
Test Case> SEED = 6
MaxNxFixed=0, ubs=Inf, glb=3169.3286429107866, agap=Inf
MaxNxFixed=1, ubs=Inf, glb=3169.328986803071, agap=Inf
MaxNxFixed=2, ubs=Inf, glb=3169.329644849912, agap=Inf
MaxNxFixed=3, ubs=3169.5165398414315, glb=3169.3296464953078, agap=0.18689334612372477
MaxNxFixed=4, ubs=3169.3974180687605, glb=3169.330048931129, agap=0.06736913763143093
MaxNxFixed=5, ubs=3169.3974180687605, glb=3169.330121400167, agap=0.06729666859337158
MaxNxFixed=6, ubs=3169.3974180687615, glb=3169.330148938401, agap=0.06726913036027327
MaxNxFixed=7, ubs=3169.3974180687615, glb=3169.3302373506285, agap=0.06718071813293136
MaxNxFixed=8, ubs=3169.397418068761, glb=3169.3302677876254, agap=0.06715028113558219
Test Case> SEED = 7
MaxNxFixed=0, ubs=Inf, glb=5074.33097468098, agap=Inf
MaxNxFixed=1, ubs=Inf, glb=5074.361725897526, agap=Inf
MaxNxFixed=2, ubs=5074.370893440971, glb=5074.361725897522, agap=0.009167543448711513
Test Case> SEED = 8
MaxNxFixed=0, ubs=Inf, glb=6807.6237938281, agap=Inf
MaxNxFixed=1, ubs=Inf, glb=6807.624187936105, agap=Inf
MaxNxFixed=2, ubs=Inf, glb=6807.625651504102, agap=Inf
MaxNxFixed=3, ubs=Inf, glb=6807.625976596096, agap=Inf
MaxNxFixed=4, ubs=Inf, glb=6807.626614813091, agap=Inf
MaxNxFixed=5, ubs=Inf, glb=6807.63139806309, agap=Inf
MaxNxFixed=6, ubs=6807.810645087717, glb=6807.635626110864, agap=0.17501897685269796
MaxNxFixed=7, ubs=6807.810645087715, glb=6807.637960079628, agap=0.17268500808677345
MaxNxFixed=8, ubs=6807.810645087711, glb=6807.6692868076325, agap=0.14135828007874807
