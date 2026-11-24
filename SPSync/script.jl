include("src/Bnb.jl");
include("src/Ops.jl");
include("src/Case2383.jl");
include("src/TS.jl");
include("src/Static.jl");
include("src/Sin_Tha_Models.jl");
include("src/Sin_Tha_Cgen.jl");
include("src/SyncAlgo.jl");
import JuMP, Random
R = 7;
function test(SEEDi)
    println("\ntest> SEEDi = $SEEDi")
    T, S = 23, 127;
    rng = Random.Xoshiro(SEEDi);
    LD = Case2383.get_load(rng, Case2383.CaD);
    GD = Case2383.get_gen(rng, Case2383.CaD);
    liD = Case2383.branch_by_zone(Case2383.CaD).P;
    li0, zi0 = TS.def_t0_indices(rng); # absolute location indices
    li1 = li0 + 2(T+1)-1; # absolute terminal of planning: The last program is at li1-static
    xH,pAH = Static.get_History!(Ops.Model(), rng, li0, zi0, Case2383.F, liD, LD, GD, Case2383.WindN);
    pAH[1][1] === NaN && return nothing
    nnt = Sin_Tha_Models.prepare_nnt(S, li0, li1, zi0, T, S, xH, pAH, GD, LD, liD; solve_once = true);
    N = nnt.e[1].N;
    ND, MaxBrDepth = Bnb.NDT(N), 3;
    ubs, Xast, Xres = Inf, fill(Inf, N), fill(Inf, N);
    ato = Threads.Atomic{Float64}(0)
    srchMod = false # (true = attempt to find ub quickly; false = attempt to improve the least lb)
    for search = 1:typemax(Int)
        (Nnow = Bnb.pick_an_open_node(ND, srchMod)) == Bnb.N_Dumb && (printstyled("Quit as all nodes are closed.\n"; color = 93); break)
        NchkCutOff = SyncAlgo.explore_node!(ND, Nnow, ubs, Xres, nnt, Ops.Model(), li0, zi0, S, xH, pAH, GD, LD, liD, ato)
        ND.o[Nnow] || (printstyled("Node $Nnow pruned, NchkCutOff= $NchkCutOff\n"; color=34); continue)
        (iBr, vBr, Nfrac), NxFixed = Bnb.select_i_X(Xres, ND.r[Nnow].b), Bnb.nXfixed(ND.r[Nnow].b)
        println("$search Nnow=$Nnow (Depth=$NxFixed), NchkCutOff=$NchkCutOff; Nfrac=$Nfrac, select=($iBr, $(round(vBr;digits=2)))")
        if iBr == 0 # This Nnow node yields Int-Feas. Sol
            @. Xres = round(Xres); ub = SyncAlgo.get_valid_ub(Xres, nnt.e, S, ato, ND.S)
            ND.o[Nnow] = srchMod = false
            ub < ubs || continue # This solution is inferior
            Xast .= Xres # good enough to update incumbent
            ubs, glb = ub, Bnb.global_lb(ND)
            agap, rgap = Bnb.log_gaps(glb, ubs)
            (0≤agap && 0≤rgap<0.009) || continue
            printstyled("Quit as a good Xast found\n"; color = :magenta)
            break # good enough to quit the algorithm
        end
        _, Nˈ2 = Bnb.branch_and_closeˈnonleaf_Nnow!(ND, Nnow, iBr, vBr, NchkCutOff > 0)
        NxFixed ≥ MaxBrDepth && ((srchMod, ND.o[Nˈ2]) = (true, false))
    end
    Bnb.log_stats(ND.S)
    Bnb.log_gaps(Bnb.global_lb(ND), ubs; color = :magenta)
end;
test_rem(R) = for i = 10000:30000
    i % 10 == R && test(i)
end;
test_rem(R)
