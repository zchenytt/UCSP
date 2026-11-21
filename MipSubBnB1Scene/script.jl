include("src/Case2383.jl");
include("src/TS.jl");
include("src/Static.jl");
include("src/Ops.jl");
include("src/QsMst.jl");
include("src/Algos.jl");
import JuMP, Random
function test()
    SEEDi = rand(Int)
    println("\ntest> SEEDi = $SEEDi")
    T, S = 23, 1;
    rng = Random.Xoshiro(SEEDi);
    LD = Case2383.get_load(rng, Case2383.CaD);
    GD = Case2383.get_gen(rng, Case2383.CaD);
    liD = Case2383.branch_by_zone(Case2383.CaD).P;
    li0, zi0 = TS.def_t0_indices(rng); # absolute location indices
    li1 = li0 + 2(T+1)-1; # absolute terminal of planning: The last program is at li1-static
    xH,pAH = Static.get_History!(Ops.Model(), rng, li0, zi0, Case2383.F, liD, LD, GD, Case2383.WindN);
    pAH[1][1] === NaN && return nothing
    ng = ( # for a particular scene
        t = QsMst.train_tp(1, Ops.Model(), Ops.Model(), li0, li1, zi0, T, S, xH, pAH, GD, LD, liD),
        e = QsMst.proto!(1, Ops.Model(), li0, li1, zi0, T, S, xH, pAH, GD, LD, liD; binary=true)
    );
    N = ng.e.N;
    ND, MaxBrDepth = Ops.NDT(N), 3;
    ubs, Xast, Xres = Inf, fill(Inf, N), fill(Inf, N);
    srchMod = false # (true = attempt to find ub quickly; false = attempt to improve the least lb)
    for search = 1:typemax(Int)
        (Nnow = Algos.pick_an_open_node(ND, srchMod)) == Ops.N_Dumb && (printstyled("Quit as all nodes are closed.\n"; color = 93); break)
        NchkCutOff = Algos.explore_node!(ND, Nnow, ubs, Xres, ng, Ops.Model(), li0, zi0, S, xH, pAH, GD, LD, liD)
        ND.o[Nnow] || (printstyled("Node $Nnow pruned, NchkCutOff= $NchkCutOff\n"; color=93); continue)
        (iBr, vBr, Nfrac), NxFixed = Algos.select_i_X(Xres, ND.r[Nnow].b), Ops.nXfixed(ND.r[Nnow].b)
        println("$search Nnow=$Nnow (Depth=$NxFixed), NchkCutOff=$NchkCutOff; Nfrac=$Nfrac, select=($iBr, $vBr)")
        if iBr == 0 # This Nnow node yields Int-Feas. Sol
            @. Xres = round(Xres); ub = Algos.get_valid_ub(ND.S, Xres, ng.e)
            ND.o[Nnow] = srchMod = false
            ub < ubs || continue # This solution is inferior
            Xast .= Xres # good enough to update incumbent
            ubs, glb = ub, Ops.global_lb(ND)
            agap, rgap = Algos.log_gaps(glb, ubs)
            (0≤agap && 0≤rgap<0.009) || continue
            printstyled("Quit as a good Xast found\n"; color = :magenta)
            break # good enough to quit the algorithm
        end
        _, Nˈ2 = Algos.branch_and_closeˈnonleaf_Nnow!(ND, Nnow, iBr, vBr, NchkCutOff > 0)
        NxFixed ≥ MaxBrDepth && ((srchMod, ND.o[Nˈ2]) = (true, false))
    end
    Ops.log_stats(ND.S)
    Algos.log_gaps(Ops.global_lb(ND), ubs; color = :magenta)
end

while true
    test()
end
