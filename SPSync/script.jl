include("Bnb.jl");
include("Ops.jl");
include("Sin_Tha_Cgen.jl");
include("Case2383.jl");
include("TS.jl"); # takes an arg
include("Static.jl"); # depends on TS
include("Sin_Tha_Models.jl"); # depends on TS
include("Mono.jl"); # depends on TS
include("SyncAlgo.jl"); # depends on Sin_Tha_Models

import JuMP, Random, JLD2
function gen_input_data(SEEDi, use_geo_load, T, PfK)
    liD = Case2383.branch_by_zone(Case2383.CaD).P;
    rng, GEO_LOAD = Random.Xoshiro(SEEDi), (use_geo_load ? true : nothing);
    r1r2c1c2i1i2, LoadYr = TS.load_patterns(GEO_LOAD)
    LD = Case2383.get_load(rng, Case2383.CaD);
    GD = Case2383.get_gen(rng, Case2383.CaD);
    li0, zi0 = TS.def_t0_indices(rng, GEO_LOAD); # absolute location indices
    xH,pAH = Static.get_History!(Ops.Model(), rng, li0, zi0, Case2383.F, liD, LD, GD, Case2383.WindN, GEO_LOAD, r1r2c1c2i1i2, LoadYr)
    pAH[1][1] === NaN && return nothing
    for i=keys(liD) liD[i] *= PfK end
    # start to build the 2SSP-Decomposition model
    li1 = li0 + 2(T+1)-1; # absolute terminal of planning: The last program is at li1-static
    li0, li1, zi0, xH, pAH, GD, LD, liD, r1r2c1c2i1i2, LoadYr
end;
function get_a_2SSP_result(S, T, li0, li1, zi0, xH, pAH, GD, LD, liD, r1r2c1c2i1i2, LoadYr)
    nnt = Sin_Tha_Models.prepare_nnt(S, T, li0, li1, zi0, xH, pAH, GD, LD, liD, r1r2c1c2i1i2, LoadYr; solve_once = true);
    N = nnt.e[1].N;
    ND, MaxBrDepth = Bnb.NDT(N), 3;
    ubs, Xast, Xres, goodBit = Inf, fill(Inf, N), fill(Inf, N), false;
    ato = Threads.Atomic{Float64}(0)
    Nnxt = Bnb.N_Root # initially, explore Root
    for search = 1:typemax(Int)
        (Nnow = Bnb.pick_an_open_node(ND, Nnxt)) == Bnb.N_Dumb && (printstyled("Quit as all nodes are closed.\n"; color = 93); break)
        NchkCutOff = SyncAlgo.explore_node!(ND, Nnow, ubs, Xres, nnt, Ops.Model(), xH, pAH, GD, ato, S)
        Nnxt = 0 # default is to strengthen least lb
        ND.o[Nnow] || (printstyled("Node $Nnow pruned, NchkCutOff= $NchkCutOff\n"; color=34); continue)
        (iBr, vBr, Nfrac), NxFixed = Bnb.select_i_X(Xres, ND.r[Nnow].b), Bnb.nXfixed(ND.r[Nnow].b)
        # println("$search Nnow=$Nnow (Depth=$NxFixed), NchkCutOff=$NchkCutOff; Nfrac=$Nfrac, select=($iBr, $(round(vBr;digits=2)))")
        if iBr == 0 # This Nnow node yields Int-Feas. Sol
            @. Xres = round(Xres); ub = SyncAlgo.get_valid_ub(Xres, nnt.e, S, ato, ND.S)
            ND.o[Nnow] = false
            ub < ubs || continue # This solution is inferior
            Xast .= Xres # good enough to update incumbent
            ubs, glb = ub, Bnb.global_lb(ND)
            agap, rgap = Bnb.log_gaps(glb, ubs)
            (0≤ agap && 0≤ rgap <0.009) || continue
            # printstyled("Quit as a good Xast found\n"; color = :magenta)
            goodBit = true; break # good enough to quit the algorithm
        end
        Nˈ1, Nˈ2 = Bnb.branch_and_closeˈnonleaf_Nnow!(ND, Nnow, iBr, vBr, NchkCutOff > 0)
        NxFixed < MaxBrDepth || ((Nnxt, ND.o[Nˈ2]) = (Nˈ1, false))
    end
    ND.S.lb = glb = Bnb.global_lb(ND)
    ND.S.ub = ubs
    goodBit || Bnb.log_gaps(glb, ubs; color = 93)
    Bnb.log_stats(ND.S)
    ND.S
end;
test!(rst) = for S = (1, 5, 127, 2000), use_geo_load = (true,false)
    GEO_LOAD = (use_geo_load ? true : nothing)
    r1r2c1c2i1i2, LoadYr = TS.load_patterns(GEO_LOAD)
    liD0 = Case2383.branch_by_zone(Case2383.CaD).P;
    for SEEDi = 0:39
        println("        S = $S, GeoLoad = $(Int(use_geo_load)), SEEDi = $SEEDi")
        T = 23
        rng = Random.Xoshiro(SEEDi)
        LD = Case2383.get_load(rng, Case2383.CaD);
        GD = Case2383.get_gen(rng, Case2383.CaD);
        li0, zi0 = TS.def_t0_indices(rng, GEO_LOAD); # absolute location indices
        li1 = li0 + 2(T+1)-1; # absolute terminal of planning: The last program is at li1-static
        xH,pAH = Static.get_History!(Ops.Model(), rng, li0, zi0, Case2383.F, liD0, LD, GD, Case2383.WindN, GEO_LOAD, r1r2c1c2i1i2, LoadYr)
        pAH[1][1] === NaN && continue
        for PfK = (.75, 1.)
            println("    PfK = $PfK")
            liD = copy(liD0)
            for i=keys(liD) liD[i] *= PfK end
            if S < 6
                t_lu = Mono.get_a_2SSP_result(S, T, li0, li1, zi0, xH, pAH, GD, LD, liD, r1r2c1c2i1i2, LoadYr)
                rst[S, use_geo_load, SEEDi, PfK, false] = t_lu
            end
            stAts = get_a_2SSP_result(S, T, li0, li1, zi0, xH, pAH, GD, LD, liD, r1r2c1c2i1i2, LoadYr)
            rst[S, use_geo_load, SEEDi, PfK, true] = stAts
        end
    end
end;

rst = Dict{Tuple{Int, Bool, Int, Float64, Bool}, Union{Tuple{Float64, Float64}, Bnb.Stats}}()
test!(rst)
JLD2.save("S12000.jld2", "d", rst)
