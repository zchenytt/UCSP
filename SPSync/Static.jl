"""
This module is only for generating the history trajectories
"""
module Static
import Random, JuMP, ..Case2383, ..TS

const UpResZones = 1:4
_0() = JuMP.AffExpr()
_1(rng) = rand(rng, 5.0:1e-4:8.0)
_3(rng) = rand(rng, 0.5:1e-4:1.5)

function Model!(m, rng, li0, zi0, F, liD, LD, GN, WN, GMax, Gmin, GEO_LOAD, r1r2c1c2i1i2, LoadYr)
    lines = keys(liD)
    LDˈi, LDˈg, LDˈn = LD.i, LD.g, LD.n
    pfe, egp, rgp = Dict(k=>_0() for k=lines), _0(), [_0() for z=1:4]
    JuMP.@variables(m, begin
        -liD[k] ≤ pf[k=lines] ≤ liD[k]
        x[z=eachindex(GN), g=eachindex(GN[z])], Bin # static problem, no turn-on/off actions
        0 ≤ pA[z=eachindex(GN), g=eachindex(GN[z])]
        0 ≤ ur[z=eachindex(GN), g=eachindex(GN[z])]
        0 ≤ ϖ[z=eachindex(WN)] ≤ 0
        0 ≤ ζ[z=eachindex(LDˈn)] ≤ 0
    end)
    JuMP.@constraint(m, [z=eachindex(GN), g=eachindex(GN[z])], 
        pA[z,g]+ur[z,g] ≤ (GMax[z][g]-Gmin[z][g])x[z,g])
    #= Load =#
        if GEO_LOAD !== nothing
            for (z,v)=enumerate(LDˈg), (g,p)=enumerate(v)
                L01ˈvec, node = r1r2c1c2i1i2[LDˈi[z][g]], LDˈn[z][g]
                load = (p)L01ˈvec[li0] # time_and_geo-variant
                upResDemand = 0.03 * load
                JuMP.add_to_expression!(egp, -load)
                JuMP.add_to_expression!(egp, load, ζ[z])
                z ∈ UpResZones && JuMP.add_to_expression!(rgp[z], -upResDemand)
                z ∈ UpResZones && JuMP.add_to_expression!(rgp[z], upResDemand, ζ[z])
                for l=lines
                    Cnl, pfl = F[l,node], pfe[l]
                    JuMP.add_to_expression!(pfl, Cnl * -load)
                    JuMP.add_to_expression!(pfl, Cnl * load, ζ[z])
                end
            end
        else # no load_geo_discrepancy
            lfra = LoadYr[li0]
            for (z,v)=enumerate(LDˈg)
                ζˈz = ζ[z] # load shed
                z ∈ UpResZones && (rˈz = rgp[z]) 
                for (g,p)=enumerate(v)
                    n = LDˈn[z][g]
                    load = lfra * p # time-variant
                    JuMP.add_to_expression!(egp, -load)
                    JuMP.add_to_expression!(egp, load, ζˈz)
                    if z ∈ UpResZones
                        upResDemand = 0.03 * load
                        JuMP.add_to_expression!(rˈz, -upResDemand)
                        JuMP.add_to_expression!(rˈz, upResDemand, ζˈz)
                    end
                    for l=lines
                        Cnl, pfl = F[l,n], pfe[l]
                        JuMP.add_to_expression!(pfl, Cnl * -load)
                        JuMP.add_to_expression!(pfl, Cnl * load, ζˈz)
                    end
                end
            end
        end
    #= wind =# for (z,v)=enumerate(Case2383.WindP)
        ze = if z == 1
            TS.ze1
        elseif z == 2
            TS.ze2
        elseif z == 3
            TS.ze3
        end
        z_real, rˈz, ϖˈz = ze[zi0], rgp[z], ϖ[z]
        wfra, μfra = TS.logistic(z_real), TS.μ01[z]
        for (g,pMa)=enumerate(v)
            pWind, pavg = wfra * pMa, μfra * pMa
            upResDemand = 0.05 * pavg
            JuMP.add_to_expression!(rˈz, -upResDemand)
            JuMP.add_to_expression!(egp, pWind)
            JuMP.add_to_expression!(egp, -pWind, ϖˈz)
            n = WN[z][g]
            for l=lines
                Cnl, pfl = F[l,n], pfe[l]
                JuMP.add_to_expression!(pfl, Cnl * pWind)
                JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖˈz)
            end
        end
    end
    for (z,v)=enumerate(Gmin)
        z ∈ UpResZones && (rˈz = rgp[z])
        for (g,pmin)=enumerate(v)
            z ∈ UpResZones && JuMP.add_to_expression!(rˈz, ur[z,g])
            JuMP.add_to_expression!(egp, pmin, x[z,g])
            JuMP.add_to_expression!(egp, pA[z,g])
            n = GN[z][g]
            for l=lines
                Cnl, pfl = F[l,n], pfe[l]
                JuMP.add_to_expression!(pfl, Cnl * pmin, x[z,g])
                JuMP.add_to_expression!(pfl, Cnl, pA[z,g])
            end
        end
    end
    JuMP.@constraints(m, begin
        egp == 0
        rgp .>= 0
        [k=lines], pfe[k] == pf[k]
    end)
    JuMP.@objective(m, Min, sum(_1(rng)i for i=x) + sum(_3(rng)i for i=pA))
end

function get_History!(m, rng, li0, zi0, F, liD, LD, GD, WN, GEO_LOAD, r1r2c1c2i1i2, LoadYr)
    GDˈUT = GD.UT
    uL = max(maximum(maximum, GDˈUT),maximum(maximum, GD.DT))-1 # uH is not only for `u` but also for `v`
    uH =  [[falses(uL  ) for _=v] for v=GDˈUT]
    xH =  [[falses(uL+1) for _=v] for v=GDˈUT]
    pAH = [[NaN          for _=v] for v=GDˈUT]
    Model!(m, rng, li0, zi0, F, liD, LD, GD.n, WN, GD.pmax, GD.pmin, GEO_LOAD, r1r2c1c2i1i2, LoadYr)
    JuMP.set_time_limit_sec(m, 20.)
    JuMP.optimize!(m)
    ter = JuMP.primal_status(m)
    ter == JuMP.FEASIBLE_POINT || (println("You should regenerate a physic case"); return xH,pAH)
    pA, x = m[:pA], m[:x]
    for (z,v)=enumerate(uH), (g,bv)=enumerate(v)
        # 🟩🟧 Attention: whenever use the output result of one solve as the input to another solve,
        # you have to correct the result with Julia functions
        # otherwise the downstream solve may become strangly infeasible
        xBval = round(Bool, JuMP.value(x[z,g]))
        pAH[z][g] = xBval * JuMP.value(pA[z,g])
        if xBval
            if rand(rng) ≤ .9
                j = rand(rng, eachindex(bv))
                xH[z][g][1:j] .= bv[j] = true
            else
                xH[z][g] .= true
            end
        else
            if rand(rng) ≤ .9
                j = rand(rng, eachindex(bv))
                xH[z][g][j+1:end] .= true
            end
        end
    end
    xH,pAH # these two is okay
end

end
