"""
This module (oracle) corresonds to the lower-level ED program in practical power system operation

But in real situations there're multiple EDs, e.g. one run every 15 minutes.
Since lower-level ED is not the main focus of this project, we just implement an easy oracle.
"""
module ED
import Random, JuMP, ..Case2383, ..Settings, ..TS
const UpResZones = 1:4
const ϖCostCoeff = 2.0
const ζ2ϖCostRatio = 25

_0() = JuMP.AffExpr()
function _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
    Ru, Rd = Gwid/rui, Gwid/rdi
    Su, Sd = clamp(r2s * Ru, Pmin, PMax), clamp(r2s * Rd, Pmin, PMax)
    Ru, Rd, Su, Sd
end
rB(x) = round(Bool, x)

function static_lp!(m, li0, zi0, xH, pAH, GD, LD, liD; ϖCostCoeff=ϖCostCoeff, ζ2ϖCostRatio=ζ2ϖCostRatio, PfK=1.0, addPfLim=true)
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci, GDˈC = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.dci, GD.C
    LDˈn, LDˈg, F, ζCostCoeff, WN = LD.n, LD.g, Case2383.F, (ζ2ϖCostRatio)ϖCostCoeff, Case2383.WindN
    Qˈy, egp, lines = _0(), _0(), keys(liD)
    pfe, rgp = Dict(l => _0() for l=lines), [_0() for z=UpResZones]
    JuMP.@variables(m, begin
        0 ≤ pA[z=eachindex(GDˈn), g=eachindex(GDˈn[z])]
        0 ≤ ur[z=eachindex(GDˈn), g=eachindex(GDˈn[z])]
        0 ≤ pI[z=eachindex(GDˈn), g=eachindex(GDˈn[z])]
        0 ≤  ϱ[z=eachindex(GDˈn), g=eachindex(GDˈn[z])]
        0           ≤   ϖ[z=eachindex(WN)  ] ≤ 1
        0           ≤   ζ[z=eachindex(LDˈn)] ≤ 1
        -liD[l]PfK  ≤  pf[l=lines]           ≤ liD[l]PfK
    end)
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            Gwid, xˈzg, xHˈzg = PMax-Pmin, xH[z][g][1], xH[z][g][2]
            if rB(xHˈzg) == rB(xˈzg)
                uˈzg = vˈzg = false
            elseif rB(xˈzg)
                uˈzg, vˈzg = true, false
            else
                uˈzg, vˈzg = false, true
            end
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC, sdC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC, GDˈdci[z][g]linC
            Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        #= Def Relation =#
            JuMP.@constraint(m, pI[z,g] + ϱ[z,g] == (Pmin)xˈzg + pA[z,g])
        #= Costs =#
            JuMP.add_to_expression!(Qˈy, linC, pI[z,g]) # Gen cost
            JuMP.add_to_expression!(Qˈy, linC, ϱ[z,g])  # Gen cost
            JuMP.add_to_expression!(Qˈy, ϖCostCoeff, ϱ[z,g]) # Gen_Curtail cost
        #= Reserve =#
            z ∈ UpResZones && JuMP.add_to_expression!(rgp[z], ur[z,g])
        #= egp / pfl =#
            JuMP.add_to_expression!(egp, pI[z,g])
            for l=lines JuMP.add_to_expression!(pfe[l], F[l,node], pI[z,g]) end
        #= Ramp / Output_Ub =#
            if MmRatio > 1.1
                Ju, Jv = Su - Pmin, Sd - Pmin - Rd
                JuMP.@constraint(m, pA[z,g] + ur[z,g] - pAHˈzg ≤ Ru * xHˈzg + Ju * uˈzg)
                JuMP.@constraint(m, pAHˈzg - pA[z,g]           ≤ Rd * xHˈzg + Jv * vˈzg)
            end
            JuMP.@constraint(m, pA[z,g] + ur[z,g] ≤ Gwid * xˈzg)
    end
    #= Load =#
        for (z,v)=enumerate(LDˈg)
            z ∈ UpResZones && (rˈz = rgp[z]) 
            for (g,p)=enumerate(v)
                node = LDˈn[z][g]
                load = TS.LoadYr[li0] * p # time-variant
                upResDemand = 0.03 * load
                JuMP.add_to_expression!(egp, -load)
                JuMP.add_to_expression!(egp, load, ζ[z])
                if z ∈ UpResZones
                    JuMP.add_to_expression!(rˈz, -upResDemand)
                    JuMP.add_to_expression!(rˈz, upResDemand, ζ[z])
                end
                for l=lines
                    Cnl, pfl = F[l,node], pfe[l]
                    JuMP.add_to_expression!(pfl, Cnl * -load)
                    JuMP.add_to_expression!(pfl, Cnl * load, ζ[z])
                end
            end
        end
        for (z,v)=enumerate(LDˈg)
            JuMP.add_to_expression!(Qˈy, ζCostCoeff * sum(v), ζ[z])
        end
    #= Wind =#
        for (z,v)=enumerate(Case2383.WindP) # Wind Static part: to Reserve, or to the immediate period
            JuMP.add_to_expression!(Qˈy, ϖCostCoeff * sum(v), ϖ[z])
            μfra, rˈz = TS.μ01[z], rgp[z]; ze = if z == 1
                TS.ze1
            elseif z == 2
                TS.ze2
            elseif z == 3
                TS.ze3
            end
            wfra = TS.logistic(ze[zi0]) # decided by (zone, Δi, scene)
            for (g,pMa)=enumerate(v)
                #= Reserve is done only once =#
                    upResDemand = 0.05 * μfra * pMa
                    JuMP.add_to_expression!(rˈz, -upResDemand)
                pWind, node = wfra * pMa, WN[z][g]
                JuMP.add_to_expression!(egp, pWind)
                JuMP.add_to_expression!(egp, -pWind, ϖ[z])
                for l=lines
                    Cnl, pfl = F[l,node], pfe[l]
                    JuMP.add_to_expression!(pfl, Cnl * pWind)
                    JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖ[z])
                end
            end
        end
    JuMP.@objective(m, Min, Qˈy)
    JuMP.@constraint(m, [z=UpResZones], rgp[z] ≥ 0)
    JuMP.@constraint(m, egp == 0)
    addPfLim && JuMP.@constraint(m, [l=lines], pfe[l] == pf[l])
    
    (; m, pA, ur, pI, ϱ, ϖ, ζ, pf, lines)
end

update_pAH!(pAH, m, GD) = for (z,Ø)=enumerate(GD.n), g=eachindex(Ø) pAH[z][g]=JuMP.value(m.pA[z,g]) end

end
