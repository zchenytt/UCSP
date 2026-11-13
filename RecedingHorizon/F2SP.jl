module F2SP
import Random, JuMP, ..Case2383, ..TS
const ı = JuMP.value
const UpResZones = 1:4
const ϖCostCoeff = 2.0
const ζ2ϖCostRatio = 25
const urShortage2linCmaxRatio = 10

rB(x) = round(Bool, x)
_0() = JuMP.AffExpr()
function _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
    Ru, Rd = Gwid/rui, Gwid/rdi
    Su, Sd = clamp(r2s * Ru, Pmin, PMax), clamp(r2s * Rd, Pmin, PMax)
    Ru, Rd, Su, Sd
end
function _hrz_T(li0, li1, T)
    if li1 - li0 < T
        HRZ = li0:li1
        T = length(HRZ)-1
    else
        HRZ = range(li0; length=1+T)
    end
    HRZ, T
end
_Dis(HRZ, S1S) = Dict((i,s) => _0() for i=HRZ for s=S1S)

"""
 m Presolve 2
 wind look-ahead scenes are completely decided by zi0 + rngWindScenes
"""
function Model!(m, li0, li1, zi0, T, S, xH, pAH, GD, LD, liD;
    urShortage2linCmaxRatio=urShortage2linCmaxRatio, ϖCostCoeff=ϖCostCoeff, ζ2ϖCostRatio=ζ2ϖCostRatio, PfK=1.0, addPfLim=true)
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci, GDˈC = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.dci, GD.C
    LDˈn, LDˈg, F, ζCostCoeff, WN = LD.n, LD.g, Case2383.F, (ζ2ϖCostRatio)ϖCostCoeff, Case2383.WindN
    HRZ, T = _hrz_T(li0, li1, T); LKA, S1S = HRZ[2:end], 1:S
    urCshort = urShortage2linCmaxRatio * maximum(maximum, GDˈC) # planning horizon
    rngW, Qˈy, uvCostNow, lines = Random.Xoshiro(zi0 + S), _0(), _0(), keys(liD)
    egp, pfe, rgp = _Dis(HRZ,S1S), Dict(l => _Dis(HRZ,S1S) for l=lines), [_Dis(HRZ,S1S) for z=UpResZones]
    JuMP.@variables(m, begin
            x[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ, s=S1S], Bin
        0 ≤ u[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ, s=S1S] ≤ 1
        0 ≤ v[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ, s=S1S] ≤ 1
        0 ≤ pA[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ, s=S1S]
        0 ≤ ur[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ, s=S1S]
        0 ≤ pI[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ, s=S1S] # p_Injection into the network
        0 ≤  ϱ[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ, s=S1S] # curtailed generation
        0 ≤ ua[z=UpResZones, i=HRZ, s=S1S]
        0           ≤   ϖ[z=eachindex(WN),   i=HRZ, s=S1S] ≤ 1
        0           ≤   ζ[z=eachindex(LDˈn), i=HRZ, s=S1S] ≤ 1
        -liD[l]PfK  ≤  pf[l=lines,i=HRZ,s=S1S]             ≤ liD[l]PfK
    end)
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
            PMax, MmRatio, Pmin, pAHˈzg, xHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g], xH[z][g]
            Gwid, xHˈzg1 = PMax-Pmin, xHˈzg[1]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC, sdC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC, GDˈdci[z][g]linC
            Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        #= NA =#
            JuMP.@constraints(m, begin
                [i=HRZ, s=2:S], x[z,g,i,s] == x[z,g,i,1]
                [i=HRZ, s=2:S], u[z,g,i,s] == u[z,g,i,1]
                [i=HRZ, s=2:S], v[z,g,i,s] == v[z,g,i,1]
            end)
        #= UVX-Poly =#
            for s=S1S, i=HRZ
                e = _0()
                JuMP.add_to_expression!(e, u[z,g,i,s])
                JuMP.add_to_expression!(e, -1., v[z,g,i,s])
                JuMP.add_to_expression!(e, -1., x[z,g,i,s])
                if i==li0
                    JuMP.add_to_expression!(e, xHˈzg1)
                else
                    JuMP.add_to_expression!(e, x[z,g,i-1,s])
                end
                JuMP.@constraint(m, e == 0)
                eu, ev, CnT = _0(), _0(), 0
                for j = range(i; step = -1, length = max(UT, DT))
                    if j ≥ li0
                        CnT < UT && JuMP.add_to_expression!(eu, u[z,g,j,s])
                        CnT < DT && JuMP.add_to_expression!(ev, -1., v[z,g,j,s])
                    else
                        ȷ = li0 - j # > 0
                        Xnr, Xfr = rB(xHˈzg[ȷ]), rB(xHˈzg[ȷ+1])
                        if Xnr == Xfr
                            nothing
                        elseif Xnr
                            CnT < UT && JuMP.add_to_expression!(eu, 1.)
                        else # Xnr = false, and Xfr = true
                            CnT < DT && JuMP.add_to_expression!(ev, -1.)
                        end
                    end
                    CnT += 1
                end
                JuMP.@constraint(m, eu ≤ x[z,g,i,s])
                JuMP.@constraint(m, x[z,g,i,s] - 1 ≤ ev)
            end
        #= Defining Injection Power =#
            JuMP.@constraint(m, [i=HRZ, s=S1S], pI[z,g,i,s] + ϱ[z,g,i,s] == x[z,g,i,s]Pmin + pA[z,g,i,s])
        #= Costs =#
            for i=HRZ, s=S1S
                JuMP.add_to_expression!(Qˈy, suC / S, u[z,g,i,s]) # Start-up
                JuMP.add_to_expression!(Qˈy, sdC / S, v[z,g,i,s]) # Shut-dn
                JuMP.add_to_expression!(Qˈy, linC / S, pI[z,g,i,s]) # Variable Generation
                JuMP.add_to_expression!(Qˈy, linC / S,  ϱ[z,g,i,s]) # Variable Generation
                JuMP.add_to_expression!(Qˈy, ϖCostCoeff / S, ϱ[z,g,i,s]) # Generation Curtailment
            end
            JuMP.add_to_expression!(uvCostNow, suC, u[z,g,li0,1])
            JuMP.add_to_expression!(uvCostNow, sdC, v[z,g,li0,1])
        #= Reserve =#
            1≤z≤4 && for i=HRZ, s=S1S JuMP.add_to_expression!(rgp[z][i,s], ur[z,g,i,s]) end
        #= egp / pfl =#
            for i=HRZ, s=S1S JuMP.add_to_expression!(egp[i,s], pI[z,g,i,s]) end
            for l=lines
                Cnl, pfl = F[l,node], pfe[l]
                for i=HRZ, s=S1S JuMP.add_to_expression!(pfl[i,s], Cnl, pI[z,g,i,s]) end
            end
        #= Ramp / Output_Ub =#
            if MmRatio > 1.1         
                Ju, Jv = Su - Pmin, Sd - Pmin - Rd
                i = li0; JuMP.@constraint(m, [s=S1S], pA[z,g,i,s] + ur[z,g,i,s] - pAHˈzg ≤ Ru * xHˈzg1 + Ju * u[z,g,i,s])
                JuMP.@constraint(m, [i=range(li0+1; length=T), s=S1S], 
                    pA[z,g,i,s] + ur[z,g,i,s] - pA[z,g,i-1,s] ≤ Ru * x[z,g,i-1,s] + Ju * u[z,g,i,s]
                )
                i = li0; JuMP.@constraint(m, [s=S1S], pAHˈzg - pA[z,g,i,s] ≤ Rd * xHˈzg1 + Jv * v[z,g,i,s])
                JuMP.@constraint(m, [i=range(li0+1; length=T), s=S1S],
                    pA[z,g,i-1,s] - pA[z,g,i,s] ≤ Rd * x[z,g,i-1,s] + Jv * v[z,g,i,s]
                )
            end
            JuMP.@constraint(m, [i=HRZ, s=S1S], pA[z,g,i,s] + ur[z,g,i,s] ≤ Gwid * x[z,g,i,s])
    end
    #= Load =#
        for (z,v)=enumerate(LDˈg)
            1≤z≤4 && (rˈz = rgp[z]) 
            for (g,p)=enumerate(v)
                node = LDˈn[z][g]
                for i=HRZ
                    load = TS.LoadYr[i] * p # time-variant
                    upResDemand = 0.03 * load
                    for s=S1S
                        JuMP.add_to_expression!(egp[i,s], -load)
                        JuMP.add_to_expression!(egp[i,s], load, ζ[z,i,s])
                        if 1≤z≤4
                            JuMP.add_to_expression!(rˈz[i,s], -upResDemand)
                            JuMP.add_to_expression!(rˈz[i,s], upResDemand, ζ[z,i,s])
                        end
                        for l=lines
                            Cnl, pfl = F[l,node], pfe[l][i,s]
                            JuMP.add_to_expression!(pfl, Cnl * -load)
                            JuMP.add_to_expression!(pfl, Cnl * load, ζ[z,i,s])
                        end
                    end
                end
            end
        end
        for s=S1S, i=HRZ, (z,v)=enumerate(LDˈg)
            JuMP.add_to_expression!(Qˈy, ζCostCoeff * sum(v)/S, ζ[z,i,s])
        end
    #= Wind =#
        Δi = 0; for (z,v)=enumerate(Case2383.WindP) # Wind Static part: to Reserve, or to the immediate period
            for s=S1S JuMP.add_to_expression!(Qˈy, ϖCostCoeff * sum(v)/S, ϖ[z, li0+Δi, s]) end
            μfra, rˈz = TS.μ01[z], rgp[z]; ze = if z == 1
                TS.ze1
            elseif z == 2
                TS.ze2
            elseif z == 3
                TS.ze3
            end
            wfra = TS.logistic(ze[zi0 + Δi]) # decided by (zone, Δi, scene)
            for (g,pMa)=enumerate(v)
                #= Reserve is done only once =#
                    upResDemand = 0.05 * μfra * pMa
                    for i=HRZ, s=S1S JuMP.add_to_expression!(rˈz[i,s], -upResDemand) end
                i, pWind, node = li0 + Δi, wfra * pMa, WN[z][g]
                for s=S1S
                    JuMP.add_to_expression!(egp[i,s], pWind)
                    JuMP.add_to_expression!(egp[i,s], -pWind, ϖ[z,i,s])
                    for l=lines
                        Cnl, pfl = F[l,node], pfe[l][i,s]
                        JuMP.add_to_expression!(pfl, Cnl * pWind)
                        JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖ[z,i,s])
                    end
                end
            end
        end
        for s=S1S
            z1, z2, z3 = map(v -> v[zi0], (TS.ze1, TS.ze2, TS.ze3)) # reload initial condition
            for i = li0+1:li0+T
                egpˈis, ϵ = egp[i,s], rand(rngW, TS.D)
                z1, z2, z3 = z123 = TS.ar(ϵ,1,z1), TS.ar(ϵ,2,z2), TS.ar(ϵ,3,z3)
                WFRA = map(TS.logistic, z123)
                for (z,v)=enumerate(Case2383.WindP)
                    wfra, ϖˈzis = WFRA[z], ϖ[z,i,s]
                    JuMP.add_to_expression!(Qˈy, ϖCostCoeff * sum(v)/S, ϖˈzis)
                    for (g,pMa)=enumerate(v)
                        pWind = wfra * pMa
                        JuMP.add_to_expression!(egpˈis, pWind)
                        JuMP.add_to_expression!(egpˈis, -pWind, ϖˈzis)
                        node = WN[z][g]
                        for l=lines
                            Cnl, pfl = F[l,node], pfe[l][i,s]
                            JuMP.add_to_expression!(pfl, Cnl * pWind)
                            JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖˈzis)
                        end
                    end
                end
            end
        end
    #= buy up reserve =#
        for i=ua JuMP.add_to_expression!(Qˈy, urCshort/S, i) end
    JuMP.@objective(m, Min, Qˈy)
    JuMP.@constraint(m, [z=UpResZones, i=HRZ, s=S1S], rgp[z][i,s] + ua[z,i,s] ≥ 0)
    JuMP.@constraint(m, [i=HRZ, s=S1S], egp[i,s] == 0)
    addPfLim && JuMP.@constraint(m, [l=lines, i=HRZ, s=S1S], pfe[l][i,s] == pf[l,i,s])
    
    (; m, x, u, v, pA, ur, pI, ϱ, ϖ, ζ, pf, HRZ, S1S, lines, uvCostNow)
end

#= Receding horizon related =#
uv_immediate_cost(m) = ı(m.uvCostNow)
update_xH_vec!(m, xH, GD) = for (z,Ø)=enumerate(GD.n), g=eachindex(Ø)
    update_xH_vec!(m, xH, z, g)
end
function update_xH_vec!(m, xH, z::Int, g)
    v, x = xH[z][g], m.x
    i, s, _ = first(m.HRZ), 1, circshift!(v, 1)
    v[1] = round(Bool, ı(x[z,g,i,s]))
end

end
