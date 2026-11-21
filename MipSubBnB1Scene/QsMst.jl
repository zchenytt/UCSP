"""
a function to generate a master at any BnB node (allocate a new one for a new node)

a function to generate a prototype subproblem (allocate only once)
"""
module QsMst
import Random, JuMP, ..Case2383, ..TS, ..Ops
const UpResZones = 1:4
const ϖCostCoeff, ζ2ϖCostRatio = 2.0, 25
const urShortage2linCmaxRatio = 10

_0() = JuMP.AffExpr()
_d(N) = Vector{Cdouble}(undef, N)
_i(N) = Vector{Cint}(undef, N)
function _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
    Ru, Rd = Gwid/rui, Gwid/rdi
    Su, Sd = clamp(r2s * Ru, Pmin, PMax), clamp(r2s * Rd, Pmin, PMax)
    Ru, Rd, Su, Sd
end
function _hrz_T(li0, li1, T) # The entry `T` is the original desired (fixed) param
    if li1 - li0 < T
        HRZ = li0:li1
        T = length(HRZ)-1
    else
        HRZ = range(li0; length=1+T)
    end
    HRZ, T # This `T` is the really-used num of look-ahead periods
end
train_tp(s, mlp, mcn, a...) = (proto!(s,mlp,a...;binary=false),proto!(s,mcn,a...;binary=true))
Ψ(k, S, s) = (k)S + s
function proto!(s, m, li0, li1, zi0, T, S, xH, pAH, GD, LD, liD;
    binary,
    urShortage2linCmaxRatio=urShortage2linCmaxRatio,
    ϖCostCoeff=ϖCostCoeff,
    ζ2ϖCostRatio=ζ2ϖCostRatio,
    PfK=1.0,
    addPfLim=true)
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci, GDˈC = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.dci, GD.C
    LDˈn, LDˈg, F, ζCostCoeff, WN = LD.n, LD.g, Case2383.F, (ζ2ϖCostRatio)ϖCostCoeff, Case2383.WindN
    HRZ, T = _hrz_T(li0, li1, T); urCshort = maximum(maximum, GDˈC)urShortage2linCmaxRatio
    LKA = range(li0+1; length=T)
    rngW, Qˈy, lines = Random.Xoshiro(Ψ(zi0,S,s)), _0(), keys(liD)
    egp,pfe,rgp = Dict(i=>_0() for i=HRZ),Dict(i=>Dict(l=>_0() for l=lines) for i=HRZ),Dict(i=>Dict(z=>_0() for z=UpResZones) for i=HRZ)
    #= The x = [xLink; other_x_01Vars] =#
        x = Dict(i => Dict{Tuple{Int,Int}, Union{Bool,JuMP.VariableRef}}() for i=HRZ)
        i = li0; x0 = x[i]; N = 0; for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
            PMax, MmRatio, Pmin, pAHˈzg, xHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g], xH[z][g]
            Gwid, xHˈzg1 = PMax-Pmin, xHˈzg[1]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC, sdC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC, GDˈdci[z][g]linC
            Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
            if xHˈzg1 && (pAHˈzg + Pmin > Sd || any(==(false), view(xHˈzg, 1:UT)))
                x0[z,g] = true
            elseif !xHˈzg1 && any(view(xHˈzg, 1:DT))
                x0[z,g] = false
            else
                x0[z,g], N = JuMP.@variable(m, binary=binary), N+1
            end
        end; Xl,rC,Ci,Cd=_d(N),_d(N),_i(N+2),_d(N+2)
        Ο = 0; for (i, Δ) = zip(LKA, 1:T)
            xm1, xNo = x[i-1], x[i]
            for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
                PMax, MmRatio, Pmin, pAHˈzg, xHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g], xH[z][g]
                Gwid, xm1ˈzg = PMax-Pmin, xm1[z,g]
                UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
                rui, rdi, r2s, suC, sdC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC, GDˈdci[z][g]linC
                Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
                if xm1ˈzg === true && (pAHˈzg + Pmin > Sd + (Δ)Rd || any(==(false), view(xHˈzg, 1:UT-Δ)))
                    xNo[z,g] = true
                elseif xm1ˈzg === false && any(view(xHˈzg, 1:DT-Δ))
                    xNo[z,g] = false
                else
                    xNo[z,g], Ο = JuMP.@variable(m, binary=binary), Ο+1
                end
            end
        end
    JuMP.@variables(m, begin
        0 ≤ u[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ]
        0 ≤ v[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ]
        0 ≤ pA[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ]
        0 ≤ ur[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ]
        0 ≤ pI[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ] # p_Injection into the network
        0 ≤  ϱ[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ] # curtailed generation
        0 ≤ ua[z=UpResZones, i=HRZ]
        0           ≤   ϖ[z=eachindex(WN),   i=HRZ] ≤ 1
        0           ≤   ζ[z=eachindex(LDˈn), i=HRZ] ≤ 1
        -liD[l]PfK  ≤  pf[l=lines,i=HRZ]             ≤ liD[l]PfK
    end)
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
            PMax, MmRatio, Pmin, pAHˈzg, xHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g], xH[z][g]
            Gwid, xHˈzg1 = PMax-Pmin, xHˈzg[1]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC, sdC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC, GDˈdci[z][g]linC
            Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        #= UVX-Poly =#
            for i=HRZ
                JuMP.@constraint(m, u[z,g,i]-v[z,g,i] == x[i][z,g] - (i>li0 ? x[i-1][z,g] : xHˈzg1))
                eu, ev, CnT = _0(), _0(), 0
                for j = range(i; step = -1, length = max(UT, DT))
                    if j ≥ li0
                        CnT < UT && JuMP.add_to_expression!(eu, u[z,g,j])
                        CnT < DT && JuMP.add_to_expression!(ev, v[z,g,j])
                    else
                        ȷ = li0 - j # > 0
                        Xnr, Xfr = xHˈzg[ȷ], xHˈzg[ȷ+1]
                        if Xnr == Xfr
                            nothing # add 0. implicitly
                        elseif Xnr # turn on
                            CnT < UT && JuMP.add_to_expression!(eu, 1.)
                        else # turn off
                            CnT < DT && JuMP.add_to_expression!(ev, 1.)
                        end
                    end
                    CnT += 1
                end
                JuMP.@constraint(m, eu ≤ x[i][z,g]) # => x must be positive
                JuMP.@constraint(m, x[i][z,g] ≤ 1 - ev) # => x must be in [0,1], and therefore also make u,v in [0,1]
            end
        #= Defining Injection Power =#
            JuMP.@constraint(m, [i=HRZ], pI[z,g,i] + ϱ[z,g,i] == x[i][z,g]Pmin + pA[z,g,i])
        #= Costs =#
            for i=HRZ
                JuMP.add_to_expression!(Qˈy, suC, u[z,g,i]) # Start-up
                JuMP.add_to_expression!(Qˈy, sdC, v[z,g,i]) # Shut-dn
                JuMP.add_to_expression!(Qˈy, linC, pI[z,g,i]) # Variable Generation
                JuMP.add_to_expression!(Qˈy, linC,  ϱ[z,g,i]) # Variable Generation
                JuMP.add_to_expression!(Qˈy, ϖCostCoeff, ϱ[z,g,i]) # Generation Curtailment
            end
        #= Reserve =#
            z ∈ UpResZones && for i=HRZ JuMP.add_to_expression!(rgp[i][z], ur[z,g,i]) end
        #= egp / pfl =#
            for i = HRZ
                pfl, pIzgi = pfe[i], pI[z,g,i]
                JuMP.add_to_expression!(egp[i], pIzgi)
                for l = lines
                    JuMP.add_to_expression!(pfl[l], F[l,node], pIzgi)
                end
            end
        #= Ramp / Output_Ub =#
            if MmRatio > 1.1
                Ju, Jv = Su - Pmin, Sd - Pmin - Rd; 
                i = li0; JuMP.@constraint(m, pA[z,g,i] + ur[z,g,i] - pAHˈzg      ≤ Ru * xHˈzg1      + Ju * u[z,g,i])
                JuMP.@constraint(m, [i=LKA], pA[z,g,i] + ur[z,g,i] - pA[z,g,i-1] ≤ Ru * x[i-1][z,g] + Ju * u[z,g,i])
                i = li0; JuMP.@constraint(m, pAHˈzg      - pA[z,g,i] ≤ Rd * xHˈzg1     + Jv * v[z,g,i])
                JuMP.@constraint(m, [i=LKA], pA[z,g,i-1] - pA[z,g,i] ≤ Rd * x[i-1][z,g] + Jv * v[z,g,i])
                if UT > 1
                    JuMP.@constraint(m, [i=HRZ], pA[z,g,i]+ur[z,g,i] ≤ x[i][z,g]Gwid-(PMax - Su)u[z,g,i]-(i==li0+T ? 0 : (PMax-Sd)v[z,g,i+1]))
                else
                    JuMP.@constraint(m, [i=HRZ], pA[z,g,i]+ur[z,g,i] ≤ x[i][z,g]Gwid-max(Sd-Su,0.)u[z,g,i]-(i==li0+T ? 0 : (PMax-Sd)v[z,g,i+1]))
                    JuMP.@constraint(m, [i=HRZ], pA[z,g,i]+ur[z,g,i] ≤ x[i][z,g]Gwid-(PMax - Su)u[z,g,i]-(i==li0+T ? 0 : max(Su-Sd,0.)v[z,g,i+1]))
                end
            else
                JuMP.@constraint(m, [i=HRZ], pA[z,g,i]+ur[z,g,i] ≤ x[i][z,g]Gwid)
            end
    end
    #= Load =#
        for i=HRZ
            rgpˈi, egpˈi, pfeˈi, Ld01ˈi = rgp[i], egp[i], pfe[i], TS.LoadYr[i]
            for (z,v)=enumerate(LDˈg)
                ζˈzi = ζ[z,i]; JuMP.add_to_expression!(Qˈy, ζCostCoeff * sum(v), ζˈzi)
                z ∈ UpResZones && (rgpˈiz = rgpˈi[z])
                for (g,p)=enumerate(v)
                    node = LDˈn[z][g]
                    load = (p)Ld01ˈi # time-variant
                    upResDemand = 0.03 * load
                    JuMP.add_to_expression!(egpˈi, -load)
                    JuMP.add_to_expression!(egpˈi, load, ζˈzi)
                    z ∈ UpResZones && JuMP.add_to_expression!(rgpˈiz, -upResDemand)
                    z ∈ UpResZones && JuMP.add_to_expression!(rgpˈiz, upResDemand, ζˈzi)
                    for l=lines
                        Cnl, pfl = F[l,node], pfeˈi[l]
                        JuMP.add_to_expression!(pfl, Cnl * -load)
                        JuMP.add_to_expression!(pfl, Cnl * load, ζˈzi)
                    end
                end
            end
        end
    #= Wind =#
        Δi = 0; for (z,v)=enumerate(Case2383.WindP) # Wind Static part: to Reserve, or to the immediate period
            JuMP.add_to_expression!(Qˈy, ϖCostCoeff * sum(v), ϖ[z, li0+Δi])
            μfra = TS.μ01[z]; ze = if z == 1
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
                    for i=HRZ JuMP.add_to_expression!(rgp[i][z], -upResDemand) end
                i, pWind, node = li0+Δi, wfra * pMa, WN[z][g]
                JuMP.add_to_expression!(egp[i], pWind)
                JuMP.add_to_expression!(egp[i], -pWind, ϖ[z,i])
                for l=lines
                    Cnl, pfl = F[l,node], pfe[i][l]
                    JuMP.add_to_expression!(pfl, Cnl * pWind)
                    JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖ[z,i])
                end
            end
        end
        z1, z2, z3 = map(v -> v[zi0], (TS.ze1, TS.ze2, TS.ze3)) # reload initial condition
        for i = LKA
            egpˈi, ϵ = egp[i], rand(rngW, TS.D)
            z1, z2, z3 = z123 = TS.ar(ϵ,1,z1), TS.ar(ϵ,2,z2), TS.ar(ϵ,3,z3)
            WFRA = map(TS.logistic, z123)
            for (z,v)=enumerate(Case2383.WindP)
                wfra, ϖˈzi = WFRA[z], ϖ[z,i]
                JuMP.add_to_expression!(Qˈy, ϖCostCoeff * sum(v), ϖˈzi)
                for (g,pMa)=enumerate(v)
                    pWind = wfra * pMa
                    JuMP.add_to_expression!(egpˈi, pWind)
                    JuMP.add_to_expression!(egpˈi, -pWind, ϖˈzi)
                    node = WN[z][g]
                    for l=lines
                        Cnl, pfl = F[l,node], pfe[i][l]
                        JuMP.add_to_expression!(pfl, Cnl * pWind)
                        JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖˈzi)
                    end
                end
            end
        end
    #= buy up reserve =#
        for i=ua JuMP.add_to_expression!(Qˈy, urCshort, i) end
    JuMP.@objective(m, Min, Qˈy) # 1. doesn't contain Const, 2. doesn't relate to xLink
    JuMP.@constraint(m, [z=UpResZones, i=HRZ], rgp[i][z] + ua[z,i] ≥ 0)
    JuMP.@constraint(m, [i=HRZ], egp[i] == 0)
    addPfLim && JuMP.@constraint(m, [l=lines, i=HRZ], pfe[i][l] == pf[l,i])
    for (ji, gi) = zip(eachindex(Ci), range(S; length=N))
        Ci[ji] = gi
    end; Ci[N+1] = Ops.ni2gi4s(s)
    (; s, m, N, x, Xl, rC, HRZ, T, Ci, Cd)
end
function mst!(m, li0, zi0, S, xH, pAH, GD, LD, liD; θ_has_obj)
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci, GDˈC = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.dci, GD.C
    LDˈn, LDˈg, F, ζCostCoeff, WN = LD.n, LD.g, Case2383.F, (ζ2ϖCostRatio)ϖCostCoeff, Case2383.WindN
    JuMP.@variable(m, θ[1:S])
        x0,N=Dict{Tuple{Int,Int},Union{Bool,JuMP.VariableRef}}(),0
        for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
            PMax, MmRatio, Pmin, pAHˈzg, xHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g], xH[z][g]
            Gwid, xHˈzg1 = PMax-Pmin, xHˈzg[1]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC, sdC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC, GDˈdci[z][g]linC
            Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
            if xHˈzg1 && (pAHˈzg + Pmin > Sd || any(==(false), view(xHˈzg, 1:UT)))
                x0[z,g] = true
            elseif !xHˈzg1 && any(view(xHˈzg, 1:DT))
                x0[z,g] = false
            else
                x0[z,g], N = JuMP.@variable(m), N+1
            end
        end; Bv=fill(Cchar('C'),N)
    chk = (X=fill(NaN,N), Θ=fill(NaN,S), ul=fill(Inf,S+1)) # the "check" solution, where (X,Θ) are loaded directly
    # while ul is an additional storage space, ul[1:S] is Qy_values evaluated, Qˈy[S+1] is LB-of-mst (ObjBound)
    # ast = (X=fill(NaN,N), ul=fill(Inf,S+1)) # store the current-best primal feasible solution
    # printstyled("master built with N=$N, S=$S, θ_has_obj=$θ_has_obj\n"; color = 30)
    JuMP.@objective(m, Min, sum(θ_has_obj/S * i for i = θ))
    (; m, N, S, Bv, chk)
end

end
