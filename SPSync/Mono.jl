module Mono
import ..Case2383, ..TS, ..Ops, ..Sin_Tha_Models, Random, JuMP
const UpResZones = 1:4
const ϖCostCoeff, ζ2ϖCostRatio = 2.0, 25
const urShortage2linCmaxRatio = 10

function get_a_2SSP_result(S, T, li0, li1, zi0, xH, pAH, GD, LD, liD, r1r2c1c2i1i2, LoadYr)
    mono = proto!(Ops.Model(), S, T, li0, li1, zi0, xH, pAH, GD, LD, liD, r1r2c1c2i1i2, LoadYr);
    JuMP.set_time_limit_sec(mono, 120.)
    JuMP.optimize!(mono); ter = JuMP.termination_status(mono)
    ter ∈ (JuMP.TIME_LIMIT, JuMP.OPTIMAL) || error("mono: $ter")
    GRB_lu = (JuMP.objective_bound(mono), JuMP.objective_value(mono))
    println("GRB_lu = $GRB_lu")
    GRB_lu
end

# adapted from the original
f0!(x0, GDˈn, GDˈpmax, GDˈMmRatio, GDˈpmin, GDˈUT, GDˈDT, GDˈC, GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci, pAH, xH) = for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
    PMax, MmRatio, Pmin, pAHˈzg, xHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g], xH[z][g]
    Gwid, xHˈzg1 = PMax-Pmin, xHˈzg[1]
    UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
    rui, rdi, r2s, suC, sdC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC, GDˈdci[z][g]linC
    Ru, Rd, Su, Sd = Sin_Tha_Models._RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
    if xHˈzg1 && (pAHˈzg + Pmin > Sd || any(==(false), view(xHˈzg, 1:UT)))
        JuMP.fix(x0[z,g], 1; force = true)
    elseif !xHˈzg1 && any(view(xHˈzg, 1:DT))
        JuMP.fix(x0[z,g], 0; force = true)
    end
end
f2!(x, LKA, T, GDˈn, GDˈpmax, GDˈMmRatio, GDˈpmin, pAH, xH, GDˈUT, GDˈDT, GDˈC, GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci) = for (i, Δ) = zip(LKA, 1:T)
    xm1, xNo = x[i-1], x[i]
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        PMax, MmRatio, Pmin, pAHˈzg, xHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g], xH[z][g]
        Gwid, xm1ˈzg = PMax-Pmin, xm1[z,g]
        UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
        rui, rdi, r2s, suC, sdC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC, GDˈdci[z][g]linC
        Ru, Rd, Su, Sd = Sin_Tha_Models._RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        if (JuMP.is_fixed(xm1ˈzg) && JuMP.fix_value(xm1ˈzg) == 1) && (pAHˈzg + Pmin > Sd + (Δ)Rd || any(==(false), view(xHˈzg, 1:UT-Δ)))
            JuMP.fix(xNo[z,g], 1; force = true)
        elseif (JuMP.is_fixed(xm1ˈzg) && JuMP.fix_value(xm1ˈzg) == 0) && any(view(xHˈzg, 1:DT-Δ))
            JuMP.fix(xNo[z,g], 0; force = true)
        end
    end
end

function proto!(m, S, T, li0, li1, zi0, xH, pAH, GD, LD, liD,
    r1r2c1c2i1i2, LoadYr; urShortage2linCmaxRatio=urShortage2linCmaxRatio,
    ϖCostCoeff=ϖCostCoeff, ζ2ϖCostRatio=ζ2ϖCostRatio,
    addPfLim=true)
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci, GDˈC = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.dci, GD.C
    LDˈn, LDˈg, LDˈi, F, ζCostCoeff, WN = LD.n, LD.g, LD.i, Case2383.F, (ζ2ϖCostRatio)ϖCostCoeff, Case2383.WindN
    HRZ, T = Sin_Tha_Models._hrz_T(li0, li1, T); urCshort = maximum(maximum, GDˈC)urShortage2linCmaxRatio
    LKA = range(li0+1; length=T)
    Qˈy, lines = _0(), keys(liD)
    egp,pfe,rgp = [Dict(i=>_0() for i=HRZ) for _=1:S],[Dict(i=>Dict(l=>_0() for l=lines) for i=HRZ) for _=1:S],[Dict(i=>Dict(z=>_0() for z=UpResZones) for i=HRZ) for _=1:S]
    x = [Dict(i => Dict((z,g) => JuMP.@variable(m, binary=true) for z=eachindex(GDˈn) for g=eachindex(GDˈn[z])) for i=HRZ) for _=1:S]
    #= NA =# JuMP.@constraint(m, [s=2:S, z=eachindex(GDˈn), g=eachindex(GDˈn[z])], x[s][li0][z,g] == x[1][li0][z,g])
    u7, v7, pA7, ur7, pI7, ϱ7, ua7, ϖ7, ζ7, pf7 = JuMP.@variables(m, begin
        [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
        [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
        [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
        [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
        [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
        [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
        [z=UpResZones, i=HRZ], (lower_bound = 0)
        [z=eachindex(WN),   i=HRZ], (lower_bound = 0, upper_bound = 1)
        [z=eachindex(LDˈn), i=HRZ], (lower_bound = 0, upper_bound = 1)
        [l=lines,i=HRZ], (lower_bound = -liD[l], upper_bound = liD[l])
    end)
    u, v, pA, ur, pI, ϱ, ua, ϖ, ζ, pf = map(x -> [x], (u7, v7, pA7, ur7, pI7, ϱ7, ua7, ϖ7, ζ7, pf7))
    for s = 2:S
            u7, v7, pA7, ur7, pI7, ϱ7, ua7, ϖ7, ζ7, pf7 = JuMP.@variables(m, begin
            [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
            [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
            [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
            [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
            [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
            [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=HRZ], (lower_bound = 0)
            [z=UpResZones, i=HRZ], (lower_bound = 0)
            [z=eachindex(WN),   i=HRZ], (lower_bound = 0, upper_bound = 1)
            [z=eachindex(LDˈn), i=HRZ], (lower_bound = 0, upper_bound = 1)
            [l=lines,i=HRZ], (lower_bound = -liD[l], upper_bound = liD[l])
        end)
        foreach((x, a) -> push!(x, a), (u, v, pA, ur, pI, ϱ, ua, ϖ, ζ, pf), (u7, v7, pA7, ur7, pI7, ϱ7, ua7, ϖ7, ζ7, pf7))
    end
    for s = 1:S
        f0!(x[s][li0], GDˈn, GDˈpmax, GDˈMmRatio, GDˈpmin, GDˈUT, GDˈDT, GDˈC, GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci, pAH, xH)
        f2!(x[s], LKA, T, GDˈn, GDˈpmax, GDˈMmRatio, GDˈpmin, pAH, xH, GDˈUT, GDˈDT, GDˈC, GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci)
        Sin_Tha_Models._gMdl!(m, HRZ, LKA, li0, T, u[s], v[s], x[s], pI[s], ϱ[s], pA[s], ur[s], Qˈy, egp[s], pfe[s], rgp[s], ϖCostCoeff, lines, F, GDˈn, GDˈpmax, GDˈMmRatio, GDˈpmin, pAH, xH, GDˈUT, GDˈDT, GDˈC, GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci)
        Sin_Tha_Models.mdld!(LoadYr, r1r2c1c2i1i2, egp[s], pfe[s], rgp[s], Qˈy, HRZ, LDˈg, LDˈn, LDˈi, F, ζ[s], ζCostCoeff, lines)
        Sin_Tha_Models.b_wd!(egp[s], pfe[s], rgp[s], Qˈy, Random.Xoshiro(Ψ(zi0,S,s)), HRZ, LKA, F, WN, lines, ϖCostCoeff, ϖ[s], li0, zi0) 
        Sin_Tha_Models._psys!(m, egp[s], pfe[s], rgp[s], Qˈy, lines, HRZ, pf[s], ua[s], urCshort, addPfLim)
    end
    JuMP.@objective(m, Min, Qˈy / S) # 1. doesn't contain Const, 2. doesn't relate to xLink
    m
end

_0() = JuMP.AffExpr()
Ψ(k, S, s) = (k)S + s

end

