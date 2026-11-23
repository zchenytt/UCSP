module Bnb
import ..Ops, JuMP

function select_i_X(Xres, NDˈrˈb)
    (lbv, ubv), (ib, Xb, dmin), Nf = NDˈrˈb, (0, 0., Inf), 0
    for (i, x)=enumerate(Xres)
        (lbv[i] != ubv[i] && Ops.isfrac(x)) || continue
        d = abs(x - .5)
        d < dmin && ((ib, Xb, dmin) = (i, x, d))
        Nf += 1
    end
    ib, Xb, Nf # if ib==0, Xres should be Int-Feasible
end
function branch_and_closeˈnonleaf_Nnow!(ND, Nnow, iBr, vBr, use_sB)
    b = round(Bool, vBr)
    N1 = new_a_node!(ND, Nnow, iBr,  b, use_sB) # small label to more promising branch
    N2 = new_a_node!(ND, Nnow, iBr, !b, use_sB) # the other possibility
    ND.o[Nnow] = ND.e[Nnow] = false
    N1, N2
end
function new_a_node!(ND, Nnow, iBr, b, use_sB)
    N = 1+Ops.currentMaxN(ND.nc) # Label a new node
    ND.r[N] = tRef = Ops.BnBNodeRef(Nnow,
        map(copy, ND.r[Nnow].b),
        map(copy, ND.r[Nnow].c)
    )
    tRef.b[2-b][iBr] = b # perform a new branching step
    ND.o[N] = ND.e[N] = true
    ND.l[N], ND.s[N] = ND.l[Nnow], use_sB & ND.s[Nnow]
    k = Ops.nXfixed(tRef.b)
    haskey(ND.nc, k) ? push!(ND.nc[k], N) : (ND.nc[k] = [N])
    N
end
function pick_an_open_node(ND, srchMod)
    f, (; nc, o) = ifelse(srchMod, reverse, identity), ND
    for l=f(0:maximum(keys(nc))), n=sort!(nc[l])
        o[n] && return n
    end;  return Ops.N_Dumb
end
function log_gaps(lb, ub; color = 52)
    agap = ub - lb
    rgap = agap / ub
    printstyled("lb= $lb, ub= $ub, agap= $agap, rgap= $rgap\n"; color = color)
    agap, rgap
end

end
