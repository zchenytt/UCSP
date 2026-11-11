"Time Series"
module TS
import JLD2, Random, LinearAlgebra, Distributions
const Wind3DIR = "jld2data/wind3_dataset.jld2"
const LoadYrDIR = "jld2data/hourly_load.jld2"
logit(F) = log(F / (1-F))
logistic(R) = 1/(1 + exp(-R))

_load01(s) = (load_MW=JLD2.load(s)["v"]; load_MW ./= maximum(load_MW))
function _z1z2z3(s)
    d = JLD2.load(s)
    W1, W2, W3 = div_k(d["w1"]), div_k(d["w2"]), div_k(d["w3"])
    E1, E2, E3 = div_k(d["e1"]), div_k(d["e2"]), div_k(d["e3"])
    w1, w2, w3 = getw01(W1), getw01(W2), getw01(W3)
    e1, e2, e3 = getw01(E1), getw01(E2), getw01(E3)
    w2z.(w1), w2z.(w2), w2z.(w3), w2z.(e1), w2z.(e2), w2z.(e3)
end
function div_k(w; k=4) # convert 15-min wind series to hourly ones
    length(w) % k == 0 || error()
    v, s = Float64[], 0.
    for (i,e)=enumerate(w)
        s += e
        if i % k == 0
            push!(v, s/k)
            s = 0.
        end
    end
    v
end
getw01(w1) = w1 / maximum(w1);
open01(x) = clamp(x, 0.0001, 0.9999);
w2z = logit ∘ open01;
"z[t+1] = bias + (phi)z[t] + eps_t = (1-phi)mean + (phi)z[t] + eps_t; where, mean is the long-term steady value"
function get_bias_phi(z) # Least Squares Estimation
    S = length(z) - 1
    r = range(1; length=S)
    xbar, ybar = sum(z[i] for i=r)/S, sum(z[i+1] for i=r)/S
    phi = sum((z[i] - xbar) * (z[i+1] - ybar) for i=r) / sum((z[i] - xbar)^2 for i=r)
    bias = ybar - (phi)xbar
    bias, phi
end;
function _cov(v, bv, phiv, u, bu, phiu) # cov(X, Y)::Float64
    S = length(v) - 1
    r = range(1; length=S)
    sum((bv + v[t]phiv - v[t+1]) * (bu + u[t]phiu - u[t+1]) for t=r) / (S - 2)
end;
realize_z_path!(s1, s2, s3, rng, z1, z2, z3, bvec, phi) = for i = range(2; length=length(z1)-1) # get a joint-wind series sample
    ϵv = rand(rng, D) # spatial correlations
    ϵv[1] += phi[1]z1[i-1]
    ϵv[2] += phi[2]z2[i-1]
    ϵv[3] += phi[3]z3[i-1] # temporal correlations
    ϵv .+= bvec
    s1[i], s2[i], s3[i] = ϵv
end;
function get_wind_params(z1, z2, z3)
    b1, phi1 = get_bias_phi(z1)
    b2, phi2 = get_bias_phi(z2)
    b3, phi3 = get_bias_phi(z3)
    cov_mat = LinearAlgebra.Symmetric([
        _cov(z1,b1,phi1,z1,b1,phi1) _cov(z1,b1,phi1,z2,b2,phi2) _cov(z1,b1,phi1,z3,b3,phi3)
        0 _cov(z2,b2,phi2,z2,b2,phi2) _cov(z2,b2,phi2,z3,b3,phi3)
        0 0 _cov(z3,b3,phi3,z3,b3,phi3)
    ])
    phi = [phi1, phi2, phi3]
    bvec = [b1, b2, b3]
    φ, b, D = phi, bvec, Distributions.MvNormal(cov_mat)
end

const LoadYr = _load01(LoadYrDIR) # 1-year time series
const z1, z2, z3, ze1, ze2, ze3 = _z1z2z3(Wind3DIR) # z's (real-number) series; We need these since AR is defined on it!
const φ, b, D = get_wind_params(z1, z2, z3) # AR(1) param on `z`
const μ01 = @. logistic(b / (1 - φ)) # (timeless mean windgen) will be used in Reserve Data input calculation
def_t0_indices(rng) = (rand(rng, 1:8760), 1) # li0, zi0 (zei0)

ar(ϵv, zone, z) = b[zone] + φ[zone]z + ϵv[zone]


end
