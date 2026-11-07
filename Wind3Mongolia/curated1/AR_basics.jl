# Basic ideas of the Vector AR(1) model incorporating spatial and temporal correlation among the 3 wind farms
import JLD2, Random, LinearAlgebra, Distributions
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
end;
realize_z_path!(s1, s2, s3, rng, z1, z2, z3, bvec, phi) = for i = range(2; length=length(z1)-1) # get a joint-wind series sample
    ϵv = rand(rng, D) # spatial correlations
    ϵv[1] += phi[1]z1[i-1]
    ϵv[2] += phi[2]z2[i-1]
    ϵv[3] += phi[3]z3[i-1] # temporal correlations
    ϵv .+= bvec
    s1[i], s2[i], s3[i] = ϵv
end;
getw01(w1) = w1 / maximum(w1);
open01(x) = clamp(x, 0.0001, 0.9999);
logit(x) = log(x / (1-x));
w2z = logit ∘ open01;

# Prepare Data
tmp = JLD2.load("wind3_dataset.jld2");
W1, W2, W3 = div_k(tmp["w1"]), div_k(tmp["w2"]), div_k(tmp["w3"]);
w1 = getw01(W1); w2 = getw01(W2); w3 = getw01(W3);
z1 = w2z.(w1); z2 = w2z.(w2); z3 = w2z.(w3);

using GLMakie # Visualization Settings
f = Figure();
txs = range(0; length=length(z1));
c1, c2, c3 = :deepskyblue2, :green1, :tomato1;
c4, c5, c6 = :blue, :green, :red3;
l1, l2, l3 = :dot, :dash, :solid;
mt = :xcross;
ax = Axis(f[1,1]; xlabel = "time periods", ylabel = "W/Wmax");
scatterlines!(ax, txs, z1; marker = mt,linestyle = l1,color = c1);
scatterlines!(ax, txs, z2; marker = mt,linestyle = l2,color = c2);
scatterlines!(ax, txs, z3; marker = mt,linestyle = l3,color = c3);

# AR(1) model Fitting
b1, phi1 = get_bias_phi(z1);
b2, phi2 = get_bias_phi(z2);
b3, phi3 = get_bias_phi(z3);
cov_mat = LinearAlgebra.Symmetric([
    _cov(z1,b1,phi1,z1,b1,phi1) _cov(z1,b1,phi1,z2,b2,phi2) _cov(z1,b1,phi1,z3,b3,phi3)
    0 _cov(z2,b2,phi2,z2,b2,phi2) _cov(z2,b2,phi2,z3,b3,phi3)
    0 0 _cov(z3,b3,phi3,z3,b3,phi3)
]);
D = Distributions.MvNormal(cov_mat);
phi = [phi1, phi2, phi3];
bvec = [b1, b2, b3];

# A basic sanity test
s1, s2, s3 = similar(z1), similar(z2), similar(z3);
s1[1], s2[1], s3[1] = z1[1], z2[1], z3[1]; # initial values
realize_z_path!(s1, s2, s3, Random.Xoshiro(11), z1, z2, z3, bvec, phi) # sampling 1
scatterlines!(ax, txs, s1; marker = mt,linestyle = l1,color = c4);
scatterlines!(ax, txs, s2; marker = mt,linestyle = l2,color = c5);
scatterlines!(ax, txs, s3; marker = mt,linestyle = l3,color = c6);
realize_z_path!(s1, s2, s3, Random.Xoshiro(-222), z1, z2, z3, bvec, phi) # sampling 2
scatterlines!(ax, txs, s1; marker = mt,linestyle = l1,color = c4);
scatterlines!(ax, txs, s2; marker = mt,linestyle = l2,color = c5);
scatterlines!(ax, txs, s3; marker = mt,linestyle = l3,color = c6);
f
