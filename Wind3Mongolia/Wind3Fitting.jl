# This file is to demonstrate how to simulate wind power series from wind farms in multi-locations
# For theories, refer to wind_farms_model.png

using GLMakie
import Random, Statistics, LinearAlgebra, Distributions
lin_coeff(v) = sum(v[i]v[i-1] for i=2:length(v)) / sum(v[i-1]^2 for i=2:length(v));
logit(x) = log(x / (1 - x));
logistic(x) = 1 / (1 + exp(-x));
simulation_fill!(zrnd, kv, ΣD, rng, SimL) = for i = 2:SimL
    tmpv = rand(rng, ΣD)
    tmpv .+= (kv .* zrnd[i-1, :])
    zrnd[i, :] .= tmpv
end;
function plot_1_scene!(zrnd, SEEDi, kv, ΣD, SimL, ax1, w1s, w2s, w3s, ini_index, z1, z2, z3, lw)
    rng = Random.Xoshiro(rand(Int));
    z0 = [z1[ini_index], z2[ini_index], z3[ini_index]]
    zrnd[1, :] .= z0 # initial-time
    simulation_fill!(zrnd, kv, ΣD, rng, SimL)
    zrnd[:, 1] .*= w1s.s;
    zrnd[:, 1] .+= w1s.m;
    zrnd[:, 2] .*= w2s.s;
    zrnd[:, 2] .+= w2s.m;
    zrnd[:, 3] .*= w3s.s;
    zrnd[:, 3] .+= w3s.m;
    @. zrnd = logistic(zrnd);
    # Finally, we still need 1 more step: `*P_MAX` to convert to real power
    xtks = 0:SimL-1
    lines!(ax1, xtks, zrnd[:,1], color = :red, linewidth = lw)
    lines!(ax1, xtks, zrnd[:,2], color = :blue, linewidth = lw)
    lines!(ax1, xtks, zrnd[:,3], color = :green, linewidth = lw)
end;

include("Wind3Data.jl") # load the source time series data, e.g. from Inner-Mongolia wind farms (at 3 different locations)
for k=(w1,w2,w3) k ./= maximum(k) end # result in [0,1] values
for k=(w1,w2,w3), (i,e)=enumerate(k) k[i]=clamp(e, 0.001, 0.999) end # result in (0,1) values
for k=(w1,w2,w3) @. k = logit(k) end # result in values in the \mathbb{R} axis
w1s = (m = Statistics.mean(w1), s = Statistics.std(w1; corrected=false));
w2s = (m = Statistics.mean(w2), s = Statistics.std(w2; corrected=false));
w3s = (m = Statistics.mean(w3), s = Statistics.std(w3; corrected=false));
z1 = (w1 .- w1s.m) ./ w1s.s;
z2 = (w2 .- w2s.m) ./ w2s.s;
z3 = (w3 .- w3s.m) ./ w3s.s; # Here we have derived z-score
z0 = [z1[1], z2[1], z3[1]]; # initial-time vector for the 3 wind farms
k1 = lin_coeff(z1);
k2 = lin_coeff(z2);
k3 = lin_coeff(z3); # slopes
kv = [k1, k2, k3]; # collect the slopes into a vector
Σ, mat = zeros(3,3), zeros(3,3)
for i = 2:length(z1)
    v = [z1[i] - k1 * z1[i-1], z2[i] - k2 * z2[i-1], z3[i] - k3 * z3[i-1]]
    LinearAlgebra.mul!(mat, v, v')
    Σ .+= mat
end
Σ ./= length(2:length(z1));
Σ .*= 0.1; # tuning
ΣD = Distributions.MvNormal(zeros(3), Σ); # we've got the model param (kv, ΣD)

f = Figure();
ax1 = Axis(f[1, 1])
SimL = 1 + 96
zrnd = zeros(SimL,3);
ini_index = rand(eachindex(z1))
linewidth = 0.15
for _ = 1:500
    SEEDi = rand(Int)
    plot_1_scene!(zrnd, SEEDi, kv, ΣD, SimL, ax1, w1s, w2s, w3s, ini_index, z1, z2, z3, linewidth)
end
f
