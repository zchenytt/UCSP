# residential, commercial, industrial: 2 curves per type
import JLD2
Dct = JLD2.load("r1r2c1c2i1i2.jld2");
(r1, r2, c1, c2, i1, i2) = first(Dct).second;
using GLMakie
f = Figure();
ax = Axis(f[1,1]);
ax2 = Axis(f[2,1]);
ts = range(0; length = length(r1));
scatterlines!(ax, ts, r1; linestyle = :dot, color = :blue, label = "R1");
scatterlines!(ax, ts, c1; linestyle = :dash, color = :green, label = "C1");
scatterlines!(ax, ts, i1; linestyle = :solid, color = :tomato, label = "I1");
scatterlines!(ax2, ts, r2; linestyle = :dot, color = :blue, label = "R2");
scatterlines!(ax2, ts, c2; linestyle = :dash, color = :green, label = "C2");
scatterlines!(ax2, ts, i2; linestyle = :solid, color = :tomato, label = "I2");
axislegend(ax, position = :lb)
axislegend(ax2, position = :lb)
f
