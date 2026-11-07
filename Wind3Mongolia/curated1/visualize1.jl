getw01(w1,e1) = (t=max(maximum(w1), maximum(e1)); (w1/t, e1/t));
import JLD2
d = JLD2.load("wind3_dataset.jld2");
w1, w2, w3, e1, e2, e3 = d["w1"], d["w2"], d["w3"], d["e1"], d["e2"], d["e3"];
using GLMakie
f = Figure();
lt, le = length(w1), length(e1);
w1, e1 = getw01(w1,e1);
w2, e2 = getw01(w2,e2);
w3, e3 = getw01(w3,e3);
txs, exs = range(0; length=lt), range(lt; length=le)
c1, c2, c3 = :deepskyblue2, :green1, :tomato1
c4, c5, c6 = :blue, :green, :red3
l1, l2, l3 = :dot, :dash, :solid
mt, me = :star4, :xcross
ax = Axis(f[1,1]; xticks = ([0,95,671,lt,lt+le-1],["Sept 25th 00:00","day","week","Oct 25th 00:00","Oct 29th 23:45"]), 
title = "15min-resolution joint wind generation of 3 wind farms (QiFeng/XiMeng/QingShuiHe) in Inner Mongolia in 2020; split into fit/test intervals",
xlabel = "time periods", ylabel = "W/Wmax");
scatterlines!(ax, txs, w1; marker = mt,linestyle = l1,color = c1);
scatterlines!(ax, txs, w2; marker = mt,linestyle = l2,color = c2);
scatterlines!(ax, txs, w3; marker = mt,linestyle = l3,color = c3);
scatterlines!(ax, exs, e1; marker = me,linestyle = l1,color = c4);
scatterlines!(ax, exs, e2; marker = me,linestyle = l2,color = c5);
scatterlines!(ax, exs, e3; marker = me,linestyle = l3,color = c6);
f
