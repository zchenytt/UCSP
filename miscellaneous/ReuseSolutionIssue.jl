# Whenever utilizing a solution from an upstream solve as the data input
# to a downstream solve, you need some proper conversion on the raw result
# retrieved via JuMP.value
# Here's an example

import JuMP, Gurobi
m = JuMP.read_from_file("model.mps");
JuMP.set_optimizer(m, Gurobi.Optimizer);
JuMP.set_attribute(m, "Threads", 1);
JuMP.optimize!(m)

GRB_i = problematic_index = 262;
X_ref, LB_ref = Ref(NaN), Ref(NaN);
o = JuMP.unsafe_backend(m);
Gurobi.GRBgetdblattrelement(o,  "X", GRB_i,  X_ref) == 0 || error();
Gurobi.GRBgetdblattrelement(o, "LB", GRB_i, LB_ref) == 0 || error();
println("LB = $(LB_ref[]), X = $(X_ref[])")

m_downstream = JuMP.Model(Gurobi.Optimizer);
JuMP.@variable(m_downstream, x[t=1:2] >= 0);
JuMP.@constraint(m_downstream, x[2] <= x[1] + 0.0); # ramp-up constraint
JuMP.fix(x[1], X_ref[]; force = true); # utilize the solution of an upstream solve as the data input to a downstream solve
JuMP.optimize!(m_downstream)
JuMP.termination_status(m_downstream) # INFEASIBLE

julia> import JuMP, Gurobi

julia> m = JuMP.read_from_file("model.mps");

julia> JuMP.set_optimizer(m, Gurobi.Optimizer);

julia> JuMP.set_attribute(m, "Threads", 1);
Set parameter Threads to value 1

julia> JuMP.optimize!(m)
Set parameter Threads to value 1
Gurobi Optimizer version 13.0.1 build v13.0.1rc0 (linux64gpu - "Ubuntu 24.04.4 LTS")

CPU model: AMD EPYC 7763 64-Core Processor, instruction set [SSE2|AVX|AVX2]
Thread count: 128 physical cores, 256 logical processors, using up to 1 threads

GPU model: NVIDIA RTX A6000, CUDA compute version 8.6, NVIDIA driver compatible with CUDA version 13

Non-default parameters:
Threads  1

Optimize a model with 224 rows, 482 columns and 22253 nonzeros (Min)
Model fingerprint: 0x87db107e
Model has 254 linear objective coefficients
Variable types: 355 continuous, 127 integer (127 binary)
Coefficient statistics:
  Matrix range     [3e-07, 6e+01]
  Objective range  [5e-01, 8e+00]
  Bounds range     [3e-01, 2e+01]
  RHS range        [2e-02, 1e+02]

Presolve removed 70 rows and 109 columns
Presolve time: 0.02s
Presolved: 154 rows, 373 columns, 6283 nonzeros
Variable types: 246 continuous, 127 integer (127 binary)
Found heuristic solution: objective 534.6186536

Root relaxation: objective 1.632712e+02, 130 iterations, 0.00 seconds (0.00 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0  163.27121    0    5  534.61865  163.27121  69.5%     -    0s
H    0     0                     193.7974747  163.27121  15.8%     -    0s
H    0     0                     168.8907811  163.27121  3.33%     -    0s
H    0     0                     164.9582913  163.27121  1.02%     -    0s
     0     0  163.62177    0    3  164.95829  163.62177  0.81%     -    0s
H    0     0                     164.9363677  163.62177  0.80%     -    0s
     0     0  163.62177    0    4  164.93637  163.62177  0.80%     -    0s
*    9     0               2     164.8985760  164.89858  0.00%   3.9    0s

Cutting planes:
  Gomory: 1
  Implied bound: 2
  MIR: 4
  StrongCG: 1
  Flow cover: 2

Explored 10 nodes (284 simplex iterations) in 0.11 seconds (0.05 work units)
Thread count was 1 (of 256 available processors)

Solution count 8: 164.899 164.936 164.958 ... 534.619

Optimal solution found (tolerance 1.00e-04)
Warning: max bound violation (1.0317e-06) exceeds tolerance
Best objective 1.648985759635e+02, best bound 1.648985759635e+02, gap 0.0000%

User-callback calls 147, time in user-callback 0.00 sec

julia> GRB_i = problematic_index = 262;

julia> X_ref, LB_ref = Ref(NaN), Ref(NaN);

julia> o = JuMP.unsafe_backend(m);

julia> Gurobi.GRBgetdblattrelement(o,  "X", GRB_i,  X_ref) == 0 || error();

julia> Gurobi.GRBgetdblattrelement(o, "LB", GRB_i, LB_ref) == 0 || error();

julia> println("LB = $(LB_ref[]), X = $(X_ref[])")
LB = 0.0, X = -1.0316601832750387e-6

julia> m_downstream = JuMP.Model(Gurobi.Optimizer);
Set parameter Username
Set parameter LicenseID to value 2781633
Academic license - for non-commercial use only - expires 2027-02-22

julia> JuMP.@variable(m_downstream, x[t=1:2] >= 0);

julia> JuMP.@constraint(m_downstream, x[2] <= x[1] + 0.0); # ramp-up constraint

julia> JuMP.fix(x[1], X_ref[]; force = true);

julia> JuMP.optimize!(m_downstream)
Gurobi Optimizer version 13.0.1 build v13.0.1rc0 (linux64gpu - "Ubuntu 24.04.4 LTS")
Infeasible model
