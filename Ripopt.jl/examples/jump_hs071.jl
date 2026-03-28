# HS071 solved with JuMP + Ripopt
#
# Problem:
#   min  x1*x4*(x1+x2+x3) + x3
#   s.t. x1*x2*x3*x4 >= 25
#        x1^2 + x2^2 + x3^2 + x4^2 = 40
#        1 <= x1, x2, x3, x4 <= 5
#
# Expected solution: x* ≈ [1.0, 4.743, 3.821, 1.379], obj ≈ 17.014
#
# Usage:
#   1. Build ripopt:  cargo build --release
#   2. Run:
#      RIPOPT_LIBRARY_PATH=target/release julia --project=Ripopt.jl examples/jump_hs071.jl

using JuMP
using Ripopt

model = Model(Ripopt.Optimizer)

# Variables with bounds
@variable(model, 1 <= x[1:4] <= 5)

# Starting point
set_start_value(x[1], 1.0)
set_start_value(x[2], 5.0)
set_start_value(x[3], 5.0)
set_start_value(x[4], 1.0)

# Objective
@NLobjective(model, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])

# Constraints
@NLconstraint(model, x[1]*x[2]*x[3]*x[4] >= 25)
@NLconstraint(model, x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 == 40)

# Solve
optimize!(model)

# Report
println("=" ^ 60)
println("HS071 with JuMP + Ripopt")
println("=" ^ 60)
println("Termination status: ", termination_status(model))
println("Primal status:      ", primal_status(model))
println("Objective value:    ", objective_value(model))
println("Solution:           ", value.(x))
println("Solve time:         ", solve_time(model), " s")
