# Rosenbrock function solved with JuMP + Ripopt
#
# Problem:
#   min  100*(x2 - x1^2)^2 + (1 - x1)^2
#   (unconstrained)
#
# Expected solution: x* = [1.0, 1.0], obj = 0.0
#
# Usage:
#   1. Build ripopt:  cargo build --release
#   2. Run:
#      RIPOPT_LIBRARY_PATH=target/release julia --project=Ripopt.jl examples/jump_rosenbrock.jl

using JuMP
using Ripopt

model = Model(Ripopt.Optimizer)

@variable(model, x[1:2])
set_start_value(x[1], -1.2)
set_start_value(x[2],  1.0)

@NLobjective(model, Min, 100*(x[2] - x[1]^2)^2 + (1 - x[1])^2)

optimize!(model)

println("=" ^ 60)
println("Rosenbrock with JuMP + Ripopt")
println("=" ^ 60)
println("Termination status: ", termination_status(model))
println("Objective value:    ", objective_value(model))
println("Solution:           ", value.(x))
println("Solve time:         ", solve_time(model), " s")
