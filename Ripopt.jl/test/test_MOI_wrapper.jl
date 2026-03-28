using Test
using Ripopt
import MathOptInterface as MOI

# ===========================================================================
# HS071 via MOI NLPBlock interface
# ===========================================================================

# Custom evaluator for HS071
struct HS071Evaluator <: MOI.AbstractNLPEvaluator end

function MOI.initialize(::HS071Evaluator, features::Vector{Symbol})
    for f in features
        if !(f in [:Grad, :Jac, :Hess])
            error("Unsupported feature: $f")
        end
    end
end

MOI.features_available(::HS071Evaluator) = [:Grad, :Jac, :Hess]

function MOI.eval_objective(::HS071Evaluator, x)
    return x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3]
end

function MOI.eval_objective_gradient(::HS071Evaluator, grad, x)
    grad[1] = x[4]*(2*x[1]+x[2]+x[3])
    grad[2] = x[1]*x[4]
    grad[3] = x[1]*x[4] + 1.0
    grad[4] = x[1]*(x[1]+x[2]+x[3])
end

function MOI.eval_constraint(::HS071Evaluator, g, x)
    g[1] = x[1]*x[2]*x[3]*x[4]
    g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
end

function MOI.jacobian_structure(::HS071Evaluator)
    # (row, col) 1-based
    return [
        (1,1), (1,2), (1,3), (1,4),
        (2,1), (2,2), (2,3), (2,4),
    ]
end

function MOI.eval_constraint_jacobian(::HS071Evaluator, J, x)
    # Row 1: d(x1*x2*x3*x4)/dx
    J[1] = x[2]*x[3]*x[4]
    J[2] = x[1]*x[3]*x[4]
    J[3] = x[1]*x[2]*x[4]
    J[4] = x[1]*x[2]*x[3]
    # Row 2: d(x1^2+x2^2+x3^2+x4^2)/dx
    J[5] = 2*x[1]
    J[6] = 2*x[2]
    J[7] = 2*x[3]
    J[8] = 2*x[4]
end

function MOI.hessian_lagrangian_structure(::HS071Evaluator)
    # Lower triangle, 1-based
    return [
        (1,1),         # 1
        (2,1),         # 2
        (2,2),         # 3
        (3,1),         # 4
        (3,2),         # 5
        (3,3),         # 6
        (4,1),         # 7
        (4,2),         # 8
        (4,3),         # 9
        (4,4),         # 10
    ]
end

function MOI.eval_hessian_lagrangian(::HS071Evaluator, H, x, σ, μ)
    H .= 0.0
    # Objective
    H[1] += σ * 2*x[4]                   # d2f/dx1dx1
    H[2] += σ * x[4]                     # d2f/dx2dx1
    H[4] += σ * x[4]                     # d2f/dx3dx1
    H[7] += σ * (2*x[1]+x[2]+x[3])      # d2f/dx4dx1
    H[8] += σ * x[1]                     # d2f/dx4dx2
    H[9] += σ * x[1]                     # d2f/dx4dx3
    # Constraint 1
    H[2] += μ[1] * x[3]*x[4]
    H[4] += μ[1] * x[2]*x[4]
    H[5] += μ[1] * x[1]*x[4]
    H[7] += μ[1] * x[2]*x[3]
    H[8] += μ[1] * x[1]*x[3]
    H[9] += μ[1] * x[1]*x[2]
    # Constraint 2
    H[1]  += μ[2] * 2.0
    H[3]  += μ[2] * 2.0
    H[6]  += μ[2] * 2.0
    H[10] += μ[2] * 2.0
end

@testset "MOI Wrapper — HS071 via NLPBlock" begin
    model = Ripopt.Optimizer()
    MOI.set(model, MOI.Silent(), true)

    # Add 4 variables with bounds [1, 5]
    x = MOI.add_variables(model, 4)
    for xi in x
        MOI.add_constraint(model, xi, MOI.GreaterThan(1.0))
        MOI.add_constraint(model, xi, MOI.LessThan(5.0))
    end

    # Set initial point
    start = [1.0, 5.0, 5.0, 1.0]
    for (i, xi) in enumerate(x)
        MOI.set(model, MOI.VariablePrimalStart(), xi, start[i])
    end

    # NLP block
    evaluator = HS071Evaluator()
    nlp_data = MOI.NLPBlockData(
        [MOI.NLPBoundsPair(25.0, Inf), MOI.NLPBoundsPair(40.0, 40.0)],
        evaluator,
        true,  # has_objective
    )
    MOI.set(model, MOI.NLPBlock(), nlp_data)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    MOI.optimize!(model)

    @test MOI.get(model, MOI.TerminationStatus()) == MOI.LOCALLY_SOLVED
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(model, MOI.ResultCount()) == 1

    obj = MOI.get(model, MOI.ObjectiveValue())
    @test isapprox(obj, 17.014, atol=0.01)

    x_sol = [MOI.get(model, MOI.VariablePrimal(), xi) for xi in x]
    @test isapprox(x_sol[1], 1.0, atol=0.05)
    @test isapprox(x_sol[2], 4.743, atol=0.05)
    @test isapprox(x_sol[3], 3.821, atol=0.05)
    @test isapprox(x_sol[4], 1.379, atol=0.05)

    # Duals
    duals = MOI.get(model, MOI.NLPBlockDual())
    @test length(duals) == 2

    # Solve stats
    @test MOI.get(model, MOI.SolveTimeSec()) >= 0.0
    @test MOI.get(model, MOI.BarrierIterations()) > 0

    # Solver name
    @test MOI.get(model, MOI.SolverName()) == "Ripopt"
end

@testset "MOI Wrapper — Rosenbrock (unconstrained)" begin
    # Unconstrained Rosenbrock via NLPBlock

    struct RosenbrockEvaluator <: MOI.AbstractNLPEvaluator end

    MOI.features_available(::RosenbrockEvaluator) = [:Grad, :Jac, :Hess]
    function MOI.initialize(::RosenbrockEvaluator, features::Vector{Symbol})
        for f in features
            @assert f in [:Grad, :Jac, :Hess]
        end
    end

    MOI.eval_objective(::RosenbrockEvaluator, x) =
        100*(x[2]-x[1]^2)^2 + (1-x[1])^2

    function MOI.eval_objective_gradient(::RosenbrockEvaluator, g, x)
        g[1] = -400*x[1]*(x[2]-x[1]^2) - 2*(1-x[1])
        g[2] = 200*(x[2]-x[1]^2)
    end

    MOI.eval_constraint(::RosenbrockEvaluator, g, x) = nothing
    MOI.jacobian_structure(::RosenbrockEvaluator) = Tuple{Int,Int}[]
    MOI.eval_constraint_jacobian(::RosenbrockEvaluator, J, x) = nothing

    MOI.hessian_lagrangian_structure(::RosenbrockEvaluator) =
        [(1,1), (2,1), (2,2)]

    function MOI.eval_hessian_lagrangian(::RosenbrockEvaluator, H, x, σ, μ)
        H[1] = σ * (-400*(x[2]-3*x[1]^2) + 2)
        H[2] = σ * (-400*x[1])
        H[3] = σ * 200.0
    end

    model = Ripopt.Optimizer()
    MOI.set(model, MOI.Silent(), true)

    x = MOI.add_variables(model, 2)
    MOI.set(model, MOI.VariablePrimalStart(), x[1], -1.2)
    MOI.set(model, MOI.VariablePrimalStart(), x[2], 1.0)

    evaluator = RosenbrockEvaluator()
    nlp_data = MOI.NLPBlockData(
        MOI.NLPBoundsPair[],  # no constraints
        evaluator,
        true,
    )
    MOI.set(model, MOI.NLPBlock(), nlp_data)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    MOI.optimize!(model)

    @test MOI.get(model, MOI.TerminationStatus()) == MOI.LOCALLY_SOLVED
    obj = MOI.get(model, MOI.ObjectiveValue())
    @test isapprox(obj, 0.0, atol=1e-6)

    x_sol = [MOI.get(model, MOI.VariablePrimal(), xi) for xi in x]
    @test isapprox(x_sol[1], 1.0, atol=1e-4)
    @test isapprox(x_sol[2], 1.0, atol=1e-4)
end

@testset "MOI Wrapper — MAX_SENSE" begin
    # Maximize -(x-2)^2 = minimize (x-2)^2 but with MAX_SENSE
    # Solution: x = 2, obj = 0

    struct NegQuadEvaluator <: MOI.AbstractNLPEvaluator end

    MOI.features_available(::NegQuadEvaluator) = [:Grad, :Jac, :Hess]
    function MOI.initialize(::NegQuadEvaluator, features::Vector{Symbol})
        for f in features; @assert f in [:Grad, :Jac, :Hess]; end
    end

    MOI.eval_objective(::NegQuadEvaluator, x) = -(x[1]-2)^2
    function MOI.eval_objective_gradient(::NegQuadEvaluator, g, x)
        g[1] = -2*(x[1]-2)
    end
    MOI.eval_constraint(::NegQuadEvaluator, g, x) = nothing
    MOI.jacobian_structure(::NegQuadEvaluator) = Tuple{Int,Int}[]
    MOI.eval_constraint_jacobian(::NegQuadEvaluator, J, x) = nothing
    MOI.hessian_lagrangian_structure(::NegQuadEvaluator) = [(1,1)]
    function MOI.eval_hessian_lagrangian(::NegQuadEvaluator, H, x, σ, μ)
        H[1] = σ * (-2.0)
    end

    model = Ripopt.Optimizer()
    MOI.set(model, MOI.Silent(), true)

    x = MOI.add_variables(model, 1)
    MOI.set(model, MOI.VariablePrimalStart(), x[1], 0.0)

    evaluator = NegQuadEvaluator()
    nlp_data = MOI.NLPBlockData(MOI.NLPBoundsPair[], evaluator, true)
    MOI.set(model, MOI.NLPBlock(), nlp_data)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)

    MOI.optimize!(model)

    @test MOI.get(model, MOI.TerminationStatus()) == MOI.LOCALLY_SOLVED
    obj = MOI.get(model, MOI.ObjectiveValue())
    @test isapprox(obj, 0.0, atol=1e-6)

    x_sol = MOI.get(model, MOI.VariablePrimal(), x[1])
    @test isapprox(x_sol, 2.0, atol=1e-4)
end

@testset "MOI Wrapper — empty model" begin
    model = Ripopt.Optimizer()
    @test MOI.is_empty(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED

    MOI.add_variable(model)
    @test !MOI.is_empty(model)

    MOI.empty!(model)
    @test MOI.is_empty(model)
end

@testset "MOI Wrapper — options" begin
    model = Ripopt.Optimizer()

    # Silent
    @test MOI.get(model, MOI.Silent()) == false
    MOI.set(model, MOI.Silent(), true)
    @test MOI.get(model, MOI.Silent()) == true

    # TimeLimitSec
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing
    MOI.set(model, MOI.TimeLimitSec(), 10.0)
    @test MOI.get(model, MOI.TimeLimitSec()) == 10.0
    MOI.set(model, MOI.TimeLimitSec(), nothing)
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing

    # Raw options
    MOI.set(model, MOI.RawOptimizerAttribute("tol"), 1e-10)
    @test MOI.get(model, MOI.RawOptimizerAttribute("tol")) == 1e-10
end
