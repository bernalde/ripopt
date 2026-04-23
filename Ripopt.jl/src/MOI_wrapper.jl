# MathOptInterface wrapper for Ripopt.
#
# This provides the `Ripopt.Optimizer` type that can be used with JuMP:
#
#   using JuMP, Ripopt
#   model = Model(Ripopt.Optimizer)

import MathOptInterface as MOI

# -------------------------------------------------------------------------
# Optimizer struct
# -------------------------------------------------------------------------

mutable struct Optimizer <: MOI.AbstractOptimizer
    # The low-level C problem handle (nothing until optimize! is called)
    inner::Union{Nothing,RipoptProblem}

    # Problem data collected from MOI
    variables::MOI.Utilities.VariablesContainer{Float64}
    nlp_data::Union{Nothing,MOI.NLPBlockData}

    # Sense: MIN_SENSE, MAX_SENSE, or FEASIBILITY_SENSE
    sense::MOI.OptimizationSense

    # Options to pass to ripopt
    options::Dict{String,Any}
    silent::Bool

    # Warm-start data
    variable_primal_start::Vector{Union{Nothing,Float64}}
    mult_x_L_start::Vector{Union{Nothing,Float64}}
    mult_x_U_start::Vector{Union{Nothing,Float64}}
    nlp_dual_start::Union{Nothing,Vector{Float64}}

    # Solve statistics (set after optimize!)
    solve_time::Float64
    iter_count::Int

    function Optimizer(; kwargs...)
        opt = new(
            nothing,
            MOI.Utilities.VariablesContainer{Float64}(),
            nothing,
            MOI.MIN_SENSE,
            Dict{String,Any}(),
            false,
            Union{Nothing,Float64}[],
            Union{Nothing,Float64}[],
            Union{Nothing,Float64}[],
            nothing,
            0.0,
            0,
        )
        for (key, val) in kwargs
            MOI.set(opt, MOI.RawOptimizerAttribute(String(key)), val)
        end
        return opt
    end
end

# -------------------------------------------------------------------------
# MOI.is_empty / MOI.empty!
# -------------------------------------------------------------------------

function MOI.is_empty(model::Optimizer)
    return MOI.is_empty(model.variables) &&
           model.nlp_data === nothing &&
           model.sense == MOI.MIN_SENSE
end

function MOI.empty!(model::Optimizer)
    model.inner = nothing
    model.variables = MOI.Utilities.VariablesContainer{Float64}()
    model.nlp_data = nothing
    model.sense = MOI.MIN_SENSE
    model.variable_primal_start = Union{Nothing,Float64}[]
    model.mult_x_L_start = Union{Nothing,Float64}[]
    model.mult_x_U_start = Union{Nothing,Float64}[]
    model.nlp_dual_start = nothing
    model.solve_time = 0.0
    model.iter_count = 0
    return
end

# -------------------------------------------------------------------------
# Solver name
# -------------------------------------------------------------------------

MOI.get(::Optimizer, ::MOI.SolverName) = "Ripopt"
MOI.get(::Optimizer, ::MOI.SolverVersion) = "0.2.0"

# -------------------------------------------------------------------------
# copy_to — required for JuMP to transfer the model to our optimizer
# -------------------------------------------------------------------------

# Declare that we support the incremental interface (add_variable, set, etc.)
# so that MOI.Utilities.default_copy_to can build the model incrementally.
MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

# -------------------------------------------------------------------------
# RawOptimizerAttribute (pass-through options)
# -------------------------------------------------------------------------

function MOI.supports(::Optimizer, attr::MOI.RawOptimizerAttribute)
    return true
end

function MOI.get(model::Optimizer, attr::MOI.RawOptimizerAttribute)
    return get(model.options, attr.name) do
        error("RawOptimizerAttribute \"$(attr.name)\" not set")
    end
end

function MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, val)
    model.options[attr.name] = val
    return
end

# -------------------------------------------------------------------------
# Silent
# -------------------------------------------------------------------------

MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

function MOI.set(model::Optimizer, ::MOI.Silent, val::Bool)
    model.silent = val
    return
end

# -------------------------------------------------------------------------
# TimeLimitSec
# -------------------------------------------------------------------------

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return get(model.options, "max_wall_time", nothing)
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, val::Union{Nothing,Real})
    if val === nothing
        delete!(model.options, "max_wall_time")
    else
        model.options["max_wall_time"] = Float64(val)
    end
    return
end

# -------------------------------------------------------------------------
# Variables
# -------------------------------------------------------------------------

function MOI.add_variable(model::Optimizer)
    vi = MOI.add_variable(model.variables)
    push!(model.variable_primal_start, nothing)
    push!(model.mult_x_L_start, nothing)
    push!(model.mult_x_U_start, nothing)
    return vi
end

function MOI.add_variables(model::Optimizer, n::Int)
    return [MOI.add_variable(model) for _ in 1:n]
end

function MOI.is_valid(model::Optimizer, vi::MOI.VariableIndex)
    return MOI.is_valid(model.variables, vi)
end

# -------------------------------------------------------------------------
# Variable bounds (VariableIndex in {LessThan, GreaterThan, EqualTo, Interval})
# -------------------------------------------------------------------------

const _SETS = Union{
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64},
}

function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.VariableIndex}, ::Type{S}
) where {S<:_SETS}
    return true
end

function MOI.add_constraint(
    model::Optimizer, vi::MOI.VariableIndex, set::S,
) where {S<:_SETS}
    return MOI.add_constraint(model.variables, vi, set)
end

function MOI.is_valid(
    model::Optimizer, ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:_SETS}
    return MOI.is_valid(model.variables, ci)
end

# -------------------------------------------------------------------------
# ObjectiveSense
# -------------------------------------------------------------------------

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.get(model::Optimizer, ::MOI.ObjectiveSense) = model.sense

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    model.sense = sense
    return
end

# -------------------------------------------------------------------------
# NLPBlock (the main nonlinear interface)
# -------------------------------------------------------------------------

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true

function MOI.set(model::Optimizer, ::MOI.NLPBlock, data::MOI.NLPBlockData)
    model.nlp_data = data
    return
end

# -------------------------------------------------------------------------
# Warm-start attributes
# -------------------------------------------------------------------------

function MOI.supports(
    ::Optimizer,
    ::MOI.VariablePrimalStart,
    ::Type{MOI.VariableIndex},
)
    return true
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    val::Union{Nothing,Float64},
)
    model.variable_primal_start[vi.value] = val
    return
end

MOI.supports(::Optimizer, ::MOI.NLPBlockDualStart) = true

function MOI.set(model::Optimizer, ::MOI.NLPBlockDualStart, vals)
    model.nlp_dual_start = (vals === nothing) ? nothing : Float64.(vals)
    return
end

# -------------------------------------------------------------------------
# NumberOfVariables
# -------------------------------------------------------------------------

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)
    return MOI.get(model.variables, MOI.NumberOfVariables())
end

# -------------------------------------------------------------------------
# Status code mapping
# -------------------------------------------------------------------------

const _STATUS_MAP = Dict{Cint,MOI.TerminationStatusCode}(
    SOLVE_SUCCEEDED               => MOI.LOCALLY_SOLVED,
    ACCEPTABLE_LEVEL              => MOI.ALMOST_LOCALLY_SOLVED,
    INFEASIBLE_PROBLEM            => MOI.LOCALLY_INFEASIBLE,
    SEARCH_DIRECTION_TOO_SMALL    => MOI.SLOW_PROGRESS,
    DIVERGING_ITERATES            => MOI.NORM_LIMIT,
    USER_REQUESTED_STOP           => MOI.INTERRUPTED,
    MAXITER_EXCEEDED              => MOI.ITERATION_LIMIT,
    RESTORATION_FAILED            => MOI.NUMERICAL_ERROR,
    ERROR_IN_STEP_COMPUTATION     => MOI.NUMERICAL_ERROR,
    MAX_WALLTIME_EXCEEDED         => MOI.TIME_LIMIT,
    NOT_ENOUGH_DEGREES_OF_FREEDOM => MOI.INVALID_MODEL,
    INVALID_PROBLEM_DEFINITION    => MOI.INVALID_MODEL,
    INVALID_NUMBER_DETECTED       => MOI.INVALID_MODEL,
    INTERNAL_ERROR                => MOI.OTHER_ERROR,
)

const _RAW_STATUS_MAP = Dict{Cint,String}(
    SOLVE_SUCCEEDED               => "Solve_Succeeded",
    ACCEPTABLE_LEVEL              => "Acceptable_Level",
    INFEASIBLE_PROBLEM            => "Infeasible_Problem_Detected",
    SEARCH_DIRECTION_TOO_SMALL    => "Search_Direction_Becomes_Too_Small",
    DIVERGING_ITERATES            => "Diverging_Iterates",
    USER_REQUESTED_STOP           => "User_Requested_Stop",
    MAXITER_EXCEEDED              => "Maximum_Iterations_Exceeded",
    RESTORATION_FAILED            => "Restoration_Failed",
    ERROR_IN_STEP_COMPUTATION     => "Error_In_Step_Computation",
    MAX_WALLTIME_EXCEEDED         => "Maximum_WallTime_Exceeded",
    NOT_ENOUGH_DEGREES_OF_FREEDOM => "Not_Enough_Degrees_Of_Freedom",
    INVALID_PROBLEM_DEFINITION    => "Invalid_Problem_Definition",
    INVALID_NUMBER_DETECTED       => "Invalid_Number_Detected",
    INTERNAL_ERROR                => "Internal_Error",
)

function _primal_status(status::Cint)
    if status == SOLVE_SUCCEEDED
        return MOI.FEASIBLE_POINT
    elseif status == ACCEPTABLE_LEVEL
        return MOI.NEARLY_FEASIBLE_POINT
    elseif status == INFEASIBLE_PROBLEM
        return MOI.INFEASIBLE_POINT
    elseif status in (SEARCH_DIRECTION_TOO_SMALL, DIVERGING_ITERATES,
                      USER_REQUESTED_STOP, MAXITER_EXCEEDED, RESTORATION_FAILED,
                      ERROR_IN_STEP_COMPUTATION, MAX_WALLTIME_EXCEEDED)
        return MOI.UNKNOWN_RESULT_STATUS
    else
        return MOI.NO_SOLUTION
    end
end

function _dual_status(status::Cint)
    if status == SOLVE_SUCCEEDED
        return MOI.FEASIBLE_POINT
    elseif status == ACCEPTABLE_LEVEL
        return MOI.NEARLY_FEASIBLE_POINT
    else
        return MOI.NO_SOLUTION
    end
end

# -------------------------------------------------------------------------
# C callback functions (module-level so @cfunction without $ works on ARM)
#
# All problem data is passed through user_data (a Ref to a NamedTuple).
# -------------------------------------------------------------------------

function _moi_eval_f(
    n_::Cint, x_ptr::Ptr{Float64}, ::Cint,
    obj_value::Ptr{Float64}, user_data::Ptr{Cvoid},
)::Cint
    x = unsafe_wrap(Array, x_ptr, Int(n_))
    data = unsafe_pointer_to_objref(user_data)[]
    val = MOI.eval_objective(data.evaluator, x)
    unsafe_store!(obj_value, data.obj_factor * val)
    return Cint(1)
end

function _moi_eval_grad_f(
    n_::Cint, x_ptr::Ptr{Float64}, ::Cint,
    grad_ptr::Ptr{Float64}, user_data::Ptr{Cvoid},
)::Cint
    x = unsafe_wrap(Array, x_ptr, Int(n_))
    grad = unsafe_wrap(Array, grad_ptr, Int(n_))
    data = unsafe_pointer_to_objref(user_data)[]
    MOI.eval_objective_gradient(data.evaluator, grad, x)
    if data.obj_factor != 1.0
        grad .*= data.obj_factor
    end
    return Cint(1)
end

function _moi_eval_g(
    n_::Cint, x_ptr::Ptr{Float64}, ::Cint,
    m_::Cint, g_ptr::Ptr{Float64}, user_data::Ptr{Cvoid},
)::Cint
    x = unsafe_wrap(Array, x_ptr, Int(n_))
    g = unsafe_wrap(Array, g_ptr, Int(m_))
    data = unsafe_pointer_to_objref(user_data)[]
    MOI.eval_constraint(data.evaluator, g, x)
    return Cint(1)
end

function _moi_eval_jac_g(
    n_::Cint, x_ptr::Ptr{Float64}, ::Cint,
    m_::Cint, nele::Cint,
    iRow_ptr::Ptr{Cint}, jCol_ptr::Ptr{Cint},
    values_ptr::Ptr{Float64}, user_data::Ptr{Cvoid},
)::Cint
    data = unsafe_pointer_to_objref(user_data)[]
    if values_ptr == C_NULL
        iRow = unsafe_wrap(Array, iRow_ptr, Int(nele))
        jCol = unsafe_wrap(Array, jCol_ptr, Int(nele))
        for (k, (r, c)) in enumerate(data.jac_struct)
            iRow[k] = Cint(r - 1)
            jCol[k] = Cint(c - 1)
        end
    else
        x = unsafe_wrap(Array, x_ptr, Int(n_))
        values = unsafe_wrap(Array, values_ptr, Int(nele))
        MOI.eval_constraint_jacobian(data.evaluator, values, x)
    end
    return Cint(1)
end

function _moi_eval_h(
    n_::Cint, x_ptr::Ptr{Float64}, ::Cint,
    obj_f::Float64,
    m_::Cint, lambda_ptr::Ptr{Float64}, ::Cint,
    nele::Cint,
    iRow_ptr::Ptr{Cint}, jCol_ptr::Ptr{Cint},
    values_ptr::Ptr{Float64}, user_data::Ptr{Cvoid},
)::Cint
    data = unsafe_pointer_to_objref(user_data)[]
    if values_ptr == C_NULL
        iRow = unsafe_wrap(Array, iRow_ptr, Int(nele))
        jCol = unsafe_wrap(Array, jCol_ptr, Int(nele))
        for (k, (r, c)) in enumerate(data.hess_struct)
            iRow[k] = Cint(r - 1)
            jCol[k] = Cint(c - 1)
        end
    else
        x = unsafe_wrap(Array, x_ptr, Int(n_))
        lambda = unsafe_wrap(Array, lambda_ptr, Int(m_))
        values = unsafe_wrap(Array, values_ptr, Int(nele))
        MOI.eval_hessian_lagrangian(
            data.evaluator, values, x,
            data.obj_factor * obj_f, lambda,
        )
    end
    return Cint(1)
end

# -------------------------------------------------------------------------
# optimize!
# -------------------------------------------------------------------------

function MOI.optimize!(model::Optimizer)
    nlp_data = model.nlp_data
    if nlp_data === nothing
        error("No NLPBlockData set. Use MOI.set(model, MOI.NLPBlock(), data).")
    end

    evaluator = nlp_data.evaluator
    n = MOI.get(model, MOI.NumberOfVariables())
    m = length(nlp_data.constraint_bounds)

    # Determine if Hessian is available
    features = MOI.features_available(evaluator)
    has_hessian = :Hess in features
    init_features = Symbol[:Grad, :Jac]
    if has_hessian
        push!(init_features, :Hess)
    end
    MOI.initialize(evaluator, init_features)

    # Get sparsity structures
    jac_struct = MOI.jacobian_structure(evaluator)
    nele_jac = length(jac_struct)

    hess_struct = if has_hessian
        MOI.hessian_lagrangian_structure(evaluator)
    else
        Tuple{Int,Int}[]
    end
    nele_hess = length(hess_struct)

    # Variable bounds
    x_L = fill(-1e30, n)
    x_U = fill(1e30, n)
    for ci in MOI.get(model.variables, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.GreaterThan{Float64}}())
        vi = MOI.VariableIndex(ci.value)
        set = MOI.get(model.variables, MOI.ConstraintSet(), ci)
        x_L[vi.value] = set.lower
    end
    for ci in MOI.get(model.variables, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.LessThan{Float64}}())
        vi = MOI.VariableIndex(ci.value)
        set = MOI.get(model.variables, MOI.ConstraintSet(), ci)
        x_U[vi.value] = set.upper
    end
    for ci in MOI.get(model.variables, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.EqualTo{Float64}}())
        vi = MOI.VariableIndex(ci.value)
        set = MOI.get(model.variables, MOI.ConstraintSet(), ci)
        x_L[vi.value] = set.value
        x_U[vi.value] = set.value
    end
    for ci in MOI.get(model.variables, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.Interval{Float64}}())
        vi = MOI.VariableIndex(ci.value)
        set = MOI.get(model.variables, MOI.ConstraintSet(), ci)
        x_L[vi.value] = set.lower
        x_U[vi.value] = set.upper
    end

    # Constraint bounds
    g_L = Float64[b.lower for b in nlp_data.constraint_bounds]
    g_U = Float64[b.upper for b in nlp_data.constraint_bounds]

    # Objective sign (ripopt minimizes; for MAX_SENSE we negate)
    obj_factor = (model.sense == MOI.MAX_SENSE) ? -1.0 : 1.0

    # ---- Build C-compatible callbacks ----
    # We store the evaluator, structures, and obj_factor in a RefValue that
    # gets passed as user_data to the C callbacks.
    callback_data = (
        evaluator = evaluator,
        jac_struct = jac_struct,
        hess_struct = hess_struct,
        obj_factor = obj_factor,
        n = n,
        m = m,
    )
    user_ref = Ref(callback_data)

    eval_f_cb = @cfunction(
        _moi_eval_f, Cint,
        (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cvoid})
    )
    eval_grad_f_cb = @cfunction(
        _moi_eval_grad_f, Cint,
        (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cvoid})
    )
    eval_g_cb = @cfunction(
        _moi_eval_g, Cint,
        (Cint, Ptr{Float64}, Cint, Cint, Ptr{Float64}, Ptr{Cvoid})
    )
    eval_jac_g_cb = @cfunction(
        _moi_eval_jac_g, Cint,
        (Cint, Ptr{Float64}, Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cvoid})
    )
    eval_h_cb = if has_hessian
        @cfunction(
            _moi_eval_h, Cint,
            (Cint, Ptr{Float64}, Cint, Float64, Cint, Ptr{Float64}, Cint,
             Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cvoid})
        )
    else
        nothing
    end

    # Create the ripopt problem
    prob = CreateRipoptProblem(
        n, x_L, x_U, m, g_L, g_U,
        nele_jac, nele_hess,
        eval_f_cb, eval_grad_f_cb, eval_g_cb, eval_jac_g_cb, eval_h_cb,
    )

    # Set options
    if model.silent
        AddRipoptIntOption(prob, "print_level", 0)
    end

    for (key, val) in model.options
        if val isa Integer
            AddRipoptIntOption(prob, key, Int(val))
        elseif val isa Real
            AddRipoptNumOption(prob, key, Float64(val))
        elseif val isa AbstractString
            AddRipoptStrOption(prob, key, String(val))
        end
    end

    # If no Hessian, use L-BFGS approximation
    if !has_hessian
        AddRipoptStrOption(prob, "hessian_approximation", "limited-memory")
    end

    # Set initial point
    for i in 1:n
        start = model.variable_primal_start[i]
        if start !== nothing
            prob.x[i] = start
        else
            # Clamp 0 to bounds
            prob.x[i] = min(max(0.0, x_L[i]), x_U[i])
        end
    end

    # GC protection: prevent user_ref from being collected during solve
    GC.@preserve user_ref callback_data begin
        # Store user_ref pointer as user_data for callbacks.
        # We need to pass this through ripopt_solve's user_data parameter.
        # However, ripopt_solve passes user_data to all callbacks —
        # we use pointer_from_objref on user_ref directly.
        #
        # Override the solve to pass our user_data
        obj_ref = Ref{Float64}(0.0)
        prob.status = @ccall libripopt.ripopt_solve(
            prob.ptr::Ptr{Cvoid},
            prob.x::Ptr{Float64},
            prob.g::Ptr{Float64},
            obj_ref::Ptr{Float64},
            prob.mult_g::Ptr{Float64},
            prob.mult_x_L::Ptr{Float64},
            prob.mult_x_U::Ptr{Float64},
            pointer_from_objref(user_ref)::Ptr{Cvoid},
        )::Cint
        prob.obj_val = obj_ref[]
    end

    model.inner = prob
    model.solve_time = GetSolveTime(prob)
    model.iter_count = GetIterCount(prob)

    return
end

# -------------------------------------------------------------------------
# Results query
# -------------------------------------------------------------------------

function _check_solved(model::Optimizer)
    if model.inner === nothing
        error("Model has not been optimized. Call MOI.optimize! first.")
    end
end

# TerminationStatus
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.inner === nothing
        return MOI.OPTIMIZE_NOT_CALLED
    end
    return get(_STATUS_MAP, model.inner.status, MOI.OTHER_ERROR)
end

# RawStatusString
function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    _check_solved(model)
    return get(_RAW_STATUS_MAP, model.inner.status, "Unknown")
end

# ResultCount
function MOI.get(model::Optimizer, ::MOI.ResultCount)
    if model.inner === nothing
        return 0
    end
    status = model.inner.status
    if status in (SOLVE_SUCCEEDED, ACCEPTABLE_LEVEL, INFEASIBLE_PROBLEM,
                  SEARCH_DIRECTION_TOO_SMALL, DIVERGING_ITERATES, USER_REQUESTED_STOP,
                  MAXITER_EXCEEDED, RESTORATION_FAILED, ERROR_IN_STEP_COMPUTATION,
                  MAX_WALLTIME_EXCEEDED)
        return 1
    end
    return 0
end

# PrimalStatus
function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index > MOI.get(model, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    return _primal_status(model.inner.status)
end

# DualStatus
function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if attr.result_index > MOI.get(model, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    return _dual_status(model.inner.status)
end

# ObjectiveValue
function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    _check_solved(model)
    # Undo the sign flip for MAX_SENSE
    sign = (model.sense == MOI.MAX_SENSE) ? -1.0 : 1.0
    return sign * model.inner.obj_val
end

# VariablePrimal
function MOI.get(
    model::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(model, attr)
    _check_solved(model)
    return model.inner.x[vi.value]
end

# NLPBlockDual
function MOI.get(model::Optimizer, attr::MOI.NLPBlockDual)
    MOI.check_result_index_bounds(model, attr)
    _check_solved(model)
    # ripopt returns multipliers with the same sign convention as Ipopt
    sign = (model.sense == MOI.MAX_SENSE) ? -1.0 : 1.0
    return sign .* model.inner.mult_g
end

# SolveTimeSec
MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = model.solve_time

# BarrierIterations
MOI.get(model::Optimizer, ::MOI.BarrierIterations) = Int64(model.iter_count)
