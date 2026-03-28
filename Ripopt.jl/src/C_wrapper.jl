# Low-level C wrapper for libripopt.
#
# This module provides direct Julia bindings to the ripopt C API defined
# in ripopt.h.  The ripopt C API uses 0-based indexing for Jacobian and
# Hessian sparsity patterns.

# -------------------------------------------------------------------------
# Return status codes (must match RipoptReturnStatus in ripopt.h)
# -------------------------------------------------------------------------

const SOLVE_SUCCEEDED               = Cint(0)
const ACCEPTABLE_LEVEL              = Cint(1)
const INFEASIBLE_PROBLEM            = Cint(2)
const MAXITER_EXCEEDED              = Cint(5)
const RESTORATION_FAILED            = Cint(6)
const ERROR_IN_STEP_COMPUTATION     = Cint(7)
const NOT_ENOUGH_DEGREES_OF_FREEDOM = Cint(10)
const INVALID_PROBLEM_DEFINITION    = Cint(11)
const INTERNAL_ERROR                = Cint(-1)

# -------------------------------------------------------------------------
# Callback signatures (matching ripopt.h)
#
#   Eval_F_CB:      (n, x, new_x, obj_value, user_data) -> Cint
#   Eval_Grad_F_CB: (n, x, new_x, grad_f, user_data) -> Cint
#   Eval_G_CB:      (n, x, new_x, m, g, user_data) -> Cint
#   Eval_Jac_G_CB:  (n, x, new_x, m, nele_jac, iRow, jCol, values, user_data) -> Cint
#   Eval_H_CB:      (n, x, new_x, obj_factor, m, lambda, new_lambda,
#                     nele_hess, iRow, jCol, values, user_data) -> Cint
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# RipoptProblem — Julia-side wrapper around the opaque C handle
# -------------------------------------------------------------------------

mutable struct RipoptProblem
    # Opaque C handle (from ripopt_create)
    ptr::Ptr{Cvoid}
    # Problem dimensions
    n::Int
    m::Int
    # Solution storage (populated after solve)
    x::Vector{Float64}
    g::Vector{Float64}
    mult_g::Vector{Float64}
    mult_x_L::Vector{Float64}
    mult_x_U::Vector{Float64}
    obj_val::Float64
    status::Cint
    # Callbacks (stored to prevent GC of Base.CFunction objects)
    eval_f::Union{Base.CFunction,Ptr{Cvoid}}
    eval_grad_f::Union{Base.CFunction,Ptr{Cvoid}}
    eval_g::Union{Base.CFunction,Ptr{Cvoid}}
    eval_jac_g::Union{Base.CFunction,Ptr{Cvoid}}
    eval_h::Union{Base.CFunction,Ptr{Cvoid},Nothing}

    function RipoptProblem(
        ptr::Ptr{Cvoid}, n::Int, m::Int,
        eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h,
    )
        prob = new(
            ptr, n, m,
            zeros(Float64, n),     # x
            zeros(Float64, m),     # g
            zeros(Float64, m),     # mult_g
            zeros(Float64, n),     # mult_x_L
            zeros(Float64, n),     # mult_x_U
            0.0,                   # obj_val
            Cint(-100),            # status (unset)
            eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h,
        )
        finalizer(prob) do p
            if p.ptr != C_NULL
                @ccall libripopt.ripopt_free(p.ptr::Ptr{Cvoid})::Cvoid
                p.ptr = C_NULL
            end
        end
        return prob
    end
end

Base.unsafe_convert(::Type{Ptr{Cvoid}}, p::RipoptProblem) = p.ptr

# -------------------------------------------------------------------------
# Problem creation
# -------------------------------------------------------------------------

# Dummy Hessian callback for L-BFGS mode (must be defined before CreateRipoptProblem
# because @cfunction resolves the name at compile time)
function _dummy_eval_h(
    n::Cint, x::Ptr{Float64}, new_x::Cint, obj_factor::Float64,
    m::Cint, lambda::Ptr{Float64}, new_lambda::Cint,
    nele_hess::Cint, iRow::Ptr{Cint}, jCol::Ptr{Cint},
    values::Ptr{Float64}, user_data::Ptr{Cvoid},
)::Cint
    return Cint(0)
end

"""
    CreateRipoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                        eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h)

Create a new ripopt problem. Callbacks must be `@cfunction` pointers with
C-compatible signatures matching the ripopt C API (0-based indexing).
"""
const _CB = Union{Base.CFunction,Ptr{Cvoid}}

function CreateRipoptProblem(
    n::Int, x_L::Vector{Float64}, x_U::Vector{Float64},
    m::Int, g_L::Vector{Float64}, g_U::Vector{Float64},
    nele_jac::Int, nele_hess::Int,
    eval_f::_CB,
    eval_grad_f::_CB,
    eval_g::_CB,
    eval_jac_g::_CB,
    eval_h::Union{_CB,Nothing},
)
    # If no Hessian callback, we need a dummy that returns 0
    eval_h_ptr = if eval_h === nothing
        @cfunction(
            _dummy_eval_h, Cint,
            (Cint, Ptr{Float64}, Cint, Float64, Cint, Ptr{Float64}, Cint,
             Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cvoid})
        )
    else
        eval_h
    end

    ptr = @ccall libripopt.ripopt_create(
        n::Cint,
        x_L::Ptr{Float64},
        x_U::Ptr{Float64},
        m::Cint,
        g_L::Ptr{Float64},
        g_U::Ptr{Float64},
        nele_jac::Cint,
        nele_hess::Cint,
        eval_f::Ptr{Cvoid},
        eval_grad_f::Ptr{Cvoid},
        eval_g::Ptr{Cvoid},
        eval_jac_g::Ptr{Cvoid},
        eval_h_ptr::Ptr{Cvoid},
    )::Ptr{Cvoid}

    if ptr == C_NULL
        error("ripopt_create returned NULL — allocation failure")
    end

    return RipoptProblem(ptr, n, m, eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h)
end

# -------------------------------------------------------------------------
# Options
# -------------------------------------------------------------------------

function AddRipoptNumOption(prob::RipoptProblem, keyword::String, val::Float64)
    ret = @ccall libripopt.ripopt_add_num_option(
        prob::Ptr{Cvoid}, keyword::Cstring, val::Float64,
    )::Cint
    ret == 0 && error("Unknown ripopt numeric option: \"$keyword\"")
    return nothing
end

function AddRipoptIntOption(prob::RipoptProblem, keyword::String, val::Integer)
    ret = @ccall libripopt.ripopt_add_int_option(
        prob::Ptr{Cvoid}, keyword::Cstring, Cint(val)::Cint,
    )::Cint
    ret == 0 && error("Unknown ripopt integer option: \"$keyword\"")
    return nothing
end

function AddRipoptStrOption(prob::RipoptProblem, keyword::String, val::String)
    ret = @ccall libripopt.ripopt_add_str_option(
        prob::Ptr{Cvoid}, keyword::Cstring, val::Cstring,
    )::Cint
    ret == 0 && error("Unknown ripopt string option: \"$keyword\" = \"$val\"")
    return nothing
end

# -------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------

"""
    RipoptSolve(prob::RipoptProblem)

Solve the NLP.  `prob.x` must contain the initial point on entry.
After the call, `prob.x`, `prob.g`, `prob.mult_g`, `prob.mult_x_L`,
`prob.mult_x_U`, `prob.obj_val`, and `prob.status` are populated.
"""
function RipoptSolve(prob::RipoptProblem)
    obj_ref = Ref{Float64}(0.0)
    prob.status = @ccall libripopt.ripopt_solve(
        prob::Ptr{Cvoid},
        prob.x::Ptr{Float64},
        prob.g::Ptr{Float64},
        obj_ref::Ptr{Float64},
        prob.mult_g::Ptr{Float64},
        prob.mult_x_L::Ptr{Float64},
        prob.mult_x_U::Ptr{Float64},
        C_NULL::Ptr{Cvoid},
    )::Cint
    prob.obj_val = obj_ref[]
    return prob.status
end

# -------------------------------------------------------------------------
# Post-solve statistics
# -------------------------------------------------------------------------

function GetIterCount(prob::RipoptProblem)
    return Int(@ccall libripopt.ripopt_get_iter_count(prob::Ptr{Cvoid})::Cint)
end

function GetSolveTime(prob::RipoptProblem)
    return @ccall libripopt.ripopt_get_solve_time(prob::Ptr{Cvoid})::Float64
end

function GetPrimalInf(prob::RipoptProblem)
    return @ccall libripopt.ripopt_get_primal_inf(prob::Ptr{Cvoid})::Float64
end

function GetDualInf(prob::RipoptProblem)
    return @ccall libripopt.ripopt_get_dual_inf(prob::Ptr{Cvoid})::Float64
end

function GetComplInf(prob::RipoptProblem)
    return @ccall libripopt.ripopt_get_compl_inf(prob::Ptr{Cvoid})::Float64
end

# -------------------------------------------------------------------------
# Log callback
# -------------------------------------------------------------------------

function SetLogCallback(prob::RipoptProblem, callback, user_data::Ptr{Cvoid}=C_NULL)
    @ccall libripopt.ripopt_set_log_callback(
        prob::Ptr{Cvoid}, callback::Ptr{Cvoid}, user_data::Ptr{Cvoid},
    )::Cvoid
end
