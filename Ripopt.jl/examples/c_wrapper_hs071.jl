# HS071 solved directly with the low-level C wrapper (no JuMP needed).
#
# This demonstrates how to use Ripopt.CreateRipoptProblem / RipoptSolve
# directly, which may be useful for performance-critical applications
# or when integrating with non-JuMP frameworks.
#
# Usage:
#   1. Build ripopt:  cargo build --release
#   2. Run:
#      RIPOPT_LIBRARY_PATH=target/release julia --project=Ripopt.jl examples/c_wrapper_hs071.jl

using Ripopt

# ---- Problem dimensions ----
n = 4              # variables
m = 2              # constraints
nele_jac = 8       # nonzeros in Jacobian (dense 2×4)
nele_hess = 10     # nonzeros in lower-triangle Hessian (4×4)

# ---- Callbacks (C-compatible signatures, 0-based indexing) ----

function eval_f(
    n::Cint, x_ptr::Ptr{Float64}, ::Cint,
    obj::Ptr{Float64}, ::Ptr{Cvoid},
)::Cint
    x = unsafe_wrap(Array, x_ptr, Int(n))
    unsafe_store!(obj, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
    return Cint(1)
end

function eval_grad_f(
    n::Cint, x_ptr::Ptr{Float64}, ::Cint,
    g_ptr::Ptr{Float64}, ::Ptr{Cvoid},
)::Cint
    x = unsafe_wrap(Array, x_ptr, Int(n))
    g = unsafe_wrap(Array, g_ptr, Int(n))
    g[1] = x[4]*(2*x[1]+x[2]+x[3])
    g[2] = x[1]*x[4]
    g[3] = x[1]*x[4] + 1.0
    g[4] = x[1]*(x[1]+x[2]+x[3])
    return Cint(1)
end

function eval_g(
    n_::Cint, x_ptr::Ptr{Float64}, ::Cint,
    m_::Cint, g_ptr::Ptr{Float64}, ::Ptr{Cvoid},
)::Cint
    x = unsafe_wrap(Array, x_ptr, Int(n_))
    g = unsafe_wrap(Array, g_ptr, Int(m_))
    g[1] = x[1]*x[2]*x[3]*x[4]
    g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
    return Cint(1)
end

function eval_jac_g(
    n_::Cint, x_ptr::Ptr{Float64}, ::Cint,
    m_::Cint, nele::Cint,
    iRow_ptr::Ptr{Cint}, jCol_ptr::Ptr{Cint},
    val_ptr::Ptr{Float64}, ::Ptr{Cvoid},
)::Cint
    if val_ptr == C_NULL
        iRow = unsafe_wrap(Array, iRow_ptr, Int(nele))
        jCol = unsafe_wrap(Array, jCol_ptr, Int(nele))
        # 0-based indexing!
        iRow[1]=0; jCol[1]=0; iRow[2]=0; jCol[2]=1
        iRow[3]=0; jCol[3]=2; iRow[4]=0; jCol[4]=3
        iRow[5]=1; jCol[5]=0; iRow[6]=1; jCol[6]=1
        iRow[7]=1; jCol[7]=2; iRow[8]=1; jCol[8]=3
    else
        x = unsafe_wrap(Array, x_ptr, Int(n_))
        v = unsafe_wrap(Array, val_ptr, Int(nele))
        v[1]=x[2]*x[3]*x[4]; v[2]=x[1]*x[3]*x[4]
        v[3]=x[1]*x[2]*x[4]; v[4]=x[1]*x[2]*x[3]
        v[5]=2*x[1]; v[6]=2*x[2]; v[7]=2*x[3]; v[8]=2*x[4]
    end
    return Cint(1)
end

function eval_h(
    n_::Cint, x_ptr::Ptr{Float64}, ::Cint,
    obj_factor::Float64,
    m_::Cint, lambda_ptr::Ptr{Float64}, ::Cint,
    nele::Cint,
    iRow_ptr::Ptr{Cint}, jCol_ptr::Ptr{Cint},
    val_ptr::Ptr{Float64}, ::Ptr{Cvoid},
)::Cint
    if val_ptr == C_NULL
        iRow = unsafe_wrap(Array, iRow_ptr, Int(nele))
        jCol = unsafe_wrap(Array, jCol_ptr, Int(nele))
        idx = 1
        for i in 0:3, j in 0:i
            iRow[idx] = Cint(i); jCol[idx] = Cint(j); idx += 1
        end
    else
        x = unsafe_wrap(Array, x_ptr, Int(n_))
        lam = unsafe_wrap(Array, lambda_ptr, Int(m_))
        v = unsafe_wrap(Array, val_ptr, Int(nele))
        v .= 0.0
        σ = obj_factor
        v[1] += σ*2*x[4]; v[2] += σ*x[4]; v[4] += σ*x[4]
        v[7] += σ*(2*x[1]+x[2]+x[3]); v[8] += σ*x[1]; v[9] += σ*x[1]
        v[2] += lam[1]*x[3]*x[4]; v[4] += lam[1]*x[2]*x[4]
        v[5] += lam[1]*x[1]*x[4]; v[7] += lam[1]*x[2]*x[3]
        v[8] += lam[1]*x[1]*x[3]; v[9] += lam[1]*x[1]*x[2]
        v[1] += lam[2]*2.0; v[3] += lam[2]*2.0
        v[6] += lam[2]*2.0; v[10] += lam[2]*2.0
    end
    return Cint(1)
end

# ---- Create @cfunction pointers ----
eval_f_cb = @cfunction($eval_f, Cint,
    (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cvoid}))
eval_grad_f_cb = @cfunction($eval_grad_f, Cint,
    (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cvoid}))
eval_g_cb = @cfunction($eval_g, Cint,
    (Cint, Ptr{Float64}, Cint, Cint, Ptr{Float64}, Ptr{Cvoid}))
eval_jac_g_cb = @cfunction($eval_jac_g, Cint,
    (Cint, Ptr{Float64}, Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cvoid}))
eval_h_cb = @cfunction($eval_h, Cint,
    (Cint, Ptr{Float64}, Cint, Float64, Cint, Ptr{Float64}, Cint,
     Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cvoid}))

# ---- Build problem ----
x_L = [1.0, 1.0, 1.0, 1.0]
x_U = [5.0, 5.0, 5.0, 5.0]
g_L = [25.0, 40.0]
g_U = [1e30, 40.0]

prob = Ripopt.CreateRipoptProblem(
    n, x_L, x_U, m, g_L, g_U,
    nele_jac, nele_hess,
    eval_f_cb, eval_grad_f_cb, eval_g_cb, eval_jac_g_cb, eval_h_cb,
)

# ---- Set options ----
Ripopt.AddRipoptNumOption(prob, "tol", 1e-8)

# ---- Initial point ----
prob.x .= [1.0, 5.0, 5.0, 1.0]

# ---- Solve ----
status = Ripopt.RipoptSolve(prob)

# ---- Report ----
println("=" ^ 60)
println("HS071 with Ripopt C Wrapper")
println("=" ^ 60)
println("Status:           ", status == Ripopt.SOLVE_SUCCEEDED ? "Optimal" : "Status $status")
println("Objective:        ", prob.obj_val)
println("Solution:         ", prob.x)
println("Constraint vals:  ", prob.g)
println("λ (constraints):  ", prob.mult_g)
println("z_L (lower bnd):  ", prob.mult_x_L)
println("z_U (upper bnd):  ", prob.mult_x_U)
println("Iterations:       ", Ripopt.GetIterCount(prob))
println("Solve time:       ", Ripopt.GetSolveTime(prob), " s")
