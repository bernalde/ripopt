using Test
using Ripopt

# ===========================================================================
# Test HS071: min x1*x4*(x1+x2+x3) + x3
#   s.t.  x1*x2*x3*x4 >= 25
#         x1^2 + x2^2 + x3^2 + x4^2 = 40
#         1 <= x1,x2,x3,x4 <= 5
#
# Solution: x* ≈ [1.0, 4.743, 3.821, 1.379], obj ≈ 17.014
# ===========================================================================

@testset "C Wrapper — HS071" begin
    n = 4
    m = 2
    nele_jac = 8   # dense 2×4 Jacobian
    nele_hess = 10  # lower-triangle of 4×4

    # --- Callbacks ---

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
            # Sparsity pattern (0-based)
            iRow = unsafe_wrap(Array, iRow_ptr, Int(nele))
            jCol = unsafe_wrap(Array, jCol_ptr, Int(nele))
            # Row 0 (constraint 1): all 4 vars
            iRow[1] = 0; jCol[1] = 0
            iRow[2] = 0; jCol[2] = 1
            iRow[3] = 0; jCol[3] = 2
            iRow[4] = 0; jCol[4] = 3
            # Row 1 (constraint 2): all 4 vars
            iRow[5] = 1; jCol[5] = 0
            iRow[6] = 1; jCol[6] = 1
            iRow[7] = 1; jCol[7] = 2
            iRow[8] = 1; jCol[8] = 3
        else
            x = unsafe_wrap(Array, x_ptr, Int(n_))
            v = unsafe_wrap(Array, val_ptr, Int(nele))
            # Constraint 1: x1*x2*x3*x4
            v[1] = x[2]*x[3]*x[4]
            v[2] = x[1]*x[3]*x[4]
            v[3] = x[1]*x[2]*x[4]
            v[4] = x[1]*x[2]*x[3]
            # Constraint 2: x1^2 + x2^2 + x3^2 + x4^2
            v[5] = 2*x[1]
            v[6] = 2*x[2]
            v[7] = 2*x[3]
            v[8] = 2*x[4]
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
            # Lower-triangle sparsity (0-based)
            iRow = unsafe_wrap(Array, iRow_ptr, Int(nele))
            jCol = unsafe_wrap(Array, jCol_ptr, Int(nele))
            idx = 1
            for i in 0:3
                for j in 0:i
                    iRow[idx] = Cint(i)
                    jCol[idx] = Cint(j)
                    idx += 1
                end
            end
        else
            x = unsafe_wrap(Array, x_ptr, Int(n_))
            lam = unsafe_wrap(Array, lambda_ptr, Int(m_))
            v = unsafe_wrap(Array, val_ptr, Int(nele))
            # Zero out
            v .= 0.0
            # Objective Hessian: f = x1*x4*(x1+x2+x3) + x3
            # d2f/dx1dx1 = 2*x4
            v[1] += obj_factor * 2*x[4]
            # d2f/dx2dx1 = x4
            v[2] += obj_factor * x[4]
            # d2f/dx3dx1 = x4
            v[4] += obj_factor * x[4]
            # d2f/dx4dx1 = 2*x1+x2+x3
            v[7] += obj_factor * (2*x[1]+x[2]+x[3])
            # d2f/dx4dx2 = x1
            v[8] += obj_factor * x[1]
            # d2f/dx4dx3 = x1
            v[9] += obj_factor * x[1]

            # Constraint 1 Hessian: g1 = x1*x2*x3*x4
            v[2] += lam[1] * x[3]*x[4]  # d2g1/dx2dx1
            v[4] += lam[1] * x[2]*x[4]  # d2g1/dx3dx1
            v[5] += lam[1] * x[1]*x[4]  # d2g1/dx3dx2
            v[7] += lam[1] * x[2]*x[3]  # d2g1/dx4dx1
            v[8] += lam[1] * x[1]*x[3]  # d2g1/dx4dx2
            v[9] += lam[1] * x[1]*x[2]  # d2g1/dx4dx3

            # Constraint 2 Hessian: g2 = x1^2+x2^2+x3^2+x4^2
            v[1]  += lam[2] * 2.0  # d2g2/dx1dx1
            v[3]  += lam[2] * 2.0  # d2g2/dx2dx2
            v[6]  += lam[2] * 2.0  # d2g2/dx3dx3
            v[10] += lam[2] * 2.0  # d2g2/dx4dx4
        end
        return Cint(1)
    end

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

    x_L = [1.0, 1.0, 1.0, 1.0]
    x_U = [5.0, 5.0, 5.0, 5.0]
    g_L = [25.0, 40.0]
    g_U = [1e30, 40.0]

    prob = Ripopt.CreateRipoptProblem(
        n, x_L, x_U, m, g_L, g_U,
        nele_jac, nele_hess,
        eval_f_cb, eval_grad_f_cb, eval_g_cb, eval_jac_g_cb, eval_h_cb,
    )

    # Set options
    Ripopt.AddRipoptIntOption(prob, "print_level", 0)

    # Initial point
    prob.x .= [1.0, 5.0, 5.0, 1.0]

    # Solve
    status = Ripopt.RipoptSolve(prob)

    @test status == Ripopt.SOLVE_SUCCEEDED

    @test isapprox(prob.obj_val, 17.014, atol=0.01)
    @test isapprox(prob.x[1], 1.0, atol=0.05)
    @test isapprox(prob.x[2], 4.743, atol=0.05)
    @test isapprox(prob.x[3], 3.821, atol=0.05)
    @test isapprox(prob.x[4], 1.379, atol=0.05)

    # Post-solve stats
    @test Ripopt.GetIterCount(prob) > 0
    @test Ripopt.GetSolveTime(prob) >= 0.0
    @test Ripopt.GetPrimalInf(prob) < 1e-6
end

@testset "C Wrapper — Rosenbrock (unconstrained)" begin
    n = 2
    m = 0
    nele_jac = 0
    nele_hess = 3  # lower triangle of 2×2

    function eval_f_rosen(
        n::Cint, x_ptr::Ptr{Float64}, ::Cint,
        obj::Ptr{Float64}, ::Ptr{Cvoid},
    )::Cint
        x = unsafe_wrap(Array, x_ptr, Int(n))
        unsafe_store!(obj, 100*(x[2]-x[1]^2)^2 + (1-x[1])^2)
        return Cint(1)
    end

    function eval_grad_f_rosen(
        n::Cint, x_ptr::Ptr{Float64}, ::Cint,
        g_ptr::Ptr{Float64}, ::Ptr{Cvoid},
    )::Cint
        x = unsafe_wrap(Array, x_ptr, Int(n))
        g = unsafe_wrap(Array, g_ptr, Int(n))
        g[1] = -400*x[1]*(x[2]-x[1]^2) - 2*(1-x[1])
        g[2] = 200*(x[2]-x[1]^2)
        return Cint(1)
    end

    function eval_g_rosen(
        ::Cint, ::Ptr{Float64}, ::Cint,
        ::Cint, ::Ptr{Float64}, ::Ptr{Cvoid},
    )::Cint
        return Cint(1)
    end

    function eval_jac_g_rosen(
        ::Cint, ::Ptr{Float64}, ::Cint,
        ::Cint, ::Cint,
        ::Ptr{Cint}, ::Ptr{Cint},
        ::Ptr{Float64}, ::Ptr{Cvoid},
    )::Cint
        return Cint(1)
    end

    function eval_h_rosen(
        n_::Cint, x_ptr::Ptr{Float64}, ::Cint,
        obj_factor::Float64,
        ::Cint, ::Ptr{Float64}, ::Cint,
        nele::Cint,
        iRow_ptr::Ptr{Cint}, jCol_ptr::Ptr{Cint},
        val_ptr::Ptr{Float64}, ::Ptr{Cvoid},
    )::Cint
        if val_ptr == C_NULL
            iRow = unsafe_wrap(Array, iRow_ptr, Int(nele))
            jCol = unsafe_wrap(Array, jCol_ptr, Int(nele))
            iRow[1] = 0; jCol[1] = 0  # (1,1)
            iRow[2] = 1; jCol[2] = 0  # (2,1)
            iRow[3] = 1; jCol[3] = 1  # (2,2)
        else
            x = unsafe_wrap(Array, x_ptr, Int(n_))
            v = unsafe_wrap(Array, val_ptr, Int(nele))
            v[1] = obj_factor * (-400*(x[2]-3*x[1]^2) + 2)
            v[2] = obj_factor * (-400*x[1])
            v[3] = obj_factor * 200.0
        end
        return Cint(1)
    end

    eval_f_cb = @cfunction($eval_f_rosen, Cint,
        (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cvoid}))
    eval_grad_f_cb = @cfunction($eval_grad_f_rosen, Cint,
        (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cvoid}))
    eval_g_cb = @cfunction($eval_g_rosen, Cint,
        (Cint, Ptr{Float64}, Cint, Cint, Ptr{Float64}, Ptr{Cvoid}))
    eval_jac_g_cb = @cfunction($eval_jac_g_rosen, Cint,
        (Cint, Ptr{Float64}, Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cvoid}))
    eval_h_cb = @cfunction($eval_h_rosen, Cint,
        (Cint, Ptr{Float64}, Cint, Float64, Cint, Ptr{Float64}, Cint,
         Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cvoid}))

    x_L = [-1e30, -1e30]
    x_U = [1e30, 1e30]
    g_L = Float64[]
    g_U = Float64[]

    prob = Ripopt.CreateRipoptProblem(
        n, x_L, x_U, m, g_L, g_U,
        nele_jac, nele_hess,
        eval_f_cb, eval_grad_f_cb, eval_g_cb, eval_jac_g_cb, eval_h_cb,
    )

    Ripopt.AddRipoptIntOption(prob, "print_level", 0)
    prob.x .= [-1.2, 1.0]

    status = Ripopt.RipoptSolve(prob)

    @test status == Ripopt.SOLVE_SUCCEEDED
    @test isapprox(prob.obj_val, 0.0, atol=1e-6)
    @test isapprox(prob.x[1], 1.0, atol=1e-4)
    @test isapprox(prob.x[2], 1.0, atol=1e-4)
end
