# Ripopt.jl

Julia interface to [ripopt](https://github.com/jkitchin/ripopt), a nonlinear
interior-point optimizer written in Rust. Provides both a low-level C wrapper
and a [JuMP](https://github.com/jump-dev/JuMP.jl) /
[MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) integration.

## Installation

### 1. Build libripopt

```bash
cd /path/to/ripopt
cargo build --release
```

This produces `target/release/libripopt.so` (Linux) or
`target/release/libripopt.dylib` (macOS).

### 2. Install Ripopt.jl

```julia
using Pkg
Pkg.develop(path="/path/to/ripopt/Ripopt.jl")
```

### 3. Set the library path

Before `using Ripopt`, tell Julia where to find the shared library:

```julia
ENV["RIPOPT_LIBRARY_PATH"] = "/path/to/ripopt/target/release"
using Ripopt
```

Or export the environment variable before starting Julia:

```bash
export RIPOPT_LIBRARY_PATH=/path/to/ripopt/target/release
julia --project=Ripopt.jl
```

## Usage with JuMP

```julia
using JuMP, Ripopt

model = Model(Ripopt.Optimizer)
set_silent(model)

@variable(model, 1 <= x[1:4] <= 5)
set_start_value(x[1], 1.0)
set_start_value(x[2], 5.0)
set_start_value(x[3], 5.0)
set_start_value(x[4], 1.0)

@NLobjective(model, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
@NLconstraint(model, x[1]*x[2]*x[3]*x[4] >= 25)
@NLconstraint(model, x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 == 40)

optimize!(model)

println("Status:   ", termination_status(model))
println("Obj:      ", objective_value(model))
println("Solution: ", value.(x))
```

## Usage with MathOptInterface

See `examples/c_wrapper_hs071.jl` for the low-level C wrapper API and
`test/test_MOI_wrapper.jl` for the MOI interface.

## Solver Options

Options are passed via JuMP's `set_optimizer_attribute` or
`MOI.RawOptimizerAttribute`:

```julia
set_optimizer_attribute(model, "tol", 1e-10)
set_optimizer_attribute(model, "max_iter", 1000)
set_optimizer_attribute(model, "mu_strategy", "adaptive")
set_optimizer_attribute(model, "print_level", 5)
```

Standard MOI attributes are also supported:

```julia
set_silent(model)                          # print_level = 0
set_time_limit_sec(model, 30.0)            # max_wall_time = 30
```

## Running Tests

```bash
cd /path/to/ripopt
cargo build --release
cd Ripopt.jl
RIPOPT_LIBRARY_PATH=../target/release julia --project=. -e 'using Pkg; Pkg.test()'
```

## Architecture

```
Ripopt.jl/
├── Project.toml           # Package metadata and dependencies
├── src/
│   ├── Ripopt.jl          # Module definition, library loading
│   ├── C_wrapper.jl       # Low-level @ccall bindings to libripopt
│   └── MOI_wrapper.jl     # MathOptInterface Optimizer implementation
├── test/
│   ├── runtests.jl        # Test entrypoint
│   ├── test_C_wrapper.jl  # Tests for low-level C API
│   └── test_MOI_wrapper.jl# Tests for MOI/JuMP integration
└── examples/
    ├── jump_hs071.jl      # JuMP example: HS071
    ├── jump_rosenbrock.jl # JuMP example: Rosenbrock
    └── c_wrapper_hs071.jl # Low-level C wrapper example
```

The package has three layers:

1. **C wrapper** (`C_wrapper.jl`): Direct `@ccall` bindings to `libripopt.so`.
   Uses `@cfunction` to create C-compatible callback function pointers.
   Handles 0-based (C) to 1-based (Julia) index conversion.

2. **MOI wrapper** (`MOI_wrapper.jl`): Implements `MOI.AbstractOptimizer`,
   bridging MathOptInterface's `AbstractNLPEvaluator` callbacks to ripopt's
   C callback interface. Supports `NLPBlock`, variable bounds, objective sense,
   warm-starting, and all standard result queries.

3. **JuMP integration**: Automatic via MOI. Use `Model(Ripopt.Optimizer)`.
