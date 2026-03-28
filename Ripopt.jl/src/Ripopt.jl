"""
    Ripopt

Julia interface to the [ripopt](https://github.com/jkitchin/ripopt) nonlinear
interior-point optimizer.  Provides both a low-level C wrapper and a
[MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) /
[JuMP](https://github.com/jump-dev/JuMP.jl) integration.

## Quick start with JuMP

```julia
using JuMP, Ripopt

model = Model(Ripopt.Optimizer)
set_silent(model)
@variable(model, x[1:2])
@NLobjective(model, Min, (x[1] - 1)^2 + (x[2] - 2)^2)
optimize!(model)
value.(x)  # [1.0, 2.0]
```

## Library path

By default, the package looks for `libripopt` in the standard library search
path.  To override, set the environment variable `RIPOPT_LIBRARY_PATH` to the
directory containing `libripopt.so` / `libripopt.dylib` before loading the
package:

```julia
ENV["RIPOPT_LIBRARY_PATH"] = "/path/to/ripopt/target/release"
using Ripopt
```
"""
module Ripopt

# Locate the shared library
const libripopt = if haskey(ENV, "RIPOPT_LIBRARY_PATH")
    joinpath(ENV["RIPOPT_LIBRARY_PATH"], Sys.iswindows() ? "ripopt.dll" :
             Sys.isapple() ? "libripopt.dylib" : "libripopt.so")
else
    "libripopt"
end

include("C_wrapper.jl")
include("MOI_wrapper.jl")

export Optimizer

end # module Ripopt
