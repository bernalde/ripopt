# Installation

## Prerequisites: Rust and Cargo

ripopt is written in Rust. Install the Rust toolchain via [rustup](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

Verify: `rustc --version && cargo --version`

## Install ripopt

```bash
git clone https://github.com/jkitchin/ripopt.git
cd ripopt
make install
```

This builds the release binary and shared library, then installs:
- `ripopt` AMPL solver binary → `~/.cargo/bin/`
- `libripopt.dylib` / `libripopt.so` → `~/.local/lib/`

Verify: `ripopt --version`

## Using ripopt as a Rust library

Add to your `Cargo.toml`:

```toml
[dependencies]
ripopt = { git = "https://github.com/jkitchin/ripopt" }
```

## Using ripopt with Python/Pyomo

```bash
pip install ./pyomo-ripopt
```

```python
from pyomo.environ import *
solver = SolverFactory('ripopt')
result = solver.solve(model, tee=True)
```

## Using ripopt with Julia/JuMP

```bash
cargo build --release
julia -e 'import Pkg; Pkg.develop(path="Ripopt.jl")'
```

```julia
ENV["RIPOPT_LIBRARY_PATH"] = "/path/to/ripopt/target/release"
using JuMP, Ripopt

model = Model(Ripopt.Optimizer)
@variable(model, 1 <= x[1:4] <= 5)
set_start_value.(x, [1.0, 5.0, 5.0, 1.0])
@NLobjective(model, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
@NLconstraint(model, x[1]*x[2]*x[3]*x[4] >= 25)
@NLconstraint(model, x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 == 40)
optimize!(model)
println(objective_value(model))  # ≈ 17.014
```

## Using ripopt with GAMS

```bash
cargo build --release
make -C gams && sudo make -C gams install
```

```gams
option nlp = ripopt;
Solve mymodel using nlp minimizing obj;
```

Options via `ripopt.opt` (same format as Ipopt):
```
tol 1e-8
max_iter 1000
print_level 5
```

## Using the C API

After `make install`, link against the shared library using `ripopt.h`:

```bash
cc my_program.c -I/path/to/ripopt -L~/.local/lib -lripopt -lm
```

## Uninstall

```bash
make uninstall
```
