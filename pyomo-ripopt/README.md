# pyomo-ripopt

Pyomo solver plugin for [ripopt](https://github.com/jkitchin/ripopt), a fast interior-point NLP solver written in Rust.

## Installation

```bash
pip install pyomo-ripopt
```

This installs the solver plugin and a bundled `ripopt` binary. No Rust toolchain needed.

## Usage

```python
import pyomo_ripopt  # registers the solver
from pyomo.environ import *

model = ConcreteModel()
model.x = Var(initialize=0.5)
model.obj = Objective(expr=(model.x - 2)**2)

solver = SolverFactory('ripopt')
result = solver.solve(model, tee=True)
print(f"x* = {value(model.x)}")  # 2.0
```

## Solver Options

Pass options the same way as Ipopt:

```python
solver = SolverFactory('ripopt')
solver.options['max_iter'] = 1000
solver.options['tol'] = 1e-10
solver.options['print_level'] = 5
```

## Building from Source

If a pre-built wheel is not available for your platform:

```bash
cargo install ripopt   # installs the ripopt binary
pip install pyomo-ripopt --no-binary :all:
```

The solver will find the `ripopt` binary on your PATH.

## License

EPL-2.0, same as ripopt.
