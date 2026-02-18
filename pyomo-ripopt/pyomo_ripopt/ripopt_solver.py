"""Ripopt solver plugin for Pyomo.

Registers 'ripopt' with Pyomo's SolverFactory. Uses the AMPL NL/SOL interface,
identical to how Pyomo integrates with IPOPT.

Usage:
    import pyomo_ripopt
    from pyomo.environ import *
    solver = SolverFactory('ripopt')
    result = solver.solve(model)
"""
import shutil

from pyomo.opt import SolverFactory
from pyomo.solvers.plugins.solvers.ASL import ASL


@SolverFactory.register("ripopt", doc="The ripopt NLP solver")
class RIPOPT(ASL):
    """Pyomo solver interface for ripopt via AMPL Solver Library protocol."""

    def __init__(self, **kwds):
        kwds["type"] = "ripopt"
        super().__init__(**kwds)
        self._metasolver = False
        self.options.solver = "ripopt"

    def _default_executable(self):
        return shutil.which("ripopt")
