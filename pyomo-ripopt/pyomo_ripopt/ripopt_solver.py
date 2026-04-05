"""Ripopt solver plugin for Pyomo.

Registers 'ripopt' with Pyomo's SolverFactory. Uses the AMPL NL/SOL interface,
identical to how Pyomo integrates with IPOPT.

Usage:
    import pyomo_ripopt
    from pyomo.environ import *
    solver = SolverFactory('ripopt')
    result = solver.solve(model)
"""
import os
import platform
import shutil
import sys

from pyomo.opt import SolverFactory
from pyomo.solvers.plugins.solvers.ASL import ASL


def _bundled_binary():
    """Find the ripopt binary bundled inside this wheel (if any)."""
    name = "ripopt.exe" if platform.system() == "Windows" else "ripopt"
    bin_dir = os.path.join(os.path.dirname(__file__), "bin")
    path = os.path.join(bin_dir, name)
    if os.path.isfile(path) and os.access(path, os.X_OK):
        return path
    return None


@SolverFactory.register("ripopt", doc="The ripopt NLP solver")
class RIPOPT(ASL):
    """Pyomo solver interface for ripopt via AMPL Solver Library protocol."""

    def __init__(self, **kwds):
        kwds["type"] = "ripopt"
        super().__init__(**kwds)
        self._metasolver = False
        self.options.solver = "ripopt"

    def _default_executable(self):
        # Prefer the binary bundled in the wheel, fall back to PATH
        return _bundled_binary() or shutil.which("ripopt")
