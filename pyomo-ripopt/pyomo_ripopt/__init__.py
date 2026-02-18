"""Pyomo solver plugin for the ripopt NLP solver.

Usage:
    import pyomo_ripopt  # registers 'ripopt' with SolverFactory
    from pyomo.environ import *
    solver = SolverFactory('ripopt')
"""
from pyomo_ripopt.ripopt_solver import RIPOPT

__all__ = ["RIPOPT"]
