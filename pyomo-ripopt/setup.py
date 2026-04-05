"""Setup script that forces platform-specific wheels.

The ripopt binary bundled in pyomo_ripopt/bin/ is platform-specific,
so wheels must be tagged per-platform. The empty Extension trick
tells setuptools this is not a pure-Python package.
"""
from setuptools import setup
from setuptools.dist import Distribution


class BinaryDistribution(Distribution):
    """Force platform-specific wheel (not pure Python)."""
    def has_ext_modules(self):
        return True


setup(distclass=BinaryDistribution)
