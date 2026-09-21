"""scEPS (single-cell Expression exPlainability Statistics).

Integrates GWAS and single-cell disease cell atlas data to identify
disease-associated cell neighborhoods.

The Python API lives in :mod:`sceps.sceps_core`::

    from sceps.sceps_core import *

Nothing heavyweight is imported here on purpose: the console scripts in
``sceps.cli`` import this package first, and pulling scanpy in at this point
would slow down ``--help`` for every command.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("sceps")
except PackageNotFoundError:  # running from a source tree without an install
    __version__ = "0.0.0.dev0"

__all__ = ["__version__"]
