# Python imports
from ._setuptools_rust_starter import __all__, __doc__
from ._setuptools_rust_starter import *
from splitp.phylogeny import *
from splitp.enums import *
from splitp import alignment
from splitp import constants
from splitp import constructions

name = "splitp"

# Rust bindings
__all__ = __all__ + ["PythonClass"]


class PythonClass:
    def __init__(self, value: int) -> None:
        self.value = value
