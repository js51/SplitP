# Python imports
from ._setuptools_rust_starter import __all__, __doc__
from ._setuptools_rust_starter import *
from splitp.phylogeny import *
from splitp.enums import *
name = "splitp"

# Rust bindings
__all__ = __all__ + ["PythonClass"]


class PythonClass:
    def __init__(self, value: int) -> None:
        self.value = value
