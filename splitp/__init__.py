name = "splitp"
from splitp.nx_tree import *
from splitp.parsers import *
from splitp.tree_helper_functions import *
from splitp.squangles import *
from splitp.enums import *
from splitp.tree_reconstruction import *
from ._setuptools_rust_starter import *
from ._setuptools_rust_starter import __all__, __doc__

__all__ = __all__ + ["PythonClass"]

class PythonClass:
    def __init__(self, value: int) -> None:
        self.value = value
