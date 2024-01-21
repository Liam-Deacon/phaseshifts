"""The subpackage is for library code, likely compiled from Fortran or C sources."""
from ._fortran_lib import (
    is_module_importable,
    compile_f2py_shared_library,
    FORTRAN_LIBS,
    LIBPHSH_MODULE,
)

from .. import settings

# NOTE: Attempt to compile libphsh library on import if not already available, useful when installing from source dist
if not is_module_importable(LIBPHSH_MODULE) and settings.COMPILE_MISSING:
    compile_f2py_shared_library(**FORTRAN_LIBS[LIBPHSH_MODULE])
