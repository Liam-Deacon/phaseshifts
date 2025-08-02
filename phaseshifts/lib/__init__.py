"""The subpackage is for library code, likely compiled from Fortran or C sources."""

import warnings

from ._fortran_lib import (
    is_module_importable,
    compile_f2py_shared_library,
    FORTRAN_LIBS,
    LIBPHSH_MODULE,
)

from .. import settings

_COMPILE_WARNING_TEMPLATE = (
    "Unable to use Fortran libphsh => Failed to compile {}: '{}'"
)

# NOTE: Attempt to compile libphsh library on import if not already available, useful when installing from source dist
if not is_module_importable(LIBPHSH_MODULE) and settings.COMPILE_MISSING:
    try:
        compile_f2py_shared_library(**FORTRAN_LIBS[LIBPHSH_MODULE])
    except Exception as err:  # pylint: disable=broad-except
        warnings.warn(
            msg=_COMPILE_WARNING_TEMPLATE.format(LIBPHSH_MODULE, err),  # type: ignore[call-overload]
            category=UserWarning,
        )

        # TODO: Add pure-python imports of critical libphsh functions here for compatibility
        raise err
