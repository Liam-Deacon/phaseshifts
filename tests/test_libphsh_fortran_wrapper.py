import sys

import pytest


def test_import_libphsh():
    """
    GIVEN a compiled libphsh.f shared object
    WHEN importing phaseshifts.lib.libphsh compiled wrapper extension
    THEN import should be successful
    """
    ext = ".so" if sys.platform != "win32" else ".pyd"
    try:
        import phaseshifts.lib.libphsh  # noqa
    except ModuleNotFoundError:
        pytest.fail(f"libphsh*{ext} has not been compiled")
    except ImportError as err:
        err_message = f"{ext}: ".join(str(err).split(f"{ext}: ")[1:])
        pytest.fail(f"Unable to import compiled libphsh due to: '{err_message}'")
