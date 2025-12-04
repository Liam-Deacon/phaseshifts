import glob
import os
import sys

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
LIB_DIR = os.path.join(PROJECT_ROOT, "phaseshifts", "lib")


def _has_local_compiled_lib():
    patterns = ("libphsh*.so", "libphsh*.pyd", "libphsh*.dll")
    return any(glob.glob(os.path.join(LIB_DIR, pattern)) for pattern in patterns)


use_source = os.environ.get("PHASESHIFTS_TEST_USE_SOURCE")

# Prefer the installed package (with compiled extension) unless explicitly
# told to use the source tree or a local compiled lib already exists.
if use_source or _has_local_compiled_lib():
    sys.path.insert(0, PROJECT_ROOT)
else:
    # Ensure the repo root does not shadow the installed wheel
    sys.path[:] = [p for p in sys.path if os.path.abspath(p or os.curdir) != PROJECT_ROOT]
