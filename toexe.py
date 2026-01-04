# setup.py
try:
    from distutils.core import setup
except Exception:
    from setuptools import setup
import os
import py2exe

# py2exe registers distutils commands via import side effects.
if not hasattr(py2exe, "__name__"):  # pragma: no cover
    raise ImportError("py2exe failed to import")

setup(
    console=[os.path.join("phaseshifts", "phsh.py")],
    options={"destdir": os.path.join("dist", "bin")},
)
