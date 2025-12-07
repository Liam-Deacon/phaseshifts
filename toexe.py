# setup.py
try:
    from distutils.core import setup
except Exception:
    from setuptools import setup
import py2exe, os

setup(
    console=[os.path.join("phaseshifts", "phsh.py")],
    options={"destdir": os.path.join("dist", "bin")},
)
