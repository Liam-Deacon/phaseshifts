# setup.py
from distutils.core import setup
import py2exe, os
import sqlite3 
setup(console=[os.path.join("phaseshifts", "phsh.py")], options={'destdir':os.path.join('dist', 'bin')})