#!/usr/bin/env python
# encoding: utf-8
"""
Created on 18 Mar 2015

@author: Liam Deacon

@contact: liam.deacon@diamond.ac.uk

@copyright: 2015 Liam Deacon

@license: MIT License

Provides a simple wrapper for executing EEASiSSS

"""

from ctypes import cdll, create_string_buffer
from ctypes.util import find_library

import os
from sys import platform

_ext = ".dll" if str(platform).startswith("win") else ".so"
_lib = os.path.join(os.path.dirname(__file__), "lib")

os.environ["PATH"] = _lib + ";" + os.environ["PATH"]
_library_path = find_library("EEASiSSS") or os.path.join(_lib, "libEEASiSSS" + _ext)
libEEASiSSS = cdll.LoadLibrary(_library_path)


def eeasisss(input_file="inputX"):
    """Wrapper function to call EEASiSSS Fortran library using ctypes"""
    libEEASiSSS.hartfock_(create_string_buffer(str(input_file)), size=255)


if __name__ == "__main__":
    # use this file as a script
    eeasisss()
