"""Helpers for locating bundled WebAssembly artifacts."""

from __future__ import absolute_import

import os


def get_dist_dir():
    """Return the on-disk location of bundled WASM assets."""
    return os.path.join(os.path.dirname(__file__), "dist")
