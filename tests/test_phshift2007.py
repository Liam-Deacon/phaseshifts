import os
import tempfile

import pytest

import phaseshifts.phshift2007

try:
    from urllib.error import URLError
except ImportError:  # pragma: no cover - Python 2 fallback
    from urllib2 import URLError

try:
    from shutil import which as find_executable
except ImportError:  # pragma: no cover - Python 2 fallback
    from distutils.spawn import find_executable
import ssl


def test_e2e():
    if os.environ.get("PHASESHIFTS_SKIP_PHSHIFT2007_DOWNLOAD"):
        pytest.skip("PHASESHIFTS_SKIP_PHSHIFT2007_DOWNLOAD is set")
    if not find_executable("gfortran"):
        pytest.skip("gfortran is not available")
    with tempfile.TemporaryDirectory() as directory:
        try:
            phaseshifts.phshift2007.do_action("all", directory)
        except (URLError, ssl.SSLError) as exc:
            pytest.skip("phshift2007 download unavailable: %s" % exc)
