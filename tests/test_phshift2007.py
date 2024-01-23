import tempfile
import phaseshifts.phshift2007


def test_e2e():
    with tempfile.TemporaryDirectory() as directory:
        phaseshifts.phshift2007.do_action("all", directory)
