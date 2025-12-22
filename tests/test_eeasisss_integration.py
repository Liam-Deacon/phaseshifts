import os

import pytest

from phaseshifts.wrappers import EEASiSSSWrapper


def _has_eeasisss_lib():
    try:
        import phaseshifts.lib.libphsh  # noqa: F401
    except Exception:
        return False

    try:
        from phaseshifts.lib.EEASiSSS.EEASiSSS import eeasisss  # noqa: F401
    except Exception:
        return False
    return True


@pytest.mark.skipif(
    os.environ.get("PHASESHIFTS_RUN_EEASISSS_TESTS") != "1",
    reason="Set PHASESHIFTS_RUN_EEASISSS_TESTS=1 to enable EEASiSSS integration tests.",
)
def test_eeasisss_integration(tmp_path):
    if not _has_eeasisss_lib():
        pytest.skip("EEASiSSS or libphsh is not available in this environment.")

    inputx = os.environ.get("PHASESHIFTS_EEASISSS_INPUTX")
    bulk = os.environ.get("PHASESHIFTS_EEASISSS_BULK")
    slab = os.environ.get("PHASESHIFTS_EEASISSS_SLAB")

    if not inputx or not bulk or not slab:
        pytest.skip(
            "Provide PHASESHIFTS_EEASISSS_INPUTX, PHASESHIFTS_EEASISSS_BULK, "
            "and PHASESHIFTS_EEASISSS_SLAB to run this test."
        )

    output = EEASiSSSWrapper.autogen_from_input(
        bulk,
        slab,
        tmp_dir=str(tmp_path),
        inputX=inputx,
        format="cleed",
        store=str(tmp_path),
    )

    assert output
