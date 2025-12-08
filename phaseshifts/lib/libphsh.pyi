"""
Type stubs for the compiled Fortran extension `libphsh`.
These are manually derived from `libphsh.pyf`.
"""

from typing import Optional, Any

def hartfock(input_file: str) -> None:
    """
    Run the Hartree-Fock (or Dirac-Fock) atomic structure solver.

    Parameters
    ----------
    input_file : str
        Path to the input file (max 255 chars).
    """
    ...

def cavpot(
    mtz_string: str,
    slab_flag: int,
    atomic_file: str,
    cluster_file: str,
    mufftin_file: str,
    output_file: str,
    info_file: str,
) -> float:
    """
    Calculate the cavity potential (muffin-tin zero).
    """
    ...

def hb(x: float, factor: float) -> float:
    """
    Hankel function related utility.
    """
    ...

def phsh_cav(
    mufftin_file: str,
    phasout_file: str,
    dataph_file: str,
    zph_file: str = ...,
) -> Any:
    """
    Calculate phase shifts using the Cavendish method.

    Parameters
    ----------
    mufftin_file : str
        Path to the muffin-tin potential input file.
    phasout_file : str
        Path to the phase shift output file.
    dataph_file : str
        Path to the data phase file.
    zph_file : str, optional
        Path to the zph output file.
    """
    ...

def phsh_wil(
    mufftin_file: str,
    phasout_file: str,
    dataph_file: str,
    *args: Any,
    **kwargs: Any,
) -> Any:
    """
    Calculate phase shifts using Williams' method.

    Parameters
    ----------
    mufftin_file : str
        Path to the muffin-tin potential input file.
    phasout_file : str
        Path to the phase shift output file.
    dataph_file : str
        Path to the data phase file.
    """
    ...

def phsh_rel(
    mufftin_file: str,
    phasout_file: str,
    dataph_file: str,
    *args: Any,
    **kwargs: Any,
) -> Any:
    """
    Calculate relativistic phase shifts.

    Parameters
    ----------
    mufftin_file : str
        Path to the muffin-tin potential input file.
    phasout_file : str
        Path to the phase shift output file.
    dataph_file : str
        Path to the data phase file.
    """
    ...
