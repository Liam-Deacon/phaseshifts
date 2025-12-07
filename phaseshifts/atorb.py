#!/usr/bin/env python
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2013-2014 Liam Deacon                             #
#                                                                            #
# License: MIT License                                                       #
#                                                                            #
# Permission is hereby granted, free of charge, to any person obtaining a    #
# copy of this software and associated documentation files (the "Software"), #
# to deal in the Software without restriction, including without limitation  #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,   #
# and/or sell copies of the Software, and to permit persons to whom the      #
# Software is furnished to do so, subject to the following conditions:       #
#                                                                            #
# The above copyright notice and this permission notice shall be included in #
# all copies or substantial portions of the Software.                        #
#                                                                            #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    #
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        #
# DEALINGS IN THE SOFTWARE.                                                  #
#                                                                            #
##############################################################################
"""
**atorb.py**

This module provides a high-level Python interface for generating input files and
executing atomic structure calculations using the 'atorb' program (part of the
Barbieri/Van Hove LEED phase shift package).

The core functionality allows users to:
1.  Determine the ground-state electronic configuration of elements.
2.  Construct formatted input files for the `atorb` Fortran solver.
3.  Calculate radial charge densities :math:`\\rho(r)` by solving the
    Dirac-Fock (relativistic) or Hartree-Fock (non-relativistic) equations
    for a free atom.

These charge densities serve as the starting potential for calculating
scattering phase shifts in Low-Energy Electron Diffraction (LEED) and
photoelectron diffraction experiments.

:See: http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/

:Requires: f2py (for libphsh fortran wrapper generation)

.. note::
   To generate libphsh fortran wrappers (libphsh.pyd) for your platform
   then use 'python setup.py' in the lib directory of this package to
   install into your python distribution. Alternatively, use::

     f2py -c -m libphsh libphsh.f

   Windows users may have to add appropriate compiler switches, e.g. ::

    # 32-bit
    f2py -c -m libphsh --fcompiler=gfortran --compiler=mingw-32 libphsh.f

    # 64-bit
    f2py -c -m libphsh --fcompiler=gfortran --compiler=mingw-64 libphsh.f

"""

import os
from collections import OrderedDict

# TODO: Clean this up when the project officially drops support for Python 2.7
try:
    from typing import Optional  # noqa: F401
    from typing_extensions import Literal

    ElementBackendType = Literal["mendeleev", "elementy", "periodictable"]
except ModuleNotFoundError:
    # we are likely on an old version of python
    ElementBackendType = str  # type: ignore

try:
    import periodictable  # type: ignore [import-untyped]
except ImportError:
    periodictable = None

try:
    import elementy  # type: ignore [import-not-found]
except ImportError:
    elementy = None

try:
    import mendeleev  # type: ignore [import-not-found]
except ImportError:
    mendeleev = None  # Need to install mendeleev

from phaseshifts.validation.atorb import (
    AtorbElectron,
    AtorbInputModel,
    coerce_model,
    render_atorb_file,
    validate_atorb_file,
)

from phaseshifts import elements

elements_dict = OrderedDict(
    [
        ("H", "Hydrogen"),
        ("He", "Helium"),
        ("Li", "Lithium"),
        ("Be", "Beryllium"),
        ("B", "Boron"),
        ("C", "Carbon"),
        ("N", "Nitrogen"),
        ("O", "Oxygen"),
        ("F", "Fluorine"),
        ("Ne", "Neon"),
        ("Na", "Sodium"),
        ("Mg", "Magnesium"),
        ("Al", "Aluminium"),
        ("Si", "Silicon"),
        ("P", "Phosphorus"),
        ("S", "Sulfur"),
        ("Cl", "Chlorine"),
        ("Ar", "Argon"),
        ("K", "Potassium"),
        ("Ca", "Calcium"),
        ("Sc", "Scandium"),
        ("Ti", "Titanium"),
        ("V", "Vanadium"),
        ("Cr", "Chromium"),
        ("Mn", "Manganese"),
        ("Fe", "Iron"),
        ("Co", "Cobalt"),
        ("Ni", "Nickel"),
        ("Cu", "Copper"),
        ("Zn", "Zinc"),
        ("Ga", "Gallium"),
        ("Ge", "Germanium"),
        ("As", "Arsenic"),
        ("Se", "Selenium"),
        ("Br", "Bromine"),
        ("Kr", "Krypton"),
        ("Rb", "Rubidium"),
        ("Sr", "Strontium"),
        ("Y", "Yttrium"),
        ("Zr", "Zirconium"),
        ("Nb", "Niobium"),
        ("Mo", "Molybdenum"),
        ("Tc", "Technetium"),
        ("Ru", "Ruthenium"),
        ("Rh", "Rhodium"),
        ("Pd", "Palladium"),
        ("Ag", "Silver"),
        ("Cd", "Cadmium"),
        ("In", "Indium"),
        ("Sn", "Tin"),
        ("Sb", "Antimony"),
        ("Te", "Tellurium"),
        ("I", "Iodine"),
        ("Xe", "Xenon"),
        ("Cs", "Cesium"),
        ("Ba", "Barium"),
        ("La", "Lanthanum"),
        ("Ce", "Cerium"),
        ("Pr", "Praseodymium"),
        ("Nd", "Neodymium"),
        ("Pm", "Promethium"),
        ("Sm", "Samarium"),
        ("Eu", "Europium"),
        ("Gd", "Gadolinium"),
        ("Tb", "Terbium"),
        ("Dy", "Dysprosium"),
        ("Ho", "Holmium"),
        ("Er", "Erbium"),
        ("Tm", "Thulium"),
        ("Yb", "Ytterbium"),
        ("Lu", "Lutetium"),
        ("Hf", "Hafnium"),
        ("Ta", "Tantalum"),
        ("W", "Tungsten"),
        ("Re", "Rhenium"),
        ("Os", "Osmium"),
        ("Ir", "Iridium"),
        ("Pt", "Platinum"),
        ("Au", "Gold"),
        ("Hg", "Mercury"),
        ("Tl", "Thallium"),
        ("Pb", "Lead"),
        ("Bi", "Bismuth"),
        ("Po", "Polonium"),
        ("At", "Astatine"),
        ("Rn", "Radon"),
        ("Fr", "Francium"),
        ("Ra", "Radium"),
        ("Ac", "Actinium"),
        ("Th", "Thorium"),
        ("Pa", "Protactinium"),
        ("U", "Uranium"),
        ("Np", "Neptunium"),
        ("Pu", "Plutonium"),
        ("Am", "Americium"),
        ("Cm", "Curium"),
        ("Bk", "Berkelium"),
        ("Cf", "Californium"),
        ("Es", "Einsteinium"),
        ("Fm", "Fermium"),
        ("Md", "Mendelevium"),
        ("No", "Nobelium"),
        ("Lr", "Lawrencium"),
        ("Rf", "Rutherfordium"),
        ("Db", "Dubnium"),
        ("Sg", "Seaborgium"),
        ("Bh", "Bohrium"),
        ("Hs", "Hassium"),
        ("Mt", "Meitnerium"),
        ("Ds", "Darmstadtium"),
        ("Rg", "Roentgenium"),
        ("Cn", "Copernicium"),
        ("Uut", "Ununtrium"),
        ("Fl", "Flerovium"),
        ("Uup", "Ununpentium"),
        ("Lv", "Livermorium"),
        ("Uus", "Ununseptium"),
        ("Uuo", "Ununoctium"),
    ]
)


def get_element(
    element, backend=None
):  # type: (str, Optional[ElementBackendType]) -> object
    """
    Retrieve an element object from the specified chemical data backend.

    This function abstracts the details of various chemistry libraries
    (mendeleev, elementy, periodictable) to provide a unified interface
    for accessing atomic properties like proton count (Z) and electronic
    configurations.

    Parameters
    ----------
    element : str or int
        The symbol (e.g., 'Fe') or atomic number (e.g., 26) of the element.
    backend : str, optional
        The preferred backend library to use ('mendeleev', 'elementy', 'periodictable').
        If None, tries available backends in order.

    Returns
    -------
    object
        An object containing element data (attributes vary by backend but usually
        include `protons`, `symbol`, `name`).

    Raises
    ------
    LookupError
        If the element cannot be found in any available backend.
    """
    ele_obj = elements.ELEMENTS.get(element)
    if mendeleev and not ele_obj and backend in ("mandeleev", None):
        elements_data = mendeleev.get_all_elements()
        elements_data = {e.protons: e for e in elements_data}
        elements_data.update({e.symbol: e for e in elements_data})
        elements_data.update({e.name: e for e in elements_data})
        ele_obj = elements_data.get(element)
    if elementy and not ele_obj and backend in ("elementy", None):
        periodic_table = elementy.PeriodicTable()
        elements_data = {e.protons: e for e in periodic_table.elements}
        elements_data.update({e.symbol: e for e in periodic_table.elements})
        elements_data.update({e.name: e for e in periodic_table.elements})
        ele_obj = elements_data.get(element)
    if periodictable and not ele_obj and backend in ("periodtable", None):
        ele_obj = getattr(periodictable, element, None) or getattr(
            periodictable, str(element).lower()
        )
    if ele_obj is None:
        raise LookupError("Unable to match element {}".format(element))
    return ele_obj


def get_electron_config(element_obj):
    """
    Extract the electronic configuration string from an element object.

    Retrieves the standard ground-state configuration (e.g., '[Ar] 4s2 3d10 4p5')
    from the backend-specific element object.

    Parameters
    ----------
    element_obj : object
        The element object returned by `get_element`.

    Returns
    -------
    str
        The electronic configuration string.
    """
    electron_config = None
    if hasattr(element_obj, "orbitals"):
        electron_config = " ".join(
            ["{orbital}{electrons}".format(**x) for x in element_obj.orbitals]
        )
    else:
        electron_config = getattr(element_obj, "eleconfig", None) or getattr(
            element_obj, "econf", None
        )
    return electron_config


class Atorb(object):
    r"""
    Wrapper for the `atorb` atomic structure solver.

    This class encapsulates the logic for interacting with the Barbieri/Van Hove
    `atorb` program. It calculates atomic orbitals and radial charge densities
    on a logarithmic grid.

    The solver numerically integrates the radial Dirac equation (or Schrödinger
    equation if relativistic effects are disabled) to find the eigenvalues
    and eigenfunctions of the bound electrons.

    Mathematical Context
    --------------------
    The calculations are performed on a logarithmic radial grid defined by:

    .. math::
        r_i = r_{min} \cdot \left(\frac{r_{max}}{r_{min}}\right)^{\frac{i}{N_R}}, \quad i=1, \dots, N_R

    where :math:`N_R` is the number of grid points. Distances are in Bohr radii
    (:math:`a_0 \simeq 0.529 \mathrm{\AA}`).

    The total spherical charge density :math:`\rho(r)` (in electrons per cubic Bohr)
    is constructed from the radial wavefunctions (stored in internal `phe` arrays):

    .. math::
        \rho(r_i) = \sum_{j=1}^{N_{orb}} \frac{occ_j \cdot |\phi_j(r_i)|^2}{4\pi r_i^2}

    where :math:`occ_j` is the occupancy of the :math:`j`-th orbital.

    For relativistic calculations (Dirac-Fock), the orbital density includes
    both large (:math:`F`) and small (:math:`G`) components:

    .. math::
        |\phi(r)|^2 = F(r)^2 + G(r)^2

    Notes
    -----
    The Breit interaction is neglected, which is a valid approximation for
    valence states dominating the scattering potential in LEED.

    Original author: Eric Shirley

    There are nr grid points, and distances are in Bohr radii
    :math:`a_0 \simeq 0.539 \mathrm{\AA}`

    :math:`r(i) = r_{min} \cdot (r_{max} / r_{min})^{(i/n_r)}`,
    :math:`i=1,2,3,...n_r-1,n_r`

    The orbitals are stored in phe(), first index goes :math:`1...n_r`, the
    second index is the orbital index (:math:`i...n_{el}`)

    Look at the atomic files after printing this out to see everything...
    Suffice it to say, that the charge density at radius :math:`r(i)`
    in units of electrons per cubic Bohr radius is given by:

    :math:`\sum_{j-1}^{n_el}{occ(j) \cdot phe(i,j)^2 / (4.0\,\pi\,{r(i)^2)}}`

    Think of the phe functions as plotting the radial wave-functions
    as a function of radius on a logarithmic mesh...

    The Dirac equation is solved for the orbitals, whereas their density
    is treated by setting :math:`phe(i,j)` to Dirac's
    :math:`\sqrt{F(i,j)^2 + G(i,j)^2}` times the sign of :math:`G(i,j)`...

    So we are doing Dirac-Fock, except that we are not treating exchange
    exactly, in terms of working with major and minor components of the
    orbitals, and the phe's give the CORRECT CHARGE DENSITY...

    The above approximation ought to be very small for valence states,
    so you need not worry about it...

    The Breit interaction has been neglected altogether...it should not
    have a huge effect on the charge density you are concerned with...
    """

    def __init__(self, **kwargs):
        """
        Initialize the Atorb wrapper.

        Parameters
        ----------
        **kwargs : dict
            Arbitrary attributes to store on the instance.
        """
        self.__dict__.update(kwargs)

    @staticmethod
    def get_quantum_info(shell):  # (str) -> Tuple[int|float|List[int|float], ...]
        r"""
        Parse quantum numbers from a subshell string (e.g., '3d6').

        Decomposes a standard spectroscopic notation into the set of quantum
        numbers required for the radial equation solver. Handles the splitting
        of shells into spin-orbit coupled states (:math:`j = l \pm 1/2`).

        Parameters
        ----------
        shell : str
            Subshell string, e.g., '1s2', '4f14', '3d5'.

        Returns
        -------
        tuple
            (n, l, j_list, occ_list)

            - **n** (int): Principal quantum number.
            - **l** (int): Azimuthal quantum number.
            - **j_list** (list[float]): Total angular momentum values :math:`j`
              associated with this shell. For :math:`l > 0`, this will be
              [:math:`l-1/2`, :math:`l+1/2`].
            - **occ_list** (list[float]): Occupancy for each :math:`j` level.
              Electrons are distributed according to statistical weight
              (:math:`2j+1`).

        Examples
        --------
        For a full '3d10' shell:

        >>> Atorb.get_quantum_info('3d10')
        (3, 2, [1.5, 2.5], [4.0, 6.0])

        Here, :math:`n=3`, :math:`l=2`. The states are :math:`d_{3/2}` (occupancy 4)
        and :math:`d_{5/2}` (occupancy 6).
        """

        subshell = "".join([s for s in shell if s.isalpha()])
        try:
            (n, nelectrons) = [
                t(s) for t, s in zip((int, int), shell.replace(subshell, " ").split())
            ]
        except ValueError:  # assume 1 electron in shell
            n = int(shell.replace(subshell, " ").split()[0])
            nelectrons = 1

        s = 0.5
        shell_info = None
        if subshell == "s":
            l = 0
            occ = [nelectrons / 1.0]
            j = [l + s]
            shell_info = (n, l, j, occ)
        elif subshell == "p":
            # 3 subshells
            l = 1
            max_occ = 6
            occ = []
            for j in [l - s, l + s]:
                occ.append(((2.0 * j) + 1) * nelectrons / max_occ)
            shell_info = (n, l, [l - s, l + s], occ)
        elif subshell == "d":
            # 5 subshells
            l = 2
            max_occ = 10
            occ = []
            for j in [l - s, l + s]:
                occ.append(((2.0 * j) + 1) * nelectrons / max_occ)
            shell_info = (n, l, [l - s, l + s], occ)
        elif subshell == "f":
            # 7 subshells!
            l = 3
            max_occ = 14
            occ = []
            for j in [l - s, l + s]:
                occ.append(((2.0 * j) + 1) * nelectrons / max_occ)
            shell_info = (n, l, [l - s, l + s], occ)
        else:
            raise NotImplementedError(
                "Exotic sub-shells beyond f-block have not been implemented"
            )
        return shell_info

    @staticmethod
    def replace_core_config(electron_config):
        """
        Expand noble gas core abbreviations into full orbital strings.

        The `atorb` solver requires explicit definitions for all orbitals
        starting from 1s. This function replaces '[Ar]', '[Xe]', etc., with
        their constituent subshells.

        Parameters
        ----------
        electron_config : str
            Electronic configuration string, potentially containing noble gas
            cores (e.g., '[Ar] 4s2').

        Returns
        -------
        str
            The fully expanded configuration string.

        Examples
        --------
        >>> Atorb.replace_core_config('[He] 2s1')
        '1s2 2s1'
        """
        cores = {
            "[He]": "1s2",
            "[Ne]": "1s2 2s2 2p6",
            "[Ar]": "1s2 2s2 2p6 3s2 3p6",
            "[Kr]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6",
            "[Xe]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2 4d10 5p6",
            "[Rn]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2 4d10 5p6" "4f14 5d10 6s2 6p6",
        }
        core = electron_config.split()[0]

        if core in cores:
            config = electron_config.replace(core, cores.get(core))
        else:
            config = electron_config
        return config

    @staticmethod
    def gen_input(element, **kwargs):
        """
        Generate the input file for the 'atorb' atomic structure solver.

        This method constructs a formatted input file required by the Barbieri/Van Hove
        'atorb' program (encapsulated in `libphsh`). It determines the electronic
        configuration of the specified element, handles noble gas core expansion,
        and sets up the radial grid and exchange-correlation parameters for the
        Dirac-Fock calculation.

        Parameters
        ----------
        element : str or int
            The chemical symbol (e.g., 'Cu') or atomic number (e.g., 29) of the target element.
        **kwargs : dict, optional
            Configuration options for the calculation:

            output : str
                Filename for the resulting charge density output (default: 'at_<symbol>.i').
            ngrid : int
                Number of points in the logarithmic radial grid (default: 1000).
                Higher values provide better numerical resolution near the nucleus.
            rel : int or bool
                Relativistic flag (default: 1).
                1 (or True): Solve Dirac equations (includes spin-orbit coupling).
                0 (or False): Solve non-relativistic Schrödinger equations.
            filename : str
                Filename for the generated input text file (default: 'atorb_<symbol>.txt').
            header : str
                Comment line at the top of the input file (default: auto-generated).
            method : str
                Exchange-correlation potential approximation (default: '0.d0' for Hartree-Fock).
                '0.d0': Hartree-Fock (HF).
                '1.d0': Local Density Approximation (LDA).
                '-alpha': X-alpha method (requires providing alpha value).
            relic : float
                Mixing parameter for the self-consistency cycle (default: 0.0).
            mixing_SCF : float
                Mixing parameter for Self-Consistent Field convergence (default: 0.5).
            tolerance : float
                Convergence tolerance for orbital eigenvalues (default: 0.0005 Hartrees).
            ech : int
                Parameter for the exchange potential (default: 100).

        Example
        -------
        >>> Atorb.gen_input('H',rel=False,filename="atorb.txt",output='at_H.i')
        >>> with open('atorb_H.txt', 'r') as f: print("".join(f.readlines())
         C*********************************************************************
         C  atorb input file: atorb_H.txt.
         C*********************************************************************
         i
         1 1000                   ! Z NR (number of points in radial grid)
         d
         0                        ! 1=rel, 0=n.r.
         x
         0.d0                     ! 0.d0=HF, 1.d0=LDA, -alfa = xalfa...
         a
         0 1 0.5 0.0005 100       ! relic,levels,mixing SCF, eigen. tol,for ech
         1 0 0 -0.5 1 1.0         ! n, l, l, -j, <1>, occupation
         w
         at_H.i
         q

        Returns
        -------
        str
            The path to the generated input file.

        Notes
        -----
        The resulting file allows `libphsh` to solve the radial Dirac equation:

        .. math::
            H \\Psi = E \\Psi

        where density functional theory (DFT) or Hartree-Fock approximations define
        the potential.
        """
        ele = get_element(element, backend=None)
        Z = ele.protons

        # get full electronic configuration
        electron_config = get_electron_config(ele)
        config = Atorb.replace_core_config(electron_config)

        # get quantum numbers & occupancy for each electronic orbital in atom
        electrons = []
        nlevels = 0
        for shell in config.split():
            (n, l, J, occ) = Atorb.get_quantum_info(shell)
            for i, j in enumerate(J):
                electrons.append((n, l, l, -j, 1, occ[i]))
                nlevels += 1

        # test kwargs and generate output arguments
        if "output" in kwargs:
            output = kwargs.get("output")
        else:
            output = "at_{0}.i".format(ele.symbol)

        if "ngrid" in kwargs:
            NR = kwargs.get("ngrid")
        else:
            NR = 1000  # default grid resolution

        if "rel" in kwargs:
            rel = int(kwargs.get("rel"))
        else:
            rel = 1  # default is relativistic

        if "filename" in kwargs:
            filename = kwargs.get("filename")
        else:
            filename = "atorb_{0}.txt".format(ele.symbol)

        if "header" in kwargs:
            header = kwargs.value("header")
        else:
            header = "  atorb input file: {0}.".format(os.path.basename(filename))

        if "method" in kwargs:
            method = str(kwargs.value("method"))
            methods_dict = {"HF": "0.d0", "LDA": "1.d0", "xalpha": "-alpha"}
            if method in methods_dict:
                method = methods_dict.get("method")
            else:
                method = "0.d0"
        else:
            method = "0.d0"

        if "relic" in kwargs:
            relic = float(kwargs.get("relic"))
        else:
            relic = 0

        if "mixing_SCF" in kwargs:
            mixing_SCF = float(kwargs.get("mixing_SCF"))
        else:
            mixing_SCF = 0.5

        if "tolerance" in kwargs:
            eigen_tol = float(kwargs.get("tolerance"))
        else:
            eigen_tol = 0.0005

        if "ech" in kwargs:
            ech = int(kwargs.get("ech"))
        else:
            ech = 100

        orbital_models = []
        for entry in electrons:
            orbital_models.append(
                coerce_model(
                    AtorbElectron,
                    {
                        "n": entry[0],
                        "l": entry[1],
                        "m": entry[2],
                        "j": entry[3],
                        "s": entry[4],
                        "occ": entry[5],
                    },
                )
            )

        atorb_model = coerce_model(
            AtorbInputModel,
            {
                "z": Z,
                "nr": NR,
                "rel": rel,
                "method": method,
                "relic": relic,
                "nlevels": nlevels,
                "mixing_scf": mixing_SCF,
                "eigen_tol": eigen_tol,
                "ech": ech,
                "orbitals": orbital_models,
                "output": output,
                "header": header,
            },
        )
        atorb_model.ensure_valid()

        render_atorb_file(atorb_model, filename=filename)

        return filename  # return output filename for further use

    @staticmethod
    def validate_input_file(input_path):
        """
        Validate a generated or user-supplied atorb input file.

        Parses the input file using the `phaseshifts.validation.atorb` logic
        to ensure it meets the strict formatting and physical requirements of
        the solver.

        Parameters
        ----------
        input_path : str
            Path to the file to validate.

        Returns
        -------
        AtorbInputModel
            Parsed and validated representation of the input file.
        """
        return validate_atorb_file(os.path.abspath(input_path))

    @staticmethod
    def calculate_Q_density(**kwargs):
        """
        Execute the atomic structure calculation to obtain radial charge densities.

        This method serves as the primary driver for the 'atorb' Fortran routine.
        It either generates a new input file based on an element symbol or validates
        an existing input file, then invokes the solver (`hartfock`) via `libphsh`.

        The solver computes the radial wavefunctions :math:`\\phi_{n,l,j}(r)` and
        constructs the total spherical charge density :math:`\\rho(r)`.

        Parameters
        ----------
        **kwargs : dict, optional
            Keyword arguments passed to :meth:`gen_input` if `element` is provided.

            element : str or int
                Target element. If provided, an input file is generated on-the-fly.
            input : str
                Path to an existing atorb input file. Ignored if `element` is provided.
            output_dir : str
                Directory where the output file (containing phase shifts or charge densities)
                should be placed. Defaults to current directory.

        Returns
        -------
        str
            Path to the output file containing the calculated atomic charge densities.

        Raises
        ------
        ValueError
            If neither `element` nor `input` is specified.
        IOError
            If the output directory cannot be created or accessed.

        Examples
        --------
        >>> Atorb.calculate_Q_density(input='atorb_C.txt')
              18.008635    -33.678535
               4.451786    -36.654271
               1.569616    -37.283660
               0.424129    -37.355634
               0.116221    -37.359816
               0.047172    -37.360317
               0.021939    -37.360435
               0.010555    -37.360464
               0.005112    -37.360471
               0.002486    -37.360473
               0.001213    -37.360473
               0.000593    -37.360473
               0.000290    -37.360474
            N L M J S OCC.
            1   0 0  -1/2   1    2.0000        -11.493862
            2   0 0  -1/2   1    2.0000         -0.788618
            2   1 1  -1/2   1    0.6667         -0.133536
            2   1 1  -3/2   1    1.3333         -0.133311
         TOTAL ENERGY =      -37.360474  -1016.638262

        >>> Atorb.calculate_Q_density(element='H')
               0.500007     -0.343752
               0.152392     -0.354939
               0.065889     -0.357254
               0.028751     -0.357644
               0.012732     -0.357703
               0.005743     -0.357711
               0.002641     -0.357712
               0.001236     -0.357713
               0.000587     -0.357713
               0.000282     -0.357713
         N L M J S OCC.
            1   0 0  -1/2   1    1.0000         -0.229756
         TOTAL ENERGY =       -0.357713     -9.733932


        """
        inp = None
        if "input" in kwargs:
            inp = os.path.abspath(kwargs.pop("input"))

        if "element" in kwargs:
            inp = os.path.abspath(
                Atorb.gen_input(element=kwargs.pop("element"), **kwargs)
            )

        current_dir = os.path.curdir
        output_dir = None
        if "output_dir" in kwargs:
            output_dir = kwargs.pop("output_dir")
            if os.path.isdir(output_dir):
                os.chdir(output_dir)
            else:
                try:
                    os.makedirs(output_dir, exist_ok=True)
                    os.chdir(output_dir)
                except (OSError, IOError) as err:
                    raise IOError(
                        "Unable to create output directory due to {!r}".format(err)
                    )  # noqa

        if not inp:
            raise ValueError("Input file not specified")

        validated = validate_atorb_file(inp)

        # do lazy loading due to documentation not needing compiled code
        import phaseshifts.lib.libphsh  # noqa

        phaseshifts.lib.libphsh.hartfock(
            inp
        )  # calculates atomic orbital charge densities for atom

        output_filename = validated.output

        os.chdir(current_dir)  # return to original directory

        return (
            os.path.join(output_dir, output_filename)
            if output_dir is not None
            else output_filename
        )
