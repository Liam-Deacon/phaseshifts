#!/usr/bin/env python

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2014 Liam Deacon                                  #
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
Provides CLEED validator and Converter classes.

The CLEEDInputValidator() class provides a method for checking
the input files for errors, whereas the Converter.import_CLEED()
method allows importing CLEED input files as a MTZ_model class

"""

import json
import os
import sys
import math
import tempfile
import re
import warnings
from glob import glob

from phaseshifts import model
from phaseshifts.elements import ELEMENTS

try:  # optional validation
    import jsonschema  # type: ignore
except ImportError:
    jsonschema = None


class CLEEDInputValidator:
    """
    Class for validation of CLEED input files
    """

    @staticmethod
    def is_cleed_file(filename):
        """
        Determine if file is a CLEED input file

        Returns
        -------
        True
            if a valid filename ending in any of .bul, .inp, .bsr, .bmin

        False
            otherwise
        """

        # expand filename string and match to first file found
        try:
            filename = glob(os.path.expanduser(os.path.expandvars(filename)))[0]
        except IndexError:  # no file
            raise OSError("filename '%s' not found" % filename)

        # check filename exists
        if not os.path.isfile(filename):
            return False

        # check extension is valid
        ext = os.path.splitext(filename)[1]
        if ext not in [".bmin", ".bul", ".inp", ".pmin"]:
            if ext != ".i":
                sys.stderr.write(
                    "Warning: %s is not a valid CLEED extension\n" % filename
                )
                sys.stderr.flush()
            return False
        else:
            return True

    def validate(self, filename, aoi=False):
        r"""
        Checks CLEED input file for errors

        Parameters
        ----------
        filename : str
            Path to input file. Should be \*.bul , \*.ctr, \*.inp or \*.bmin
        aoi : bool
            Check for angle of incidence parameters

        """
        filename = glob(os.path.expanduser(os.path.expandvars(filename)))[0]

        basename, ext = os.path.splitext(filename)
        if isinstance(ext, int):
            aoi = True
            basename, ext = os.path.splitext(basename)

        BULK_INPUT_FILE_EXTENSIONS = {".bul", ".bmin", ".bsr"}
        CONTROL_INPUT_FILE_EXTENSIONS = {".ctr"}
        SURFACE_INPUT_FILE_EXTENSIONS = {".inp", ".pmin", ".par"}
        MAKEIV_INPUT_FILE_EXTENSIONS = {".mkiv"}

        if ext in BULK_INPUT_FILE_EXTENSIONS:  # bulk input file
            filetype = "bulk"
        elif ext in SURFACE_INPUT_FILE_EXTENSIONS:  # surface input file
            filetype = "surface"
        elif ext in CONTROL_INPUT_FILE_EXTENSIONS:  # control file
            filetype = "control"
        elif ext in MAKEIV_INPUT_FILE_EXTENSIONS:  # mkiv file
            filetype = "make_iv"
        else:  # try to determine filetype from file data
            with open(filename) as f:
                # get input lines, stripping left whitespace and comments
                lines = "".join(
                    [line.lstrip() for line in f if not line.lstrip().startswith("#")]
                )
                # check for control parameters
                if (
                    "ef=" in lines
                    and "ti=" in lines
                    and "id=" in lines
                    and "wt=" in lines
                ):
                    # probably control file
                    filetype = "control"
                elif "pb:" in lines and (
                    "ei:" in lines
                    or "ef:" in lines
                    or "es:" in lines
                    or "ep:" in lines
                    or "vi:" in lines
                    or "vr:" in lines
                    or "lm:" in lines
                    or "m1:" in lines
                    or "m2:" in lines
                ):  # probably bulk file
                    filetype = "bulk"
                elif "po:" in lines and (
                    "a1:" in lines or "a2:" in lines or "m1:" in lines or "m2:" in lines
                ):  # surface file?
                    filetype = "surface"
                elif "po:" in lines and not (
                    "a1:" in lines or "a2:" in lines or "m1:" in lines or "m2:" in lines
                ):
                    filetype = "parameter"
                else:
                    raise OSError("unknown CLEED format for file '%s'" % filename)

        with open(filename) as f:
            # get input lines, stripping left whitespace and comments
            lines = [line for line in f]

        with open(filename, "w") as f:
            for i, line in enumerate(lines):
                line = line.lstrip().rstrip()
                if not line.startswith("#") or line.startswith(""):
                    # process line based on file type
                    if line.startswith("po:"):
                        if filetype not in ["surface", "parameter"]:
                            # comment out line
                            sys.stderr.write("warning: commented out line " "%i\n" % i)
                            sys.stderr.flush()
                            line = "#" + line

                    else:  # line unknown - comment out
                        sys.stderr.write("warning: commented out line " "%i\n" % i)
                        sys.stderr.flush()
                        line = "#" + line
                    if filetype == "bulk" or filetype == "surface":
                        pass
                f.write(line + "\r\n")  # playing it safe with line endings
                sys.stderr.flush()
        # unit cell
        try:
            basis = [
                line
                for line in lines
                if line.startswith("a") and int(line[1]) <= 3 and int(line[1]) > 0
            ]

        except ValueError:
            raise ValueError("'%s' is not a valid basis vector" % line[:2])

        # try to get minimum radii for mufftin input
        try:
            minimum_radii = [
                "".join(line.split(":")[1].split()[:4])
                for line in lines
                if line.startswith("rm:")
            ]
        except ValueError:
            raise ValueError("'%s' is not a valid minimum radius input" % line[:2])
        radii_dict = {}
        for radius in minimum_radii:
            tag, value = zip((str, float), radius.split()[:2])
            radii_dict[tag] = value

        try:
            overlayer_atoms = [
                "".join(line.split(":")[1].split()[:4])
                for line in lines
                if line.startswith("po:")
            ]
        except ValueError:
            raise ValueError("'%s' is not a valid overlayer atom input" % line[:2])

        try:
            bulk_atoms = [
                "".join(line.split(":")[1].split()[:4])
                for line in lines
                if line.startswith("pb:")
            ]
        except ValueError:
            raise ValueError("'%s' is not a valid overlayer atom input" % line[:2])


_CLEEDPY_SCHEMA = {
    "type": "object",
    "properties": {
        "system_name": {"type": "string"},
        "unit_cell": {
            "type": "object",
            "required": ["a1", "a2", "a3"],
            "properties": {
                "a1": {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 3,
                    "maxItems": 3,
                },
                "a2": {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 3,
                    "maxItems": 3,
                },
                "a3": {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 3,
                    "maxItems": 3,
                },
            },
        },
        "superstructure_matrix": {
            "type": "object",
            "properties": {
                "m1": {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 2,
                    "maxItems": 2,
                },
                "m2": {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 2,
                    "maxItems": 2,
                },
            },
        },
        "overlayers": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["phase_file", "position"],
                "properties": {
                    "phase_file": {"type": "string"},
                    "position": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                    },
                },
            },
        },
        "bulk_layers": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["phase_file", "position"],
                "properties": {
                    "phase_file": {"type": "string"},
                    "position": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                    },
                },
            },
        },
        "minimum_radius": {
            "type": "object",
            "additionalProperties": {"type": "number"},
        },
        "optical_potential": {"type": "array", "items": {"type": "number"}},
        "energy_range": {
            "type": "object",
            "properties": {
                "initial": {"type": "number"},
                "final": {"type": "number"},
                "step": {"type": "number"},
            },
        },
        "maximum_angular_momentum": {"type": "integer"},
    },
    "required": ["unit_cell", "bulk_layers", "overlayers"],
}


class Converter:
    """
    Convert different input into phaseshift compatible input

    """

    @staticmethod
    def import_CLEED(filename, **kwargs):
        """
        Imports CLEED input file and converts model to muffin-tin input.

        It assumes the following:
          * the basis vectors :math:`a_1`, :math:`a_2`, & :math:`a_3` are
            :math:`x`,:math:`y`,:math:`z` cartezian coordinates
          * if no :math:`a_3` is found, the maximum :math:`z` distance between
            atoms multiplied by four is given.
          * the unitcell is converted from cartezian to fractional coordinates
          * atom sites are converted from Angstrom to Bohr units
          * additional info from the phase shift filename is provided by
            splitting the '_' chars:

            1. First string segment is element or symbol, e.g. Ni
            2. Second string segment is the oxidation (valence), e.g. +2

          * lines with :code:`rm:` provide the radii dictionary of the
            individual phase shifts, whereas lines starting with
            :code:`lmax:` provide the lmax dictionary for individual phase
            shifts. Additionally, the valency can be given in lines
            starting with :code:`ox:`. Both bulk and surface input files
            will be searched for these.
          * if no :code:`rm:` found for that species, the atomic radius is used
            for zero valence, otherwise the covalent radius is used.
          * if no lmax values are found for the specific phase shift, then the
            global value will be used instead.
          * if no valence values are found for the specific phase shift, then
            the guessed oxidation value from the phase shift filename is used
            instead. However, if the oxidation state cannot be parsed using
            the filename then a default oxidation state of zero is used.

        Additional information can, however, be provided using 'phs:' at the
        start of a line within the input file and may have the following
        formats:

          1. "**phs:** c *<float>* nh *<int>* nform *<int>* exchange *<float>*"
          2. "**phs:** *<phase_shift>* valence *<float>* radius *<float>*"

        The identifiers ``exchange``, ``nform``, ``valence`` and ``radius`` may
        be abbreviated to ``exc``, ``nf``, ``val`` and ``rad``, respectively.

        Information given in this manner overrides any other input.


        Parameters
        ----------
        filename : str
            Path to input file.

        Returns
        -------
        phaseshifts.model.MTZ_model

        Raises
        ------
        IOError : filename invalid
        ValueError : bad formatting of input

        """

        try:
            verbose = kwargs["verbose"]
        except KeyError:
            verbose = False

        if not isinstance(filename, str):
            # careful: passing file object not string
            raise OSError("'%s' is not a string" % filename)

        filename = glob(os.path.expanduser(os.path.expandvars(filename)))[0]

        if os.path.isdir(filename):
            # check directory for unique .bul or .bmin files
            raise OSError("'%s' is a directory, not a file" % filename)

        elif os.path.isfile(filename):
            # determine type of file
            name, ext = os.path.splitext(filename)

            if ext in [".bul", ".bmin"]:  # bulk input file
                other_input = name + ".inp"
            elif ext in [".inp", ".pmin", ".par"]:  # surface input file
                other_input = name + ".bul"

            try:
                with open(other_input) as f:
                    other_lines = [line for line in f]

            except OSError:
                other_lines = []

        else:
            raise OSError("cannot open '%s'" % filename)

        with open(filename) as f:
            # get input lines, stripping left whitespace and comments
            lines = [line.lstrip() for line in f if not line.lstrip().startswith("#")]

        # initialise model
        atoms = []
        radii_dict = {}
        lmax_dict = {}
        oxidation_dict = {}
        unitcell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        z_dist = 4.0  # calculate 'c' for slab from atom coordinates
        z_min = sys.float_info.max
        z_max = sys.float_info.min
        title = "/".join(
            [line.split(":")[1].lstrip() for line in lines if line.startswith("c:")]
        )

        # Hartree-Fock exchange term alpha
        try:
            alpha = float(
                [
                    line.split(":")[1]
                    for line in lines
                    if line.startswith("alpha:")
                    or line.startswith("exchange:")
                    or line.startswith("exc:")
                ][0]
            )
        except (ValueError, IndexError):
            alpha = None

        # unit cell
        basis = [
            line
            for line in lines
            if line.startswith("a") and int(line[1]) <= 3 and int(line[1]) > 0
        ]

        a1 = a2 = a3 = [0.0, 0.0, 0.0]  # intialise basis vectors
        for vector_line in basis:
            line = " ".join(vector_line.replace(":", "").split()[:4])

            try:
                (vector_str, x, y, z) = (
                    t(s) for (t, s) in zip((str, float, float, float), line.split())
                )
            except ValueError:
                raise ValueError(
                    "'%s' is not a valid basis vector" % vector_line.replace("\n", "")
                )

            exec("{} = [{:f}, {:f}, {:f}]".format(vector_str, x, y, z))

        zmin = sys.float_info.max
        a = zmax = sys.float_info.min
        for vector in [a1, a2, a3]:
            # compare magnitude of xy components
            len_a = math.sqrt(math.pow(vector[0], 2) + math.pow(vector[1], 2))
            if len_a > a:
                a = len_a

            # find min and max z
            if vector[2] < zmin:
                zmin = vector[2]

            if vector[2] > zmax:
                zmax = vector[2]

        # try to get minimum radii for mufftin input
        try:
            minimum_radii = [
                " ".join(line.split(":")[1].split()[:4])
                for line in lines + other_lines
                if line.startswith("rm:")
            ]

            for radius in minimum_radii:
                tag, value = zip((str, float), radius.split()[:2])
                radii_dict[tag] = value

        except ValueError:
            raise ValueError("'%s' is not a valid minimum radius input" % line[:2])

        # try to get lmax values for each phaseshift type
        try:
            lmax_values = [
                " ".join(line.split(":")[1].split()[:4])
                for line in lines + other_lines
                if line.startswith("lmax:")
            ]

            for lmax_str in lmax_values:
                tag, value = (t(s) for t, s in zip((str, int), lmax_str.split()[:2]))
                lmax_dict[tag] = value

        except ValueError:
            raise ValueError("'%s' is not a valid lmax input" % line[:2])

        # try to get oxidation values for each phaseshift type
        try:
            oxidation_values = [
                " ".join(line.split(":")[1].split()[:4])
                for line in lines + other_lines
                if line.startswith("val:")
            ]

            for oxi_str in oxidation_values:
                tag, value = (t(s) for t, s in zip((str, float), oxi_str.split()[:2]))
                oxidation_dict[tag] = value

        except ValueError:
            raise ValueError("'%s' is not a valid valency (oxidation state)" % line[:2])

        # try to get the overlayer atoms
        overlayer_atoms = [
            " ".join(line.split(":")[1].split()[:4])
            for line in lines
            if line.startswith("po:")
        ]

        # extract details for each overlayer atom
        for overlayer_atom in overlayer_atoms:
            try:
                id_str, x, y, z = (
                    t(s)
                    for t, s in zip((str, float, float, float), overlayer_atom.split())
                )

            except ValueError:
                raise ValueError("'%s' line input is invalid")

            # extract further information from id string
            id_str = os.path.basename(
                os.path.expanduser(os.path.expandvars(id_str))
            )  # ensure path not included
            element = id_str.split("_")[0]  # assume element name / symbol
            element = "".join([ch for ch in element if ch.isalpha()])

            if element not in ELEMENTS:
                raise NameError(
                    "Unknown element '%s' from phase shift name"
                    "'%s'" % (element, id_str)
                )

            # try to determine oxidation state from string
            oxidation = None
            if id_str in oxidation_dict:
                oxidation = oxidation_dict[id_str]
            else:
                try:
                    oxidation = str(id_str.split("_")[1])  # assume valence 2nd
                    if oxidation.startswith("-") or oxidation.startswith("+"):
                        oxidation = float(oxidation)
                    elif oxidation.endswith("-") or oxidation.endswith("+"):
                        if oxidation.endswith("-"):
                            oxidation = float(oxidation.rstrip("-"))
                        else:
                            oxidation = float(oxidation.rstrip("+"))
                    else:
                        oxidation = 0.0  # oxidation not specified
                except (IndexError, ValueError):  # invalid oxidation substring
                    pass
                finally:
                    if oxidation is None:
                        oxidation = 0.0  # default oxidation if not specified

            rm = None
            try:
                rm = radii_dict[id_str]
            except KeyError:
                pass

            lmax = None
            try:
                lmax = lmax_dict[id_str]
            except KeyError:
                pass

            if lmax is None:
                try:  # to get lmax from id string
                    substrings = [
                        s for s in id_str.lower().split("_") if s.startswith("lmax")
                    ]
                    if substrings:
                        lmax = int("".join([s for s in substrings[0] if s.isdigit()]))
                except ValueError:
                    pass

            if rm and lmax:
                atom = model.Atom(
                    element,
                    [x, y, z],
                    tag=id_str,
                    valence=oxidation,
                    radius=rm,
                    lmax=lmax,
                )
            elif rm and not lmax:
                atom = model.Atom(
                    element, [x, y, z], tag=id_str, valence=oxidation, radius=rm
                )
            elif not rm and lmax:
                atom = model.Atom(
                    element, [x, y, z], tag=id_str, valence=oxidation, lmax=lmax
                )
            else:
                atom = model.Atom(element, [x, y, z], tag=id_str, valence=oxidation)

            # update z limits
            if z < z_min:
                z_min = z
            if z > z_max:
                z_max = z

            # add atom to list
            atoms.append(atom)

        # try to get the bulk atoms
        bulk_atoms = [
            " ".join(line.split(":")[1].split()[:4])
            for line in lines
            if line.startswith("pb:")
        ]
        if verbose:
            print("CLEED model")
            print("\tbulk atoms: %s" % [s for s in bulk_atoms])

        # extract details for each bulk atom
        for bulk_atom in bulk_atoms:
            # split line data
            try:
                (id_str, x, y, z) = (
                    t(s) for t, s in zip((str, float, float, float), bulk_atom.split())
                )
                # extract further information from id string
                id_str = os.path.basename(
                    os.path.expanduser(os.path.expandvars(id_str))
                )  # ensure path not included
                element = id_str.split("_")[0]  # assume element name / symbol
                element = "".join([ch for ch in element if ch.isalpha()])

                if element not in ELEMENTS:
                    raise NameError(
                        "Unknown element '%s' from phase shift name"
                        "'%s'" % (element, id_str)
                    )

            except ValueError:
                print("'%s' line input is invalid")
                raise ValueError("'%s' line input is invalid")

            # try to determine oxidation state from string
            oxidation = None
            if id_str in oxidation_dict:
                oxidation = oxidation_dict[id_str]
            else:
                try:
                    oxidation = str(id_str.split("_")[1])  # assume valence 2nd
                    if oxidation.startswith("-") or oxidation.startswith("+"):
                        oxidation = float(oxidation)
                    elif oxidation.endswith("-") or oxidation.endswith("+"):
                        if oxidation.endswith("-"):
                            oxidation = float(oxidation.rstrip("-"))
                        else:
                            oxidation = float(oxidation.rstrip("+"))
                    else:
                        oxidation = 0.0  # oxidation not specified
                except (IndexError, ValueError):  # invalid oxidation substring
                    pass
                finally:
                    if oxidation is None:
                        oxidation = 0.0  # default oxidation if not specified

            rm = None
            if id_str in radii_dict:
                rm = radii_dict[id_str]

            lmax = None
            if id_str in lmax_dict:
                lmax = lmax_dict[id_str]

            if rm and lmax:
                atom = model.Atom(
                    element,
                    [x, y, z],
                    tag=id_str,
                    valence=oxidation,
                    radius=rm,
                    lmax=lmax,
                )
            elif rm and not lmax:
                atom = model.Atom(
                    element, [x, y, z], tag=id_str, valence=oxidation, radius=rm
                )
            elif not rm and lmax:
                atom = model.Atom(
                    element, [x, y, z], tag=id_str, valence=oxidation, lmax=lmax
                )
            else:
                atom = model.Atom(element, [x, y, z], tag=id_str, valence=oxidation)

            # update z limits
            if z < z_min:
                z_min = z
            if z > z_max:
                z_max = z

            # add atom to list
            atoms.append(atom)

        # change basis to fractional coordinates (SPA units)
        # i.e. in terms of magnitude 'a'
        for vector in [a1, a2, a3]:
            vector = [x / a for x in vector]

        if a3 == [0.0, 0.0, 0.0]:  # no a3 in input file => slab
            z_dist = (z_max - z_min) / float(a) / 0.529
            c = z_dist * 5  # slab calculation
            a3[2] = c
        else:
            c = a3[2] / float(a) / 0.529

        # create unitcell
        unitcell = model.Unitcell(a, c, [a1, a2, a3])

        mtz_model = model.MTZ_model(unitcell, atoms)
        mtz_model.name = title

        if alpha:
            mtz_model.set_exchange(alpha)

        # read additional phaseshift statement lines 'phs:'
        phs_lines = [
            line.lower().lstrip("phs:").rstrip("#")
            for line in lines
            if line.lower().startswith("phs:")
        ]

        for line in phs_lines:
            strings = line.split()
            phase_shift = strings[0]
            phs_atoms = [
                atom.tag for atom in mtz_model.atoms if atom.tag == phase_shift
            ]
            if phase_shift not in phs_atoms:
                for i in range(len(strings)):
                    try:
                        if strings[i] == "c":
                            c = strings[i] / float(a) / 0.529
                            mtz_model.unitcell.set_c(c)

                        elif strings[i] == "nh":
                            mtz_model.set_nh(strings[i + 1])

                        elif strings[i] == "exc" or strings[i] == "exchange":
                            mtz_model.set_exchange(strings[i + 1])

                        elif strings[i] == "nf" or strings[i] == "nform":
                            mtz_model.set_nform(strings[i + 1])

                    except IndexError:
                        break
            else:
                # apply commands to atom
                for i in range(1, len(strings)):
                    try:
                        if strings[i] == "val" or strings[i] == "valence":
                            for atom in phs_atoms:
                                atom.valence = float(strings[i + 1])

                        elif strings[i] == "rad" or strings[i] == "radius":
                            for atom in phs_atoms:
                                atom.radius = float(strings[i + 1])

                    except IndexError:
                        break

        return mtz_model  # converted muffin-tin potential model

    @staticmethod
    def _load_structured_input(filename, yaml_loader=None):
        """Load structured input (JSON or YAML)."""
        ext = os.path.splitext(filename)[1].lower()
        with open(filename, "r") as f:
            content = f.read()

        if ext == ".json":
            return json.loads(content)

        loader_func = None
        if yaml_loader is not None:
            loader_func = yaml_loader
        else:
            try:
                import yaml  # type: ignore

                loader_func = yaml.safe_load
            except ImportError:
                loader_func = None

        if loader_func:
            return loader_func(content)

        # fallback: try JSON parsing even if extension not json
        try:
            return json.loads(content)
        except ValueError:
            raise ImportError(
                "PyYAML is required to read YAML input. Install pyyaml or provide JSON."
            )

    @staticmethod
    def _validate_structured_input(data):
        """Validate cleedpy-style input if jsonschema is available."""
        if jsonschema is None:
            warnings.warn(
                "jsonschema not installed; skipping structured input validation.",
                RuntimeWarning,
            )
            return
        try:
            jsonschema.validate(instance=data, schema=_CLEEDPY_SCHEMA)
        except jsonschema.exceptions.ValidationError as exc:  # type: ignore[attr-defined]
            raise ValueError(f"Invalid input structure: {exc.message}") from exc

    @staticmethod
    def _resolve_element(tag):
        """Extract element symbol from cleedpy phase_file."""
        if tag is None:
            raise NameError("Phase file tag is missing.")

        basename = os.path.splitext(os.path.basename(str(tag)))[0]
        for part in [p for p in re.split(r"[._-]", basename) if p]:
            letters = []
            for ch in part:
                if ch.isalpha():
                    letters.append(ch)
                else:
                    break
            if not letters:
                continue
            symbol = (letters[0].upper() + "".join(letters[1:]).lower())[:2]
            if symbol in ELEMENTS:
                return symbol

        raise NameError("Unknown element parsed from '%s'" % (tag,))

    @staticmethod
    def _vector_length(vec):
        return math.sqrt(sum(v**2 for v in vec))

    @staticmethod
    def _apply_superstructure(a1, a2, m1, m2):
        base_a1 = list(a1)
        base_a2 = list(a2)
        new_a1 = [m1[0] * base_a1[i] + m1[1] * base_a2[i] for i in range(3)]
        new_a2 = [m2[0] * base_a1[i] + m2[1] * base_a2[i] for i in range(3)]
        return new_a1, new_a2

    @staticmethod
    def _build_unitcell(data, superstructure):
        a1 = data["unit_cell"]["a1"]
        a2 = data["unit_cell"]["a2"]
        a3 = data["unit_cell"]["a3"]

        if superstructure and "superstructure_matrix" in data:
            m1 = data["superstructure_matrix"].get("m1")
            m2 = data["superstructure_matrix"].get("m2")
            if m1 and m2:
                a1, a2 = Converter._apply_superstructure(a1, a2, m1, m2)

        a_mag = max(
            Converter._vector_length(a1[:2] + [0]),
            Converter._vector_length(a2[:2] + [0]),
        )
        c_mag = Converter._vector_length(a3)
        return model.Unitcell(a_mag, c_mag, [a1, a2, a3])

    @staticmethod
    def _build_atoms(layers, minimum_radius, lmax_global):
        atoms = []
        for layer in layers:
            element = Converter._resolve_element(layer.get("phase_file"))
            radius = minimum_radius.get(element, None)
            position = layer.get("position") or [0.0, 0.0, 0.0]
            atom = model.Atom(
                element, coordinates=position, tag=layer.get("phase_file")
            )
            if radius:
                atom.set_mufftin_radius(radius)
            if lmax_global:
                atom.lmax = int(lmax_global)
            atoms.append(atom)
        return atoms

    @staticmethod
    def import_cleedpy_input(filename, superstructure=True, yaml_loader=None, **kwargs):
        """
        Parse a cleedpy-style structured input file (JSON or YAML) into bulk and slab MTZ models.

        Parameters
        ----------
        filename : str
            Path to cleedpy YAML input.
        superstructure : bool, optional
            If True, apply the superstructure_matrix (m1/m2) to the in-plane
            basis vectors to construct the surface cell. Defaults to True.

        Returns
        -------
        tuple(model.MTZ_model, model.MTZ_model, dict)
            (bulk_model, slab_model, metadata)

        Notes
        -----
        This provides a bridge from cleedpy's YAML geometry to the Barbieri/Van
        Hove phase-shift workflow by constructing equivalent muffin-tin models.
        """
        data = Converter._load_structured_input(filename, yaml_loader=yaml_loader)
        if not isinstance(data, dict):
            raise ValueError("Structured input must be a mapping/object.")
        Converter._validate_structured_input(data)
        unitcell = Converter._build_unitcell(data, superstructure)

        minimum_radius = {
            str(k).title(): float(v) for k, v in data.get("minimum_radius", {}).items()
        }
        lmax_global = data.get("maximum_angular_momentum")

        bulk_atoms = Converter._build_atoms(
            data.get("bulk_layers", []), minimum_radius, lmax_global
        )
        slab_atoms = Converter._build_atoms(
            data.get("overlayers", []) + data.get("bulk_layers", []),
            minimum_radius,
            lmax_global,
        )

        bulk_model = model.MTZ_model(unitcell, bulk_atoms)
        slab_model = model.MTZ_model(unitcell, slab_atoms)

        metadata = {
            "system_name": data.get("system_name"),
            "optical_potential": data.get("optical_potential"),
            "energy_range": data.get("energy_range"),
            "maximum_angular_momentum": lmax_global,
        }

        return bulk_model, slab_model, metadata

    @staticmethod
    def cleedpy_to_inputs(filename, tmp_dir=None, yaml_loader=None, **kwargs):
        """
        Convert cleedpy YAML to temporary bulk/slab .i files for phase-shift generation.

        Returns
        -------
        tuple(str, str, dict)
            Paths to generated bulk and slab input files, plus parsed metadata.
        """
        if tmp_dir is None:
            tmp_dir = tempfile.gettempdir()
        if not os.path.exists(tmp_dir):
            try:
                os.makedirs(tmp_dir)
            except OSError:
                pass

        bulk_model, slab_model, metadata = Converter.import_cleedpy_input(
            filename, yaml_loader=yaml_loader, **kwargs
        )

        system_name = (
            metadata.get("system_name")
            or os.path.splitext(os.path.basename(filename))[0]
        )
        bulk_file = os.path.join(tmp_dir, "%s_bulk.i" % system_name)
        slab_file = os.path.join(tmp_dir, "%s_slab.i" % system_name)

        header = "Generated from cleedpy YAML: %s" % os.path.basename(filename)
        bulk_model.gen_input(filename=bulk_file, header=header, bulk=True)
        slab_model.gen_input(filename=slab_file, header=header, bulk=False)

        return bulk_file, slab_file, metadata


class CSearch:
    """class for csearch related data exchange"""

    def __init__(self, model_name, leed_command=None):
        self.setModel(model_name)
        self._getResults()

    def setModel(self, name):
        """sets the model name"""
        self.model = str(name)

    def getIteration(self, iteration):
        try:
            with open(self.model + ".log") as f:
                return [line for line in f if line.startswith("#") and "par" in line][
                    int(iteration)
                ]
        except (IndexError, ValueError, OSError):
            return None

    def getRFactor(self, iteration):
        pass

    def getTimeStamp(self, iteration):
        pass

    def _getResults(self):
        pass
