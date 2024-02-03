"""Module defining adaptors to CLEED."""
import io
import glob
import os.path
import sys
from math import sqrt

from phaseshifts.core import model
from phaseshifts.core import elements


class CleedInputFileToMuffinTinPotentialModelConverter:
    """Convert different input into phase shift compatible input."""

    @staticmethod
    def import_cleed(filename, **kwargs):
        """
        Imports CLEED input file and converts model to muffin-tin input.

        It assumes the following:
          * the basis vectors :math:`a_1`, :math:`a_2`, & :math:`a_3` are
            :math:`x`,:math:`y`,:math:`z` cartesian coordinates
          * if no :math:`a_3` is found, the maximum :math:`z` distance between
            atoms multiplied by four is given.
          * the unitcell is converted from cartesian to fractional coordinates
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

        filename = glob.glob(os.path.expanduser(os.path.expandvars(filename)))[0]

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
                with io.open(other_input, mode="r", encoding="ascii") as input_file_ptr:
                    other_lines = input_file_ptr.readlines()

            except OSError:
                other_lines = []

        else:
            raise OSError("cannot open '%s'" % filename)

        with io.open(filename, mode="r", encoding="ascii") as input_file_ptr:
            # get input lines, stripping left whitespace and comments
            lines = [line.lstrip() for line in input_file_ptr if not line.lstrip().startswith("#")]

        # initialise model
        atoms = []
        radii_dict = {}
        lmax_dict = {}
        oxidation_dict = {}
        unitcell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        z_dist = 4.0  # calculate 'c' for slab from atom coordinates
        z_min = sys.float_info.max
        z_max = sys.float_info.min

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

        a1 = a2 = a3 = [0.0, 0.0, 0.0]  # initialise basis vectors
        for vector_line in basis:
            line = " ".join(vector_line.replace(":", "").split()[:4])

            try:
                (vector_str, x, y, z) = (
                    t(s) for (t, s) in zip((str, float, float, float), line.split())
                )
            except ValueError:
                raise ValueError(  # noqa
                    "'%s' is not a valid basis vector" % vector_line.replace("\n", "")
                )  # noqa

            exec("{} = [{:f}, {:f}, {:f}]".format(vector_str, x, y, z))  # nosec: B102  # pylint: disable=exec-used

        zmin = sys.float_info.max
        a = zmax = sys.float_info.min
        for vector in [a1, a2, a3]:
            # compare magnitude of xy components
            len_a = sqrt(pow(vector[0], 2) + pow(vector[1], 2))
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
            raise ValueError("'%s' is not a valid minimum radius input" % line)

        # try to get lmax values for each phase shift type
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
            raise ValueError("'%s' is not a valid lmax input" % line)

        # try to get oxidation values for each phase shift type
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
            raise ValueError("'%s' is not a valid valency (oxidation state)" % line)

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
                raise ValueError("'%s' line input is invalid" % overlayer_atom)

            # extract further information from id string
            id_str = os.path.basename(
                os.path.expanduser(os.path.expandvars(id_str))
            )  # ensure path not included
            element = id_str.split("_")[0]  # assume element name / symbol
            element = "".join([ch for ch in element if ch.isalpha()])

            if element not in elements.ELEMENTS:
                raise NameError(
                    "Unknown element '%s' from phase shift name"
                    "'%s'" % (element, id_str)
                )

            # try to determine oxidation state from string
            oxidation = oxidation_dict.get(id_str)  # type: float | None
            if oxidation is None:
                try:
                    oxidation_str = str(id_str.split("_")[1])  # assume valence 2nd
                    if oxidation_str.startswith("-") or oxidation_str.startswith("+"):
                        oxidation = float(oxidation_str)
                    elif oxidation_str.endswith("-") or oxidation_str.endswith("+"):
                        if oxidation_str.endswith("-"):
                            oxidation = float(oxidation_str.rstrip("-"))
                        else:
                            oxidation = float(oxidation_str.rstrip("+"))
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

                if element not in elements.ELEMENTS:
                    raise NameError(
                        "Unknown element '%s' from phase shift name"
                        "'%s'" % (element, id_str)
                    )

            except ValueError:
                raise ValueError("'%s' line input is invalid")

            # try to determine oxidation state from string
            oxidation = None
            if id_str in oxidation_dict:
                oxidation = oxidation_dict[id_str]
            else:
                try:
                    oxidation_str = str(id_str.split("_")[1])  # assume valence 2nd
                    if oxidation_str.startswith("-") or oxidation_str.startswith("+"):
                        oxidation = float(oxidation_str)
                    elif oxidation_str.endswith("-") or oxidation_str.endswith("+"):
                        if oxidation_str.endswith("-"):
                            oxidation = float(oxidation_str.rstrip("-"))
                        else:
                            oxidation = float(oxidation_str.rstrip("+"))
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

        if alpha:
            mtz_model.set_exchange(alpha)

        # read additional phase shift statement lines 'phs:'
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