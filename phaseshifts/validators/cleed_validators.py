"""Provides CLEED input file validation."""

# cspell:ignore makeiv mkiv

import glob
import re
import os
import sys


class CleedInputFileValidator(object):
    """Class for validation of CLEED input files.

    The CleedInputFileValidator() class provides a method for checking the input files for errors.
    """

    BULK_INPUT_FILE_EXTENSIONS = {".bul", ".bmin", ".bsr"}
    CONTROL_INPUT_FILE_EXTENSIONS = {".ctr"}
    SURFACE_INPUT_FILE_EXTENSIONS = {".inp", ".pmin", ".par"}
    MAKEIV_INPUT_FILE_EXTENSIONS = {".mkiv"}

    BULK_PARAMETER_PATTERN = re.compile(r"(ei|ef|es|ep|vi|vr|lm|m[12]):")

    def _validate_bulk_atoms(self, input_lines):
        line = "-"
        try:
            _bulk_atoms = [
                "".join(line.split(":")[1].split()[:4])
                for line in input_lines
                if line.startswith("pb:")
            ]
        except ValueError:
            raise ValueError("'%s' is not a valid overlayer atom input" % line)  # noqa

    def _validate_minimum_radii(self, input_lines):
        # try to get minimum radii for mufftin input
        line = "-"
        try:
            minimum_radii = [
                "".join(line.split(":")[1].split()[:4])
                for line in input_lines
                if line.startswith("rm:")
            ]
        except ValueError:
            raise ValueError("'%s' is not a valid minimum radius input" % line)  # noqa
        radii_dict = {}
        for radius in minimum_radii:
            tag, value = zip((str, float), radius.split()[:2])
            radii_dict[tag] = value
        return radii_dict

    def _validate_overlayer_atoms(self, input_lines):
        line = "-"
        try:
            _overlayer_atoms = [
                "".join(line.split(":")[1].split()[:4])
                for line in input_lines
                if line.startswith("po:")
            ]
        except ValueError:
            raise ValueError("'%s' is not a valid overlayer atom input" % line)  # noqa

    def _validate_unit_cell(self, input_lines):
        # unit cell
        line = "-"
        try:
            _basis = [
                line
                for line in input_lines
                if line.startswith("a") and int(line[1]) <= 3 and int(line[1]) > 0
            ]
        except ValueError:
            raise ValueError("'%s' is not a valid basis vector" % line)  # noqa

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

        is_cleed = False

        # expand filename string and match to first file found
        try:
            filename = glob.glob(os.path.expanduser(os.path.expandvars(filename)))[0]
        except IndexError:  # no file
            raise OSError("filename '%s' not found" % filename)  # noqa

        # check filename exists
        ext = os.path.splitext(filename)[1]
        if os.path.isfile(filename):
            if ext not in [".bmin", ".bul", ".inp", ".pmin"] and ext != ".i":
                sys.stderr.write(
                    "Warning: %s is not a valid CLEED extension\n" % filename
                )
                sys.stderr.flush()
            else:
                is_cleed = True
        return is_cleed

    @classmethod
    def file_extension_to_filetype(cls, ext):
        """Attempt to match file ``ext`` to list of known CLEED file types."""
        if ext in cls.BULK_INPUT_FILE_EXTENSIONS:  # bulk input file
            filetype = "bulk"
        elif ext in cls.SURFACE_INPUT_FILE_EXTENSIONS:  # surface input file
            filetype = "surface"
        elif ext in cls.CONTROL_INPUT_FILE_EXTENSIONS:  # control file
            filetype = "control"
        elif ext in cls.MAKEIV_INPUT_FILE_EXTENSIONS:  # mkiv file
            filetype = "make_iv"
        else:
            filetype = None
        return filetype

    def guess_filetype_from_contents(self, filepath):
        # type: (str) -> str
        """Guess the filetype from the contents of input file given by ``filepath``."""
        with open(filepath, mode="r", encoding="ascii") as input_file_ptr:
            # get input lines, stripping left whitespace and comments
            lines = "".join(
                [line.lstrip() for line in input_file_ptr if not line.lstrip().startswith("#")]
            )
            # check for control parameters
            if all(x in lines for x in ["ef=", "ti=", "id=", "wt="]):
                # => probably control file
                filetype = "control"
            elif "pb:" in lines and any(map(self.BULK_PARAMETER_PATTERN.search, lines)):
                # => probably bulk file
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
                raise OSError("unknown CLEED format for file '%s'" % filepath)
        return filetype

    def validate(self, filename, aoi=False):  # pylint: disable=unused-argument
        r"""
        Checks CLEED input file for errors

        Parameters
        ----------
        filename : str
            Path to input file. Should be \*.bul , \*.ctr, \*.inp or \*.bmin
        aoi : bool
            Check for angle of incidence parameters

        Raises
        ------
        ValueError
            If file is not a valid CLEED input file
        """
        filename = glob.glob(os.path.expanduser(os.path.expandvars(filename)))[0]

        basename, ext = os.path.splitext(filename)
        if isinstance(ext, int):
            _aoi = True  # angle of incidence
            basename, ext = os.path.splitext(basename)

        filetype = self.file_extension_to_filetype(ext) or self.guess_filetype_from_contents(filename)

        with open(filename, mode="r", encoding="ascii") as input_file_ptr:
            # get input lines, stripping left whitespace and comments
            input_lines = input_file_ptr.readlines()

        with open(filename, mode="w", encoding="ascii") as output_file_ptr:
            for i, line in enumerate(input_lines):
                line = line.lstrip().rstrip()
                if not line.startswith("#") or line.startswith(""):
                    # process line based on file type
                    if line.startswith("po:"):
                        if filetype not in ["surface", "parameter"]:
                            # comment out line
                            sys.stderr.write("warning: commented out line %i\n" % i)
                            sys.stderr.flush()
                            line = "#" + line

                    else:  # line unknown - comment out
                        sys.stderr.write("warning: commented out line %i\n" % i)
                        sys.stderr.flush()
                        line = "#" + line
                    if filetype in ("bulk", "surface"):
                        pass
                output_file_ptr.write(line + "\r\n")  # playing it safe with line endings
                sys.stderr.flush()

        self._validate_unit_cell(input_lines)
        self._validate_minimum_radii(input_lines)
        self._validate_bulk_atoms(input_lines)
        self._validate_overlayer_atoms(input_lines)
