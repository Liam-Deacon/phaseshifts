"""A collection of utilities to assist with extracting phshift2007 in a cross-platform manner."""

import argparse
import io
import os.path
import sys
import subprocess  # nosec
import urllib.request
import zipfile

PHSHIFT2007_DOWNLOAD_URL = (
    "https://www.icts.hkbu.edu.hk/VanHove_files/leed/phshift2007.zip"
)
DEFAULT_DIRPATH = "."


def download_phshift2007(dirpath=DEFAULT_DIRPATH, force=False):
    # type: (str, bool) -> str
    """Download phshift2007 from the web."""
    if not os.path.exists(dirpath):
        os.makedirs(dirpath, **({"exist_ok": True} if sys.version_info >= (3, 2) else {}))
    zipfile_path = os.path.join(dirpath, "phshift2007.zip")
    if not os.path.exists(zipfile_path) or force:
        urllib.request.urlretrieve(PHSHIFT2007_DOWNLOAD_URL, zipfile_path)  # nosec: B310
    return zipfile_path


def unzip_phshift2007(zipfile_path=None, dirpath=DEFAULT_DIRPATH):
    # type: (str|None, str) -> None
    """Unzip phshift2007 contents from the zipfile."""
    zipfile_path = zipfile_path or os.path.join(dirpath, "phshift2007.zip")
    with zipfile.ZipFile(zipfile_path, mode="r") as zip_ref:
        zip_ref.extractall(dirpath)


def _generate_notice_lines():
    """Generate the notice lines for the extracted Fortran programs."""
    return [
        "!* -------------------------------------------------------------------\n"
        "! > **IMPORTANT:** This file was automatically extracted from the\n",
        "! original phshift2007.zip archive using the phaseshifts python\n",
        "! package and has been modified for processing with the fortran\n",
        "! documentation tool, FORD.\n",
        "!\n",
        "! It is intended for documentation only. Please obtain copies of the\n",
        "! origin source from the _LEED Calculation Home Page_ site at: \n",
        "! https://www.icts.hkbu.edu.hk/VanHove_files/leed/phshift2007.zip\n",
        "!\n",
        "!---------------------------------------------------------------------\n",
        "\n",
    ]

def _extract_fortran_programs_from_ab_file(ab_filepath, output_dirpath=None):  # pylint: disable=too-many-locals
    """Extracts the Fortran programs from the .ab files in the phshift2007.zip archive."""
    with io.open(ab_filepath, encoding="ascii") as fp_in:
        lines = fp_in.readlines()[3:]  # skip first three comment lines
        prog_marker = "C  program "
        program_start_lines = [
            i for i, line in enumerate(lines) if line.startswith(prog_marker)
        ]
        programs = {}
        for i, marker in enumerate(program_start_lines):
            prog_name = (
                lines[marker].lstrip(prog_marker).lower().strip("\n").strip("\r")
            )
            filename = os.path.join(
                output_dirpath or os.path.dirname(ab_filepath), prog_name
            )
            with io.open(filename, mode="w", encoding="ascii") as fp_out:
                if i == len(program_start_lines) - 1:
                    program_lines = lines[marker - 1 :].copy()
                else:
                    program_lines = lines[marker - 1 : program_start_lines[i + 1] - 1].copy()
                header_lines = [line.lstrip(" \t").lower() for line in program_lines[:100]]
                if not any(line.startswith("program ") for line in header_lines):
                    # find the first non-comment line
                    index = 0
                    for j, line in enumerate(header_lines):
                        if line[0] not in ("c", "!"):
                            index = j
                            break
                    program_lines.insert(index, "      PROGRAM %s\n" % prog_name.upper().strip(".FOR"))
                else:
                    # Rename the existing program line to match the PHSH*.FOR program comment header.
                    # This avoids conflicts between different versions of similarly named programs/subroutines
                    # when comparing documentation.
                    index = next(
                        i for i, line in enumerate(header_lines) if line.startswith("program ")
                    )
                    program_lines[index] = "      PROGRAM %s\n" % prog_name.upper().strip(".FOR")
                # Write a notice to let users know that this is not the original source
                fp_out.writelines(_generate_notice_lines())
                # NOTE: FORD currently does not support tabs in the source code, so replace with spaces
                fp_out.writelines(map(lambda x: x.replace("\t", "      "), program_lines))
            programs[prog_name] = filename
        return programs


def extract_phshift2007(dirpath=DEFAULT_DIRPATH):
    """Extracts the phshift2007 Fortran programs from the phshift2007.zip archive."""
    for filename in ("psprog.ab3", "psprog.ab4"):
        filepath = os.path.join(dirpath, filename)
        print(_extract_fortran_programs_from_ab_file(filepath, dirpath))


def _compile_fortran_program(source_filepath, output_filepath=None):
    # type: (str, str|None) -> str
    """Compiles the Fortran program from the .for file in the phshift2007.zip archive."""
    if not output_filepath:
        output_filepath = os.path.join(
            os.path.dirname(source_filepath),
            os.path.basename(source_filepath).replace(".for", ""),
        )
    args = [
        "gfortran",
        "-pie",
        "-Wall",
        "-static-libgcc",
        "-static-libgfortran",
        "-frecursive",
        "-fcheck=bounds",
        "-std=legacy",
        "-o",
        output_filepath,
        source_filepath,
    ]
    print(" ".join(args))
    with subprocess.Popen(args) as proc:  # nosec
        status = proc.wait()
        if status != 0:
            raise subprocess.CalledProcessError(status, " ".join(map(_as_string, proc.args)))  # type: ignore[arg-type]
    return output_filepath


def _as_string(string_or_bytes):
    # type: (str|bytes) -> str
    return string_or_bytes.decode("utf-8") if isinstance(string_or_bytes, bytes) else string_or_bytes

def compile_phshift2007(dirpath=DEFAULT_DIRPATH):
    """Compiles the phshift2007 Fortran programs."""
    for program in ("phsh0", "phsh1", "phsh2cav", "phsh2wil", "phsh2rel", "phsh3"):
        source_filepath = os.path.join(dirpath, program + ".for")
        print("Compiled %s" % _compile_fortran_program(source_filepath))


def do_action(action, output_dirpath=DEFAULT_DIRPATH, force=False):
    """Perform the specified action."""
    if action in ("download", "all"):
        zipfile_path = download_phshift2007(output_dirpath, force=force)
        print("Downloaded phshift2007 to {}".format(zipfile_path))

    if action in ("unzip", "all"):
        unzip_phshift2007(dirpath=output_dirpath)
        print("Unzipped phshift2007 to {}".format(output_dirpath))

    if action in ("extract", "all"):
        extract_phshift2007(output_dirpath)
        print("Extracted phshift2007 Fortran programs to {}".format(output_dirpath))

    if action in ("compile", "all"):
        compile_phshift2007(output_dirpath)
        print("Compiled phshift2007 Fortran programs to {}".format(output_dirpath))


def get_parser():
    # type: () -> argparse.ArgumentParser
    """Get the parser object for the script."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dirpath",
        type=str,
        default=DEFAULT_DIRPATH,
        help="The directory path to download and extract phshift2007 to.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force the download of phshift2007 even if it already exists.",
    )
    parser.add_argument(
        "--action",
        dest="actions",
        type=str,
        default="all",
        nargs="+",
        choices=["download", "unzip", "extract", "compile", "all"],
        help="The action(s) to perform [default %(default)s].",
    )
    return parser


def main(argv=None):
    """Main entry point for the script."""
    parser = get_parser()
    args, _ = parser.parse_known_args(argv or sys.argv)
    for action in args.actions:
        do_action(action=action, output_dirpath=args.output_dirpath, force=args.force)


if __name__ == "__main__":
    main()
