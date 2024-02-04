"""Main module for the phaseshifts.gui package.

This module provides a command line interface for launching the GUI.
It is mostly intended for development and debugging purposes.
"""

# pylint: disable=consider-using-f-string

import argparse
import importlib
import sys

from qtpy import QtWidgets


def get_parser():
    """Get parser object for command line arguments."""
    parser = argparse.ArgumentParser(description="Update dialog for phaseshifts")
    parser.add_argument(
        "-W",
        "--widget",
        dest="widget",
        default="MainWindow",
        help="Widget to load, by default %(default)s",
    )
    return parser

def main(argv=None):
    """Entry point if executing as standalone."""
    if argv is None:
        argv = sys.argv

    parser = get_parser()
    args, _ = parser.parse_known_args(argv)

    app = QtWidgets.QApplication(argv)

    # Load the widget
    module = importlib.import_module("phaseshifts.gui.{}".format(args.widget))
    _window = getattr(module, args.widget)()

    return app.exec_()


# Execute main function if running as standalone module
if __name__ == "__main__":
    main()
