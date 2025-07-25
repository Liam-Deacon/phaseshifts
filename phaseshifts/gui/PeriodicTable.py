#!/usr/bin/env python3
# encoding: utf-8
"""
atorb.py

@author: Liam Deacon

@contact: liam.m.deacon@gmail.com

@copyright: Copyright (C) 2014 Liam Deacon

@license: MIT License (see LICENSE file for details)

@summary: Periodic Table using Qt

"""

# pylint: disable=invalid-name,super-with-arguments,consider-using-f-string

# Import standard libraries
import os
import re
import sys
from collections import OrderedDict

from qtpy import QtGui, QtCore, QtWidgets, uic

from phaseshifts.elements import ELEMENTS

# Load resources & local modules
from .. import elements

ELEMENTS_DICT = OrderedDict([(element.symbol, element.name) for element in ELEMENTS])


# Create a class for our main window
class PeriodicTableDialog(QtWidgets.QFrame):
    """Periodic table dialog class"""

    def __init__(self, parent=None, toggle=True):  # pylint: disable=unused-argument
        super(PeriodicTableDialog, self).__init__(parent)

        # Or more dynamically
        self.ui = uic.loadUi(
            os.path.join(os.path.dirname(__file__), "PeriodicTable.ui"), self
        )
        self.ui.show()

        self.selectedElement = "H"  # default is Hydrogen

        # self.init()
        self.initUi()

    # Overload exit event to write configuration before exiting app
    def closeEvent(self, event):  # pylint: disable=unused-argument
        print(self.selectedElement)
        sys.exit(0)

    # Setup extra UI elements
    def initUi(self):
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)

        # Setup slots

        # set defaults for each element
        for i in range(1, len(elements.ELEMENTS)):
            symbol = ELEMENTS_DICT.keys()[i - 1]
            element = elements.ELEMENTS[symbol]  # type: elements.Element

            config = re.sub(r"([spdf])", r"\1<sup>", element.eleconfig)
            config = config.replace(" ", "</sup>&nbsp;") + "</sup>"
            tooltip = ""
            tooltip = (
                tooltip
                + """
                <html>
                    <span style=" font-weight:600;">{name}</span><br/>
                    Z={protons}<br/>
                    {mass}&nbsp;amu<br/>
                    {config}<br/>
                    T<sub>melt</sub>={tmelt}&nbsp;K<br/>
                    T<sub>boil</sub>={tboil}&nbsp;K<br/>
                    &#961;={density}&nbsp;g/L<br/>
                    &#967;={eleneg}
                </html>""".format(
                    protons=element.protons,
                    name=element.name,
                    tmelt=element.tmelt,
                    tboil=element.tboil,
                    config=config,
                    mass=element.mass,
                    density=element.density,
                    eleneg=element.eleneg,
                )
            )

            # pylint: disable=eval-used
            eval(
                "self.element_%i.clicked.connect(self.buttonClick)" % i,
                {"self", self},
                {},
            )
            eval("""self.element_%i.setToolTip(tooltip)""" % i, {"self", self}, {})
            # pylint: enable=eval-used

    def buttonClick(self):
        self.selectedElement = elements.ELEMENTS[
            ELEMENTS_DICT.keys().index(self.sender().text()) + 1
        ]
        self.ui.labelMass.setText(
            """<html><sup>%.1f</sup></html>""" % float(self.selectedElement.mass)
        )
        self.ui.labelZ.setText(
            """<html><sup>%s</sup></html>""" % self.selectedElement.protons
        )
        self.ui.labelElement.setText(
            """<html><head/><body><p><span style="
                                        font-size:12pt;">%s</span></p></body>
                                        </html>"""
            % self.selectedElement.symbol
        )
