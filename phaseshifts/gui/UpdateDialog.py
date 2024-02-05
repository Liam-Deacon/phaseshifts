"""
Created on 31 Jan 2014

@author: Liam Deacon

@contact: liam.deacon@diamond.ac.uk

@copyright: Copyright 2014 Liam Deacon

@license: MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to
do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

# pylint: disable=invalid-name,super-with-arguments,consider-using-f-string

import os
import sys
import subprocess  # nosec

from qtpy import QtWidgets, uic


class UpdateDialog(QtWidgets.QDialog):
    """
    Dialog class for updating sequences
    """

    def __init__(
        self,
        parent=None,
        package=None,
        server="http://pypi.python.org/pypi",
        proxy=None,
    ):
        super(UpdateDialog, self).__init__(parent)

        # set dictionary
        self.package = package or "phaseshifts"
        self.parent = parent
        self.server = server
        self.proxy = proxy

        # dynamically load ui
        self.ui = uic.loadUi(
            os.path.join(os.path.dirname(__file__), "UpdateDialog.ui"), self
        )
        self.initUi()

        self.ui.show()

    def initUi(self):
        self.progressBar.hide()
        self.textEdit.hide()

        # Setup slots and signals
        self.ui.checkBoxDetails.stateChanged.connect(self.updateUi)
        self.ui.pushButtonAbort.clicked.connect(self.close)
        self.ui.checkBoxDetails.setChecked(False)
        self.minWidth = self.ui.sizeHint().width()
        self.minHeight = self.ui.minimumSizeHint().height()
        self.ui.resize(self.ui.sizeHint().width(), self.ui.minimumSizeHint().height())

    def updateUi(self):
        if self.sender().isChecked():
            self.textEdit.show()
        else:
            self.textEdit.hide()
            self.ui.resize(self.minWidth, self.minHeight)

    def doUpdate(self):
        if self.package is None:
            QtWidgets.QMessageBox.information(self, "Error", "No package name!")
            self.ui.closeEvent()

        try:
            subprocess.check_call(  # nosec: B603
                [
                    sys.executable,
                    "-m",
                    "pip",
                    "install",
                    "--upgrade",
                    self.package,
                ]
            )

        except subprocess.CalledProcessError as err:
            QtWidgets.QMessageBox.information(
                self,
                "Update %s" % self.package,
                "Unable to update package %s due to: %s" % (self.package, str(err)),
            )
