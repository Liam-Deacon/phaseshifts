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

from PyQt4 import QtGui, uic
import res_rc


class UpdateDialog(QtGui.QDialog):
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
        self.package = package
        self.parent = parent
        self.server = server
        self.proxy = proxy

        # dynamically load ui
        self.ui = uic.loadUi("gui/UpdateDialog.ui", self)
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
            QtGui.QMessageBox.information(self, "Error", "No package name!")
            self.ui.closeEvent()
        import pip

        try:
            dists = [
                dist
                for dist in pip.get_installed_distributions()
                if dist.key == self.package
            ]

        except Exception:
            QtGui.QMessageBox.information(
                self, "Update %s" % self.package, "Failed to update!"
            )

        from phaseshifts.contrib.xmlrpc_urllib2_transport import ProxyTransport

        try:
            pypi = xmlrpclib.ServerProxy(server, transport=ProxyTransport())
            dist = dists[0]
            available = pypi.package_releases(dist.project_name)
            if not available:
                # Try to capitalise package name
                available = pypi.package_releases(dist.project_name.capitalize())

            if not available:
                msg = "no releases at pypi"
            elif available[0] != dist.version:
                msg = "{} available".format(available[0])
            else:
                msg = "up to date"
            pkg_info = "{dist.project_name} {dist.version}".format(dist=dist)
            print("{pkg_info:40} {msg}".format(pkg_info=pkg_info, msg=msg))

        except Exception:
            self.logger.error("Unable to contact PyPI server for updates")
            QtGui.QMessageBox.information(
                self,
                "Update %s" % self.package,
                "Unable to contact PyPI server for updates",
            )

            status = "No update available"
            QtGui.QMessageBox.information(
                self, "Update {package}".format(package=self.package), status
            )
