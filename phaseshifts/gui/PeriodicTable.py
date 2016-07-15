#!/usr/bin/env python
#encoding: utf-8
#
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.m.deacon@gmail.com                                           #
#                                                                            #
# Copyright: Copyright (C) 2014-2016 Liam Deacon                             #
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
PeridoicTable.py

@author: Liam Deacon

@contact: liam.m.deacon@gmail.com

@copyright: Copyright (C) 2014-2016 Liam Deacon

@license: MIT License (see LICENSE file for details)

@summary: Periodic Table using Qt

"""
from __future__ import absolute_import, unicode_literals

# Import standard libraries
import logging
import ntpath
import re
import os
import platform
import sys

try:
    from . MainWindow import qt_api
    from . import res_rc
except ValueError:
    from phaseshifts.gui.MainWindow import qt_api
    import phaseshifts.gui.res_rc as res_rc        
    
# Import Qt modules
from qtsix import uic, QtCore
from qtsix.Qt import Qt
from qtsix.QtGui import QIcon
from qtsix.QtCore import Signal, Slot
from qtsix.QtWidgets import (QApplication, QFrame, QDialog,
                             QDialogButtonBox, QVBoxLayout)


try:
    import elements
except ImportError:
    try:
        from .. import elements
    except (ValueError, ImportError):
        try:
            from phaseshifts import elements
        except ImportError:
            sys.path.insert(0, os.path.abspath(os.path.dirname(
                            os.path.dirname(os.path.dirname(__file__)))))
            import elements


import re
from collections import OrderedDict

__APP_NAME__ = 'PeriodicTable'

elements_dict = OrderedDict([
('H', 'Hydrogen'), 
('He', 'Helium'), 
('Li', 'Lithium'), 
('Be', 'Beryllium'), 
('B', 'Boron'), 
('C', 'Carbon'), 
('N', 'Nitrogen'), 
('O', 'Oxygen'), 
('F', 'Fluorine'), 
('Ne', 'Neon'), 
('Na', 'Sodium'), 
('Mg', 'Magnesium'), 
('Al', 'Aluminium'), 
('Si', 'Silicon'), 
('P', 'Phosphorus'), 
('S', 'Sulfur'), 
('Cl', 'Chlorine'), 
('Ar', 'Argon'), 
('K', 'Potassium'), 
('Ca', 'Calcium'), 
('Sc', 'Scandium'), 
('Ti', 'Titanium'), 
('V', 'Vanadium'), 
('Cr', 'Chromium'), 
('Mn', 'Manganese'), 
('Fe', 'Iron'), 
('Co', 'Cobalt'), 
('Ni', 'Nickel'), 
('Cu', 'Copper'), 
('Zn', 'Zinc'), 
('Ga', 'Gallium'), 
('Ge', 'Germanium'), 
('As', 'Arsenic'), 
('Se', 'Selenium'), 
('Br', 'Bromine'), 
('Kr', 'Krypton'), 
('Rb', 'Rubidium'), 
('Sr', 'Strontium'), 
('Y', 'Yttrium'), 
('Zr', 'Zirconium'), 
('Nb', 'Niobium'), 
('Mo', 'Molybdenum'), 
('Tc', 'Technetium'), 
('Ru', 'Ruthenium'), 
('Rh', 'Rhodium'), 
('Pd', 'Palladium'), 
('Ag', 'Silver'), 
('Cd', 'Cadmium'), 
('In', 'Indium'), 
('Sn', 'Tin'), 
('Sb', 'Antimony'), 
('Te', 'Tellurium'), 
('I', 'Iodine'), 
('Xe', 'Xenon'), 
('Cs', 'Cesium'), 
('Ba', 'Barium'), 
('La', 'Lanthanum'), 
('Ce', 'Cerium'), 
('Pr', 'Praseodymium'), 
('Nd', 'Neodymium'), 
('Pm', 'Promethium'), 
('Sm', 'Samarium'), 
('Eu', 'Europium'), 
('Gd', 'Gadolinium'), 
('Tb', 'Terbium'), 
('Dy', 'Dysprosium'), 
('Ho', 'Holmium'), 
('Er', 'Erbium'), 
('Tm', 'Thulium'), 
('Yb', 'Ytterbium'), 
('Lu', 'Lutetium'), 
('Hf', 'Hafnium'), 
('Ta', 'Tantalum'), 
('W', 'Tungsten'), 
('Re', 'Rhenium'), 
('Os', 'Osmium'), 
('Ir', 'Iridium'), 
('Pt', 'Platinum'), 
('Au', 'Gold'), 
('Hg', 'Mercury'), 
('Tl', 'Thallium'), 
('Pb', 'Lead'), 
('Bi', 'Bismuth'), 
('Po', 'Polonium'), 
('At', 'Astatine'), 
('Rn', 'Radon'), 
('Fr', 'Francium'), 
('Ra', 'Radium'), 
('Ac', 'Actinium'), 
('Th', 'Thorium'), 
('Pa', 'Protactinium'), 
('U', 'Uranium'), 
('Np', 'Neptunium'), 
('Pu', 'Plutonium'), 
('Am', 'Americium'), 
('Cm', 'Curium'), 
('Bk', 'Berkelium'), 
('Cf', 'Californium'), 
('Es', 'Einsteinium'), 
('Fm', 'Fermium'), 
('Md', 'Mendelevium'), 
('No', 'Nobelium'), 
('Lr', 'Lawrencium'), 
('Rf', 'Rutherfordium'), 
('Db', 'Dubnium'), 
('Sg', 'Seaborgium'), 
('Bh', 'Bohrium'), 
('Hs', 'Hassium'), 
('Mt', 'Meitnerium'), 
('Ds', 'Darmstadtium'), 
('Rg', 'Roentgenium'), 
('Cn', 'Copernicium'), 
('Uut', 'Ununtrium'), 
('Fl', 'Flerovium'), 
('Uup', 'Ununpentium'), 
('Lv', 'Livermorium'), 
('Uus', 'Ununseptium'), 
('Uuo', 'Ununoctium'), 
])


# Create a class for our main window
class PeriodicTable(QFrame):
    """Periodic table dialog class"""
    
    selectedElementChanged = Signal(object)
    
    def __init__(self, parent=None, toggle=True, element='H'):
        super(self.__class__, self).__init__(parent)
 
        # Or more dynamically
        self.ui = uic.loadUi(os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                                          "PeriodicTable.ui")), self)
        self.ui.show()
        
        self.selectedElement = element  # default is Hydrogen
        
        #self.init()
        self.initUi()
    
    # Overload exit event to write configuration before exiting app
    def closeEvent(self, evnt):
        pass # print(self.selectedElement)
        

    # Setup extra UI elements
    def initUi(self):
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        
        # Setup slots
        
        # set defaults for each element
        self.elements = []
        for i in range(1, len(elements.ELEMENTS)):
            symbol = elements_dict.keys()[i - 1]
            element = elements.ELEMENTS[symbol]
            
            config = re.sub(r'([spdf])', r'\1<sup>', element.eleconfig)
            config = config.replace(' ', '</sup>&nbsp;') + '</sup>'
            element.description
            #eleaffin=0.0,
            #covrad=0.0, atmrad=0.0, vdwrad=0.0,
            tooltip = """
                <html>
                    <div style="width: 300px">
                    <span style=" font-weight:600;font-size:18px;">{name}</span><br/>
                    <span style="font-size:14px;">
                    Z={protons}<br/>
                    {mass}&nbsp;amu<br/>
                    {config}<br/>
                    T<sub>melt</sub>={tmelt}&nbsp;K<br/>
                    T<sub>boil</sub>={tboil}&nbsp;K<br/>
                    &#961;={density}&nbsp;g/L<br/>
                    &#967;={eleneg}<br/>
                    r<sub>atomic</sub>={atmrad}&nbsp;&#8491;<br/>
                    r<sub>covalent</sub>={covrad}&nbsp;&#8491;<br/>
                    r<sub>van-der-waal</sub>={vdwrad}&nbsp;&#8491;<br/>
                    <br/>
                    </span>
                    </div>                   
                </html>""".format(protons=element.protons, name=element.name,
                                  tmelt=element.tmelt, tboil=element.tboil,
                                  config=re.sub("([0-9]+)([spdf])([0-9]+)", 
                                                "\\1\\2<sup>\\3</sup>", config), 
                                  mass=element.mass, 
                                  density=element.density, 
                                  eleneg=element.eleneg,
                                  atmrad=element.atmrad,
                                  covrad=element.covrad,
                                  vdwrad=element.vdwrad)
            
            eval('self.element_%i.clicked.connect(self.buttonClick)' % i)
            eval('''self.element_%i.setToolTip(tooltip)''' % i)
            self.elements.append(element)
            
    def buttonClick(self):
        i = elements_dict.keys().index(self.sender().text()) + 1
        if elements.ELEMENTS[i] == self.selectedElement:
            return  # do not update
        self.selectedElement = elements.ELEMENTS[i]
        self.ui.labelMass.setText('''<html><sup>%.1f</sup></html>''' 
                                  % float(self.selectedElement.mass))
        self.ui.labelZ.setText('''<html><sup>%s</sup></html>''' 
                                  % self.selectedElement.protons)
        self.ui.labelElement.setText('''<html><head/><body><p><span style=" 
                                        font-size:12pt;">%s</span></p></body>
                                        </html>''' 
                                        % self.selectedElement.symbol)
        self.updatedSelectedElement(i)
        
    @Slot(int)
    def updatedSelectedElement(self, element):
        self.selectedElement = element
        self.selectedElementChanged.emit(element)
        
    @property
    def selectedElement(self):
        return self._element
    
    @selectedElement.setter
    def selectedElement(self, element):
        try:
            self._element = elements.ELEMENTS[element]
        except KeyError:
            pass


class PeriodicTableDialog(QDialog):
    """ Dialog class for PeriodicTable frame """
    selectedElementChanged = Signal(object)
    
    def __init__(self, parent=None, element='H'):
        super(self.__class__, self).__init__(parent)
        
        table = PeriodicTable(None, element)

        layout = QVBoxLayout()
        
        layout.addWidget(table)
        self.table = table
        
        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | 
                                   QDialogButtonBox.Cancel,
                                   Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        
        layout.addWidget(buttons)
        self.setLayout(layout)
        
        self.table.selectedElementChanged.connect(self.updateSelectedElement)
    
    @Slot(int)
    def updateSelectedElement(self, element):
        self.selectedElementChanged.emit(element)
        
    @property
    def selectedElement(self):
        return self.table.selectedElement

def main():
    # Again, this is boilerplate, it's going to be the same on
    # almost every app you write
    app = QApplication(sys.argv)
    icon = QIcon(os.path.join(os.path.dirname(__file__), 
                                    'res', 
                                    'periodictable_32x32.png'))
    app.setWindowIcon(icon)
    window = PeriodicTableDialog()
    
    # Platform specific setup
    if platform.system() is 'Windows':
        from ctypes import windll
        # Tell Windows Python is merely hosting the application (taskbar icon fix)
        windll.shell32.SetCurrentProcessExplicitAppUserModelID(__APP_NAME__)
    
    window.show()
    
    # It's exec_ because exec is a reserved word in Python
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
