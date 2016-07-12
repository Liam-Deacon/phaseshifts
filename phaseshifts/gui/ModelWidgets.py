'''
Created on 10 Jul 2016

@author: Liam Deacon
'''
import os

try:
    from cyordereddict import OrderedDict
except ImportError:
    from collections import OrderedDict

from qtsix.QtWidgets import (QLabel, QWidget, QDoubleSpinBox, 
                             QHBoxLayout, QVBoxLayout, QWidgetItem,
                             QTableWidget, QTableWidgetItem,
                             QDialogButtonBox, QInputDialog)
from qtsix.QtCore import Slot, Signal
from qtsix.QtGui import QIcon
from qtsix import uic

try:
    from . import res_rc
    from ..model import Atom, Unitcell
except ValueError:
    try:
        from phaseshifts.gui import res_rc
        from phaseshifts.model import Atom, Unitcell
    except ImportError:
        import res_rc
        Atom = Unitcell = None
    

class UnitCellWidget(QWidget):
    UNITCELL_CLASS = Unitcell
    unitCellChanged = Signal()
    
    def __init__(self, parent=None, unitcell=None):
        super(self.__class__, self).__init__(parent)
        
        self.vbox = QVBoxLayout(self)
        
        self.unitcell = unitcell or Unitcell()
        
        for vec in ('a', 'b', 'c'):
            hbox = QHBoxLayout()
            exec('self.label_{0} = QLabel("{0}"); '
                 'hbox.addWidget(self.label_{0})'.format(vec))
            for i in range(3):
                exec('self.{v}{i} = QDoubleSpinBox(self);'
                     'hbox.addWidget(self.{v}{i});'
                     'self.{v}{i}.valueChanged.connect(self.updateUnitCell)'
                     ''.format(v=vec, i=i))
            self.vbox.addLayout(hbox)
        
    
    @Slot()
    def updateUnitCell(self, unitcell=None):
        unitcell = unitcell if unitcell else self.unitcell
        vec = [[self.a1.value(), self.a2.value(), self.a3.value()],
               [self.b1.value(), self.b2.value(), self.b3.value()],
               [self.c1.value(), self.c2.value(), self.c3.value()],]
        self.unitCellChanged.emit()
        try:
            from pymatgen.core import Lattice
            import numpy as np
            return Lattice(np.array(vec))
        except ImportError:
            return unitcell.__class__(vec)


class AtomsTable(QTableWidget):
    
    column_headers = OrderedDict((('tag', {'tooltip': 'Phaseshifts tag for species'}), 
                                  ('x', {'tooltip': 'x-position in Angstroms'}), 
                                  ('y', {'tooltip': 'y-position in Angstroms'}), 
                                  ('z', {'tooltip': 'z-position in Angstroms'}), 
                                  ('Q', {'tooltip': 'Valency of specie e.g. -2'}),
                                  ('r', {'tooltip': 'Muffin-Tin radius'})))
    
    def __init__(self, parent=None, atoms=[]):
        super(self.__class__, self).__init__(parent)
        self.setObjectName("AtomsTable")
        
        self.atoms = atoms
        
        self.setColumnCount(len(self.column_headers))
        for (i, header) in enumerate(self.column_headers):
            #self.insertColumn(i)
            item = QTableWidgetItem(header)
            item.setToolTip(self.column_headers[header]['tooltip'])
            icon = self.column_headers[header].get('icon', None)
            if icon:
                item.setIcon(QIcon(icon))
            self.setVerticalHeaderItem(i, item)
        self.setHorizontalHeaderLabels(list(self.column_headers.keys()))
        
        self.atom_items = []
        
        self.setToolTip('List of atoms in model')
    
    @property
    def atom_items(self):
        return self._atom_items
    
    @atom_items.setter
    def atom_items(self, items):
        self._atom_items = items
    
    @property
    def atoms(self):
        return self._atoms
    
    @atoms.setter
    def atoms(self, atoms):
        self.clear()
        self._atoms = atoms
    
    def updateToolTip(self):
        pass
        
    def addAtom(self, atom=None):
        if not atom:
            item, ok = QInputDialog.getText(self, "Add atom", 
                                            "Enter element name, number or phaseshifts tag", 
                                            text='C')
            if ok and item:
                #atom = Atom(**Atom.tag_info(item))
                row = self.rowCount()
                self.insertRow(row)
                self.setVerticalHeaderLabels([item])
                #self.ui.table.atom_items
        print('adding: {}'.format(atom))
    
    def deleteAtom(self, atom=None):
        print('deleting: {}'.format(atom))


class BulkCrystalDialog(QWidget):
    modelChanged = Signal(object)
    
    def __init__(self, parent=None, model=None):
        super(self.__class__, self).__init__(parent)
        
        self.model = model
        
        ui_filename = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                   "CreateBulkCrystalDialog.ui"))
        self.ui = uic.loadUi(ui_filename, self)
        
        
        self._init_ui()
        
    def _init_ui(self):
        # create atoms table
        table = AtomsTable()
        self.ui.table.hide()  # 'delete' vanilla QTableWidget
        self.ui.basisVbox.addWidget(table)
        self.table = table
        
        # connect actions
        self.ui.buttonBox.button(QDialogButtonBox.Apply).clicked.connect(self.apply)
        self.ui.buttonBox.button(QDialogButtonBox.Cancel).clicked.connect(self.cancel)
        self.ui.buttonBox.button(QDialogButtonBox.Reset).clicked.connect(self.reset)
        self.ui.buttonBox.button(QDialogButtonBox.Ok).clicked.connect(self.ok)
        
        self.modelChanged.connect(lambda x: sys.stdout.write('{}'.format(x)))
        
        self.ui.addButton.clicked.connect(self.table.addAtom)
        self.ui.removeButton.clicked.connect(self.table.deleteAtom)
        
    def _doButtonClick(self, i):
        if i == QDialogButtonBox.Ok:
            print('Ok')
        print(i)
        
    def ok(self):
        self.apply()
        self.close()
    
    def cancel(self):
        self.close()
    
    def reset(self):
        self.ui.a = 1.
        self.ui.table.clearContents()
    
    def apply(self):
        ''' Apply value '''
        self.modelChanged.emit(self.model)
        

if __name__ == '__main__':
    from qtsix.Qt import QApplication
    import sys
    app = QApplication(sys.argv)
    
    widget = BulkCrystalDialog()
    widget.show()
    
    
    app.exec_()