'''
Created on 10 Jul 2016

@author: Liam Deacon
'''
import os
import sys

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

import ase.data
from ase.spacegroup import Spacegroup

try:
    from pymatgen.matproj.rest import MPRester
except ImportError:
    MPRester = None

crystal_definitions = [
    ('Spacegroup', 1, True, [1, 1, 1], [3.0, 3.0, 3.0, 90.0, 90.0, 90.0], 
        [0, 0, 0, 0, 0, 0], [True, True, True, True, True, True], [['', '', '', '']]), 
    ('fcc', 225, False, [1, 1, 1], [3.0, 3.0, 3.0, 90.0, 90.0, 90.0], 
        [0, 1, 1, 3, 3, 3], [False, False, False, False, False, False], [['', '', '', '']]), 
    ('bcc', 229, False, [1, 1, 1], [3.0, 3.0, 3.0, 90.0, 90.0, 90.0], 
        [0, 1, 1, 3, 3, 3], [False, False, False, False, False, False], [['', '', '', '']]), 
    ('diamond', 227, False, [1, 1, 1], [3.0, 3.0, 3.0, 90.0, 90.0, 90.0], 
        [0, 1, 1, 3, 3, 3], [False, False, False, False, False, False], [['', '', '', '']]), 
    ('hcp', 194, False, [1, 1, 1], [3.0, 3.0, 3.0, 90.0, 90.0, 120.0], 
        [0, 1, 0, 3, 3, 3], [False, False, False, False, False, False], [['', '1./3.', '2./3.', '3./4.']]), 
    ('graphite', 186, False, [1, 1, 1], [3.0, 3.0, 3.0, 90.0, 90.0, 120.0], 
        [0, 1, 0, 3, 3, 3], [False, False, False, False, False, False], [['', '0', '0', '0'], ['', '1./3.', 
         '2./3.', '0']]), 
    ('rocksalt', 225, False, [1, 1, 1], [3.0, 3.0, 3.0, 90.0, 90.0, 90.0], 
        [0, 1, 1, 3, 3, 3], [False, False, False, False, False, False], [['', '0', '0', '0'], ['', '0.5', 
         '0.5', '0.5']]), 
    ('rutile', 136, False, [1, 1, 1], [3.0, 3.0, 3.0, 90.0, 90.0, 90.0], 
        [0, 1, 0, 3, 3, 3], [False, False, False, False, False, False], [['', '0', '0', '0'], ['O', '0.3', 
         '0.3', '0']]),
    ('wurtzite', 186, False, [1, 1, 1], [3.0, 3.0, 3.0, 90.0, 90.0, 120.0], 
        [0, 1, 1, 3, 3, 3], [False, False, False, False, False, False], 
        [['', '0', '0', '0'], ['', '0.5', '0.5', '0.5']])]

try:
    from . import res_rc
    from ..model import Atom, Unitcell
except ValueError:
    try:
        from phaseshifts.gui import res_rc
        from phaseshifts.model import Atom, Unitcell
    except ImportError:
        import res_rc
        sys.path.insert(0, os.path.abspath(os.path.join("..", os.path.dirname(__file__))))
        try:
            from model import Atom, Unitcell
        except ImportError as err:
            from warnings import warn
            sys.stderr.write("ImportWarning: {}\n".format(err.message))
            try:
                from ase import Atom
            except:
                class DummyAtom(object):
                    def __init__(self, **kwargs):
                        self.__dict__.update(kwargs)
                Atom = DummyAtom
            Unitcell = None
    

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
        self.row_data = []
        
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
        from collections import Iterable
        if isinstance(atoms, Iterable):
            self._atoms = atoms
        else:
            self._atoms.append(atoms)
    
    def updateToolTip(self):
        pass
        
    def addAtom(self, atom=None):
        if not atom:
            item, ok = QInputDialog.getText(self, "Add atom", 
                                            "Enter element name, number or phaseshifts tag", 
                                            text='C')
            if ok and item:
                symbol = item.split('_')[0]
                atom = Atom(symbol=symbol, position=[0., 0., 0.], tag=item, charge=0)
                valence = 0.
                radius = None
                try:
                    atom = Atom(**Atom.tag_info(item))
                    symbol = atom.symbol
                    valence = atom.valence
                    radius = atom.radius
                except AttributeError:
                    valence = atom.charge
                finally:
                    row = self.rowCount()
                    self.insertRow(row)
                    self.setVerticalHeaderLabels([i['label'] for i in self.row_data] + [symbol])
                    row_data = OrderedDict((('tag', QTableWidgetItem(atom.tag)),
                                            ('x', QTableWidgetItem(str(atom.position[0]))),
                                            ('y', QTableWidgetItem(str(atom.position[1]))),
                                            ('z', QTableWidgetItem(str(atom.position[2]))),
                                            ('Q', QTableWidgetItem(str(valence))),
                                            ('r', QTableWidgetItem(str(radius))),
                                            ('label', symbol),
                                            ('atom', atom)))
                    self.row_data.append(row_data)
                    
                    for (col, name) in enumerate(self.column_headers):
                        self.setItem(row, col, row_data[name])
                        
                    self.atoms.append(atom)
                    
                #self.ui.table.atom_items
        print('adding: {}'.format(atom))
    
    def deleteAtom(self, atom=None):
        row = self.currentRow()
        if row == -1:
            return  # invalid row
        print('deleting: {}'.format(row))
        self.atoms.pop(row)
        row_data = self.row_data.pop(row)
        self.removeRow(row)
        


class BulkCrystalDialog(QWidget):
    crystal_definitions = crystal_definitions
    modelChanged = Signal(object)
    
    def __init__(self, parent=None, model=None):
        super(self.__class__, self).__init__(parent)
        
        self.model = model
        
        ui_filename = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                   "CreateBulkCrystalDialog.ui"))
        self.ui = uic.loadUi(ui_filename, self)
        
        self.ui.MAPIEdit.setText(os.environ.get('MAPI_KEY', '<Not Set>'))
        
        
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
        
        self.ui.getFromDatabaseButton.clicked.connect(self.getFromDatabase)
        self.ui.pymatgenRadio.setEnabled(MPRester is not None)
        
        self.ui.addButton.clicked.connect(self.table.addAtom)
        self.ui.removeButton.clicked.connect(self.table.deleteAtom)
    
    def updateModel(self):
        
        return self.model
    
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
    
    def getFromDatabase(self, *args):
        if self.ui.aseRadio.isChecked():
            element = self.elements[0][0].get_text()
            z = ase.data.atomic_numbers[self.legal_element]
            ref = ase.data.reference_states[z]
            lattice = ref['symmetry']
            index = 0
            while index < len(crystal_definitions) and crystal_definitions[index][0] != lattice:
                index += 1
            if index == len(crystal_definitions) or not self.legal_element:
                QInputDialog.error(_("Can't find lattice definition!"))
                return False
            self.structinfo.set_active(index)
            self.lattice_lbuts[0].set_value(ref['a'])
            if lattice == 'hcp':
                self.lattice_lbuts[2].set_value(ref['c/a']*ref['a'])
            self.elements[0][0].set_text(element)
            if lattice in ['fcc', 'bcc', 'diamond']:
                self.elements[0][1].set_text('0')
                self.elements[0][2].set_text('0')
                self.elements[0][3].set_text('0')
        elif self.ui.pymatgenRadiu.isChecked():
            print('use pymatgen here')
            compound, ok = QInputDialog.getText(self, "Search Materials Project", 
                                                "Enter element or compound:")
            if not ok:
                return
            search_text = compound.replace(" ", "") + " " if compound else ''
            space_group = SpaceGroup(self.space_group).symbol.replace(" ", "")
            search_text += space_group
            with MPRester("USER_API_KEY") as m:
            
                # Get the formulas and energies of materials with materials_id mp-1234
                # or with formula FeO.
                results = m.query(space_group, ['structure'])
                print(results)


if __name__ == '__main__':
    from qtsix.Qt import QApplication
    import sys
    app = QApplication(sys.argv)
    
    
    
    widget = BulkCrystalDialog()
    widget.show()
    
    
    app.exec_()