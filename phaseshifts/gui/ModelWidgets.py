'''
Created on 10 Jul 2016

@author: Liam Deacon
'''
import os
import sys
from __builtin__ import str


try:
    from cyordereddict import OrderedDict
except ImportError:
    from collections import OrderedDict

from qtsix.QtWidgets import (QLabel, QWidget, QDoubleSpinBox, 
                             QHBoxLayout, QVBoxLayout, QWidgetItem,
                             QTableWidget, QTableWidgetItem,
                             QDialogButtonBox, QInputDialog, 
                             QMessageBox)
from qtsix.QtCore import Slot, Signal
from qtsix.QtGui import QIcon
from qtsix import uic

import ase.data
from ase.spacegroup import Spacegroup, crystal
from ase.spacegroup.spacegroup import SpacegroupNotFoundError

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
                atom = Atom(symbol, position=[0., 0., 0.], tag=item, charge=0)
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
    spacegroupChanged = Signal(object)
    latticeChanged = Signal(object)
    
    def __init__(self, parent=None, model=None):
        super(self.__class__, self).__init__(parent)
        
        self.model = model
        self.spacegroup = Spacegroup(1)
        
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
        self.spacegroupChanged.connect(self.updateSpaceGroup)
        self.ui.spaceGroupEdit.editingFinished.connect(
                lambda: self.spacegroupChanged.emit(self.ui.spaceGroupEdit.text()))
        
        self.ui.getFromDatabaseButton.clicked.connect(self.getFromDatabase)
        self.ui.pymatgenRadio.setEnabled(MPRester is not None)
        
        self.ui.addButton.clicked.connect(self.table.addAtom)
        self.ui.removeButton.clicked.connect(self.table.deleteAtom)
        
        self.structinfo = self.ui.spaceGroupCombo
        self.structinfo.clear() 
        self.structinfo.addItems([c[0] for c in self.crystal_definitions] + ["custom"])
        self.structures = dict((c[0], c) for c in self.crystal_definitions)
        self.ui.spaceGroupCombo.currentIndexChanged.connect(self.updateSpaceGroup)
        
        self.ui.spaceGroupCheck.stateChanged.connect(
                lambda i: self.updateSpaceGroup(bool(i)))

    
    @Slot()
    @Slot(bool)
    @Slot(int, name='updateSpaceGroupNumber')
    @Slot(str, name='updateSpaceGroupSymbol')
    def updateSpaceGroup(self, group=None):
        if not group or isinstance(group, bool):
            group = self.spacegroup.no
        if isinstance(group, int):
            try:
                group = self.crystal_definitions[group][1]
            except IndexError:
                group = self.ui.spaceGroupEdit.text()
        try:
            self.spacegroup = Spacegroup(str(group))
        except SpacegroupNotFoundError:
            try:
                self.spacegroup = Spacegroup(int(group))
            except (ValueError, SpacegroupNotFoundError):
                QMessageBox.critical(self, "Invalid Space Group",
                                     "Space group '{}' is not valid".format(group))
                self.ui.spaceGroupEdit.setText(self.spacegroup.symbol)
        finally:
            self.ui.spaceGroupEdit.setText(str(self.spacegroup.symbol if 
                self.ui.spaceGroupCheck.isChecked() else self.spacegroup.no))
            number = self.spacegroup.no
            spacegroups = [c[1] for c in self.crystal_definitions]
            if number in spacegroups:
                self.ui.spaceGroupCombo.setCurrentIndex(spacegroups.index(number))
            else:
                self.ui.spaceGroupCombo.setCurrentIndex(self.ui.spaceGroupCombo.count()-1)
            
            self.updateLattice(self.spacegroup)
            
    @Slot()
    @Slot(Spacegroup)
    def updateLattice(self, lattice=None):
        if isinstance(lattice, Spacegroup):
            index = [c[1] for c in self.crystal_definitions].index(lattice.no)
            cell_pars = self.crystal_definitions[index][4]
            self.ui.aSpinBox.setValue(cell_pars[0])
            self.ui.bSpinBox.setValue(cell_pars[1])
            self.ui.cSpinBox.setValue(cell_pars[2])
            self.ui.alphaSpinBox.setValue(cell_pars[3])
            self.ui.betaSpinBox.setValue(cell_pars[4])
            self.ui.gammaSpinBox.setValue(cell_pars[5])
        else:
            print("To implement")
        
                
    def setLatticeType(self, *args):
        """ set defaults from original """
        self.clearing_in_process = True
        self.clear_lattice()
        lattice = crystal_definitions[self.structinfo.get_active()]
        self.spacegroup.set_text(str(lattice[1]))
        self.spacegroup.set_sensitive(lattice[2])
        for s, i in zip(self.size,lattice[3]):
            s.set_value(i)
        self.lattice_lbuts[0].set_value(lattice[4][0])
        self.lattice_lbuts[1].set_value(lattice[4][1])
        self.lattice_lbuts[2].set_value(lattice[4][2])
        self.lattice_abuts[0].set_value(lattice[4][3])
        self.lattice_abuts[1].set_value(lattice[4][4])
        self.lattice_abuts[2].set_value(lattice[4][5])
        self.lattice_lequals[0].set_active(lattice[5][0])
        self.lattice_lequals[1].set_active(lattice[5][1])
        self.lattice_lequals[2].set_active(lattice[5][2])
        self.lattice_aequals[0].set_active(lattice[5][3])
        self.lattice_aequals[1].set_active(lattice[5][4])
        self.lattice_aequals[2].set_active(lattice[5][5])
        self.lattice_lequals[0].set_sensitive(lattice[6][0])
        self.lattice_lequals[1].set_sensitive(lattice[6][1])
        self.lattice_lequals[2].set_sensitive(lattice[6][2])
        self.lattice_aequals[0].set_sensitive(lattice[6][3])
        self.lattice_aequals[1].set_sensitive(lattice[6][4])
        self.lattice_aequals[2].set_sensitive(lattice[6][5])
        for n, at in enumerate(lattice[7]):
            l = 0
            if n > 0:
                l = len(self.elements)
                self.add_basis_atom()
            for i, s in enumerate(at):
                self.elements[l][i].set_text(s)
        self.clearing_in_process = False
        self.update()
    
    def update(self):
        self.modelChanged.emit(self.model)
        return self.model
    
    def ok(self):
        if self.apply():
            self.close()
    
    def cancel(self):
        self.close()
    
    def reset(self):
        self.ui.a = 1.
        self.ui.table.clearContents()
    
    def apply(self):
        ''' Apply value '''
        self.update()
        if self.atoms:
            return True
        else:
            QMessageBox.critical(self, 'No valid atoms.',
                                 'You have not (yet) specified a '
                                 'consistent set of parameters.')
            return False

    @property
    def atoms(self):
        return self.ui.table.atoms
        
    
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