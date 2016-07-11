# encoding: utf-8
"""setupwindow.py - Window base class for setup modules.
"""

from qtsix.QWtWidgets import (QLabel, QWidget)
from gettext import gettext as _
from ase.gui.widgets import pack
import ase

def pack(vbox, widgets, end=False, bottom=False, expand=False, padding=0):
    if not isinstance(widgets, list):
        widgets.show()
        if bottom:
            vbox.pack_end(widgets, expand, expand, padding)
        else:
            vbox.pack_start(widgets, expand, expand, padding)
        return widgets
    hbox = gtk.HBox(0, 0)
    hbox.show()
    if bottom:
        vbox.pack_end(hbox, expand, expand, padding)
    else:
        vbox.pack_start(hbox, expand, expand, padding)
    for widget in widgets:
        if type(widget) is gtk.Entry: # isinstance does not work here
            widget.set_size_request(widget.get_max_length() * 9, 24)
        widget.show()
        if end and widget is widgets[-1]:
            hbox.pack_end(widget, expand, expand, padding)
        else:
            hbox.pack_start(widget, expand, expand, padding)
    
    return widgets


class SetupWindow(gtk.Window):
    "Base class for ase.gui setup windows."
    # __init__ inherited from gtk.Window

    def packtext(self, vbox, text, label=None):
        "Pack an text frame into the window."
        pack(vbox, QLabel(""))
        txtframe = gtk.Frame(label)
        txtlbl = QLabel(text)
        txtframe.add(txtlbl)
        txtlbl.show()
        pack(vbox, txtframe)
        pack(vbox, QLabel(""))

    def update_element(self, *args):
        "Called when a new element may have been entered."
        # Assumes the element widget is self.element and that a label
        # to keep updated is self.elementinfo.  The chemical symbol is
        # placed in self.legalelement - or None if the element is
        # invalid.
        elem = self.element.get_text()
        if not elem:
            self.invalid_element(_("  No element specified!"))
            return False
        try:
            z = int(elem)
        except ValueError:
            # Probably a symbol
            try:
                z = ase.data.atomic_numbers[elem]
            except KeyError:
                self.invalid_element()
                return False
        try:
            symb = ase.data.chemical_symbols[z]
        except KeyError:
            self.invalid_element()
            return False
        name = ase.data.atomic_names[z]
        ref = ase.data.reference_states[z]
        if ref is None:
            struct = _("No crystal structure data")
        else:
            struct = ref['symmetry']
            if struct == 'fcc' or struct == 'bcc':
                struct = "%s (a=%.3f �)" % (struct, ref['a'])
        
        txt = "  %s: %s, Z=%i, %s" % (name, symb, z, struct)
        self.elementinfo.set_text(txt)
        self.legal_element = symb
        return True
        
    def invalid_element(self, txt=_("  ERROR: Invalid element!")):
        self.legal_element = False
        self.elementinfo.set_text(txt)
        
if __name__ == "__main__":
    from qtsix.Qt import QApplication
    import sys
    app = QApplication(sys.argv)
    
    window = SetupWindow()
    window.show()
    
    sys.exit(app.exec_())