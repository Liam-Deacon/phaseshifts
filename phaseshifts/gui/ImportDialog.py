'''
Created on 10 Feb 2014

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
'''
from PyQt4 import QtGui, uic
import res_rc


class ImportDialog(QtGui.QDialog):
    '''
    Dialog class for updating sequences 
    '''
    def __init__(self, parent=None, model=None):
        super(ImportDialog, self).__init__(parent)
        
        # set dictionary
        self.action = None
        
        # dynamically load ui
        self.ui = uic.loadUi("gui/ImportDialog.ui", self)
        self.initUi()
        
        if isinstance(model, str):
            model = model.lower()
            if model == 'slab':
                self.ui.radioSlab.setChecked(True)
        
        self.ui.show()
            
    def initUi(self):
        # Setup slots and signals
        self.ui.buttonBox.clicked[
                        QtGui.QAbstractButton].connect(self.buttonPress)
    
    def buttonPress(self, button):
        '''Deal with user interaction of button group'''
        action = str(button.text()).lower()
        if action == 'cancel':
            # do not apply settings & close dialog
            self.action = action
            self.ui.close()
        
        elif action == 'ok':
            self.action = action 
            self.ui.close()
