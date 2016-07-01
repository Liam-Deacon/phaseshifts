'''
Created on 1 Jul 2016

@author: Liam Deacon

@contact: liam.m.deacon@gmail.com

@copyright: Copyright 2016 Liam Deacon

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
from __future__ import (absolute_import, division, 
                        print_function, with_statement, unicode_literals)
from qtsix import uic
from qtsix.QtWidgets import (QAbstractButton, QDialog)
import os

try:
    import res_rc
except ImportError:
    try:
        import sys
        sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
        import res_rc
    except ImportError:
        import phaseshifts.gui.res_rc as res_rc   


class ModelBuilderDialog(QDialog):
    '''
    Dialog class for updating sequences 
    '''
    def __init__(self, parent=None, model=None):
        super(ModelBuilderDialog, self).__init__(parent)
        
        # dynamically load ui
        path = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                            os.path.basename(__file__).replace('.py', '.ui')))
        self.ui = uic.loadUi(path, self)
        self.initUi()
        
        self.ui.show()
            
    def initUi(self):
        # Setup slots and signals
        pass
