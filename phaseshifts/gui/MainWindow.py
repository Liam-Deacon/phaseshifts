'''
Created on 30 Jan 2014

@author: Liam Deacon

@contact: liam.deacon@diamond.ac.uk

@copyright: 2014 Liam Deacon

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

# Python version compatibility
from __future__ import print_function, with_statement

# Import standard library modules
import logging
import ntpath
import os
import platform
import sys
from collections import OrderedDict

# Import Qt modules
import PyQt4
from PyQt4 import QtCore, QtGui, uic
import res_rc  # note this requires compiled resource file res_rc.py
__QT_TYPE__ = 'PyQt4' 

# other modules
from settings import Settings
from ImportDialog import ImportDialog
from phaseshifts.model import MTZ_model, Unitcell, Atom

# Define globals
__APP_AUTHOR__ = 'Liam Deacon'
__APP_COPYRIGHT__ = '\xa9' + '2013 {0}'.format(__APP_AUTHOR__)
__APP_DESCRIPTION__ = ('A simple Python-based program \n '
                        'for generation of phase shifts')
__APP_DISTRIBUTION__ = 'phaseshifts'
__APP_EMAIL__ = 'liam.m.deacon@diamond.ac.uk'
__APP_LICENSE__ = 'MIT License'
__APP_NAME__ = 'Phase Shifts'
__APP_VERSION__ = '0.1-alpha'
__PYTHON__ = "{0}.{1}.{2}".format(sys.version_info.major,
                                         sys.version_info.minor,
                                         sys.version_info.micro, 
                                         sys.version_info.releaselevel)
__UPDATE_URL__ = ""


# Platform specific setup
if platform.system() is 'Windows':
    from ctypes import windll
    # Tell Windows Python is merely hosting the application (taskbar icon fix)
    windll.shell32.SetCurrentProcessExplicitAppUserModelID(__APP_NAME__)


#==============================================================================
# BEGIN GUI WIDGET PROGRAMMING
#==============================================================================

class MainWindow(QtGui.QMainWindow):
    '''Class for main application widget'''
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        
        # dynamically load ui
        uiFile = "gui/MainWindow.ui"  # change to desired relative ui file path
        self.ui = uic.loadUi(uiFile, self) 
        self.ui.show()
        
        self.init()
        self.initUi()

    # initialise class
    def init(self):
        '''Class to initialise logging and non-gui objects''' 

        ######################################
        # APP LOGGING
        ######################################
        
        # create logger with 'spam_application'
        self.logger = logging.getLogger(__APP_NAME__)
        self.logger.setLevel(logging.DEBUG)
        
        # create file handler which logs all messages
        fh = logging.FileHandler(os.path.join(os.environ['TEMP'], __APP_NAME__ 
                + str('.log')))  # temp directory is emptied on system reboot
        formatter = logging.Formatter(
                        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        fh.setLevel(logging.INFO)  # change to taste
        
        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.WARNING)
        ch.setFormatter(formatter)
        
        # add the handlers to the logger
        self.logger.addHandler(ch)
        self.logger.addHandler(fh)
        
        ######################################
        # SETTINGS
        ######################################
        self.lastpath = ''  # last used file/directory path
        self.settings = Settings()
        
        # tree dictionaries
        #self.treeRootDict = self.ui.getChildItems(self.ui.treeWidgetBulk)
        bulkTree = self.ui.treeWidgetBulk.invisibleRootItem()
        slabTree = self.ui.treeWidgetSlab.invisibleRootItem()
        #bulk = OrderedDict([str(bulkTree.child(i).text(0)) for i in range(bulkTree.childCount())])
        #slab = OrderedDict([str(slabTree.child(i).text(0)) for i in range(slabTree.childCount())])
        #print(bulk)
        #od = OrderedDict([('bulk', bulk),
        #                  ('slab', slab)])
        
        #print(od)
        
    # initialise UI
    def initUi(self):
        '''Class to initialise the Qt Widget and setup slots and signals'''

        # Setup slots
        
        # actions
        self.ui.actionAbout.triggered.connect(self.about)
        self.ui.actionAboutQt.triggered.connect(self.aboutQt)
        self.ui.actionContact.triggered.connect(self.contactDeveloper)
        self.ui.actionExport.triggered.connect(self.exportModel)
        self.ui.actionHelp.triggered.connect(self.help)
        self.ui.actionImport.triggered.connect(self.importModel)
        self.ui.actionModelBuilder.triggered.connect(self.modelBuilderDialog)
        self.ui.actionOpen.triggered.connect(self.importDialog)
        self.ui.actionSettings.triggered.connect(self.settingsDialog)
        self.ui.actionTextView.triggered.connect(self.changeModelView)
        self.ui.actionTreeView.triggered.connect(self.changeModelView)
        self.ui.actionUpdate.triggered.connect(self.getUpdate)
        
        # main widget
        self.ui.tabWidget.currentChanged.connect(self.changeMainTab)
        self.ui.stackedWidgetBulk.currentChanged.connect(self.changeModelView)
        self.ui.stackedWidgetSlab.currentChanged.connect(self.changeModelView)
        self.ui.pushBulkToText.pressed.connect(self.changeModelView)
        self.ui.pushBulkToTree.pressed.connect(self.changeModelView)
        self.ui.pushSlabToText.pressed.connect(self.changeModelView)
        self.ui.pushSlabToTree.pressed.connect(self.changeModelView)
    
    # Show about dialog
    def about(self):
        """Display 'About' dialog"""
        text = __APP_DESCRIPTION__
        text += '\n\nAuthor: {0} \nEmail: {1}'.format(__APP_AUTHOR__, 
                                                      __APP_EMAIL__)
        text += '\n\nApp version: {0}'.format(__APP_VERSION__)
        text += '\n{0}'.format(__APP_COPYRIGHT__)
        text += '\n' + __APP_LICENSE__
        text += '\n\nPython: {0}'.format(__PYTHON__)
        text += '\nGUI frontend: {0} {1}'.format(__QT_TYPE__, 
                                                 QtCore.QT_VERSION_STR)

        msg = QtGui.QMessageBox.about(self, self.ui.windowTitle(), text)
    
    # Display about dialog
    def aboutQt(self):
        """Display Qt dialog"""
        QtGui.QApplication.aboutQt()
    
    # Report bug / email devs
    def contactDeveloper(self):
        QtGui.QDesktopServices.openUrl(QtCore.QUrl(
                str("mailto: {email}?subject={name} feedback&body=").format(
                            email=__APP_EMAIL__, name=__APP_NAME__)))
    
    # check for update
    def getUpdate(self):
        """Check for app updates"""
        from UpdateDialog import UpdateDialog
        updateDialog = UpdateDialog(parent=self)
        updateDialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        updateDialog.exec_()
    
    # change main tab
    def changeMainTab(self):
        '''Change main tab selection'''
        tabText = str(self.ui.tabWidget.tabText(
                             self.ui.tabWidget.currentIndex())).lower()
        if tabText == 'bulk':  # bulk
            self.model = 'bulk'
        elif tabText == 'slab':  # slab
            self.model = 'slab'
        else:
            self.model = None
        
    # change view of model
    def changeModelView(self):
        '''Change model view'''
        if self.sender() is self.ui.actionTreeView:
            self.ui.actionTextView.setChecked(False)
            self.ui.stackedWidgetBulk.setCurrentIndex(0)
            self.ui.stackedWidgetSlab.setCurrentIndex(0)
        elif self.sender() is self.ui.actionTextView:
            self.ui.actionTextView.setChecked(False)
            self.ui.stackedWidgetBulk.setCurrentIndex(1)
            self.ui.stackedWidgetSlab.setCurrentIndex(1)
        elif (self.sender() is self.ui.pushBulkToText or
                self.sender() is self.ui.pushSlabToText):
            self.ui.actionTextView.setChecked(True)
            self.ui.actionTreeView.setChecked(False)
            self.ui.stackedWidgetBulk.setCurrentIndex(1)
            self.ui.stackedWidgetSlab.setCurrentIndex(1)
        elif (self.sender() is self.ui.pushBulkToTree or 
                self.sender() is self.ui.pushSlabToTree):
            self.ui.actionTextView.setChecked(False)
            self.ui.actionTreeView.setChecked(True)
            self.ui.stackedWidgetBulk.setCurrentIndex(0)
            self.ui.stackedWidgetSlab.setCurrentIndex(0)
    
    # export model as text file
    def exportModel(self):
        '''Export model as text file'''
        pass
    
    def importDialog(self):
        '''Open dialog and radio options'''
        importDialog = ImportDialog(parent=self, 
                            model=str(self.ui.tabWidget.tabText(
                                    self.ui.tabWidget.currentIndex())).lower())
        importDialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        importDialog.finished.connect(self.parseInput)
        importDialog.exec_()
    
    # import model from text file
    def importModel(self):
        '''Import model from text file'''
        pass
    
    def help(self):
        """Display help"""
        try:
            helpDialog = Help.HelpWidget(parent=self)
            helpDialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
            helpDialog.exec_()
            
        except NameError:
            QtGui.QMessageBox.information(self, 'Help', 
                                    'Help is not currently available')
            self.logger.error('unable to create Help dialog')

    # model builder
    def modelBuilderDialog(self):
        '''Start new instance of model builder wizard'''
        pass

    def getInputFile(self, startpath=str(
                        QtGui.QDesktopServices.storageLocation(
                            QtGui.QDesktopServices.HomeLocation)), model=None):
        '''returns file path of input'''
        if model == None:
            model = ''
        else:
            model += ' '
        
        model = model.capitalize()
        
        # start at last known directory
        if os.path.exists(self.lastpath):
            if os.path.isfile(self.lastpath):
                startpath = os.path.dirname(self.lastpath)
            else:
                startpath = self.lastpath
        
        filepath = str(QtGui.QFileDialog.getOpenFileName(parent=None, 
                     caption='Open %sInput File' % model, directory=startpath))
        
        return filepath
        
    # check type of input file
    def parseInput(self):
        if isinstance(self.sender(), ImportDialog):
            # check user did not abort
            if self.sender().action == 'cancel':
                print('cancel') 
                return
            
            # determine file type
            if self.sender().ui.radioBulk.isChecked():
                model = 'bulk'
            else:
                model = 'slab'
        
        else:  # guess from active tab
            tabText = str(self.ui.tabWidget.tabText(
                                self.ui.tabWidget.currentIndex())).lower()
            if tabText == 'bulk' or tabText == 'slab':
                model = tabText
            else:  # unknown
                return self.importDialog()  # start dialog
        
        filename = self.getInputFile(model=model)

        if not os.path.exists(filename):
            return  # user aborted
        else:
            self.lastpath = filename
        
        try:
            atom = Atom('H')  # dummy atom
            uc = Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) 
            mtz = MTZ_model(uc, atoms=[atom])  # initialise muffin-tin model
            mtz.load_from_file(filename)  # load file
            exec('self.%s = mtz' % model)
            self.updateModelUi(model)
            
        except IOError:
            self.logger.error("IOError: Unable to open input file '%s'" 
                              % filename)
        
    def settingsDialog(self):
        """Launch settings dialog"""
        from SettingsDialog import SettingsDialog
        settingsDialog = SettingsDialog(parent=self)
        settingsDialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        settingsDialog.finished.connect(self.updateSettings)
        settingsDialog.exec_()
    
    def updateModelUi(self, model=None):
        """Update model in gui""" 
        if isinstance(model, str):
            model = model.lower()
            mtz = MTZ_model(Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]), 
                            atoms=[Atom('H')])
            if model == 'bulk':
                tree = self.ui.treeWidgetBulk
                mtz = eval('self.%s' % model)
            elif model == 'slab':
                tree = self.ui.treeWidgetSlab
            else:
                return
            
            root = tree.invisibleRootItem()
            trunk = self.getChildItemsDict(tree)
            
            for branch in trunk:
                branch = branch.replace(' ', '').split('=')[0]
                item = root.child(self.treeRootDict.get(model).get(branch))
                if 'Name':
                    item.setText = mtz.header
                elif branch == 'Unitcell':
                    for i in range(3):
                        item.child(i).setText('a%i = %s' 
                                            % (i, mtz.coordinates[i]))
                elif branch == 'Atoms':
                    item.clear()
                    
                elif branch == 'Parameters':
                    params = self.getChildItemsDict(
                                    tree.topLevelItem(trunk.get(branch)))
                    parent = root.child(trunk.get(branch))
                    for param in params:
                        node = item.child(self.treeRootDict.get(model).get(branch))
                        if param == 'nh':  # update nh
                            node.setText(0, 'nh = %s' % str(mtz.nh))
                        elif param == 'Output':  # update output type
                            node.setText(0, 'Output = %s' % str(mtz.nform))
                        elif param == 'Exchange':  # update alpha
                            node.setText(0, 'Exchange = %s' % str(mtz.exchange))
    
    def updateSettings(self):
        '''update the application settings'''
        self.settings = self.sender().settings
        print(self.settings.__dict__)

    def getChildItemsDict(self, obj):
        try:
            if isinstance(obj, QtGui.QTreeWidget):
                root = obj.invisibleRootItem()
            elif isinstance(obj, QtGui.QTreeWidgetItem):
                root = obj
            child_count = root.childCount()
            topLevelDict = {}
            for i in range(child_count):
                item = root.child(i)
                var = str(item.text(0))
                exec('%s = i' % var)
                topLevelDict.update({var: eval(var)})
            return topLevelDict
        except any as e:
            self.logger.error(e.msg)
            
    def getChildItemHandle(self, obj, name=str):
        if isinstance(obj, QtGui.QTreeWidget):
            root = obj.invisibleRootItem()
        elif isinstance(obj, QtGui.QTreeWidgetItem):
            root = obj
        
        if isinstance(name, int):
            return root.child(name)
        elif isinstance(name, str):
            for i in range(root.childCount()):
                item = root.child(i)
                if str(item.text(0)) == name:
                    return item 
            
# boilerplate function - should be applicable to most applications
def main(argv=None):
    '''Entry point if executing as standalone''' 
    if argv is None:
        argv = sys.argv
        
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    
    return app.exec_()

# Execute main function if running as standalone module
if __name__ == '__main__':
    main()
