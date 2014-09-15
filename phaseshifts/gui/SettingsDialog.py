'''
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
'''
import os
import sys
from time import gmtime, strftime
from PyQt4 import QtGui, uic
import res_rc

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

from settings import Settings


class SettingsDialog(QtGui.QDialog):
    '''
    Dialog class for updating sequences 
    '''
    def __init__(self, parent=None):
        super(SettingsDialog, self).__init__(parent)

        # dynamically load ui
        self.ui = uic.loadUi("gui/SettingsDialog.ui", self)
        
        # settings default
        self.settings = Settings()
        
        # initialise ui
        self.initUi()
        
        self.ui.show()
            
    def initUi(self):
        '''Initialise user interface'''
        # Setup slots and signals
        self.ui.spinGrid.valueChanged.connect(self.updateGrid)
        self.ui.spinMixing.valueChanged.connect(self.updateMixing)
        self.ui.spinTolerence.valueChanged.connect(self.updateTolerence)
        self.ui.spinEch.valueChanged.connect(self.updateEch)
        self.ui.spinRelic.valueChanged.connect(self.updateRelic)
        self.ui.spinExchange.valueChanged.connect(self.updateExchange)
        self.ui.checkRelativistic.stateChanged.connect(self.updateRelativistic)
        self.ui.comboMethod.currentIndexChanged.connect(self.updateMethod)
        self.ui.checkWriteFiles.stateChanged.connect(self.updateWriteFiles)
        self.ui.checkDelete.stateChanged.connect(self.updateDeleteFiles)
        self.ui.checkLog.stateChanged.connect(self.updateLogging)
        self.ui.comboPath.currentIndexChanged.connect(self.updateWritePath)
        self.ui.lineWritePath.editingFinished.connect(self.updateWritePath)
        self.ui.toolWritePath.pressed.connect(self.getWritePath)
        self.ui.lineLogFile.editingFinished.connect(self.updateLogPath)
        self.ui.toolLogFile.pressed.connect(self.getLogPath)
        self.ui.buttonBoxSave.clicked[
                            QtGui.QAbstractButton].connect(self.adjustSettings)
        self.ui.buttonBoxApply.clicked[
                            QtGui.QAbstractButton].connect(self.adjustSettings)
                            
        self.ui.spinAlpha.valueChanged.connect(self.updateAlpha)
        self.ui.spinMuffinTinZero.valueChanged.connect(self.updateNh)
        self.ui.radioCavendish.clicked.connect(self.updatePhaseMethod)
        self.ui.radioWilliams.clicked.connect(self.updatePhaseMethod)
        self.ui.radioRelativistic.clicked.connect(self.updatePhaseMethod)
        self.ui.comboFormat.currentIndexChanged.connect(self.updateFormat)
        
        self.resetSettings()
        
    def adjustSettings(self, button):
        '''determine whether to apply settings and respond accordingly'''
        action = str(button.text()).lower()
        if action == 'apply':
            # apply settings, but do not close dialog
            self.applySettings()
        
        elif action == 'cancel':
            # do not apply settings & close dialog
            self.ui.close()
        
        elif action == 'reset':
            self.resetSettings()
        
        elif action == 'ok':    
            self.ui.close()
            
        elif action == 'open':
            # open settings ini file
            self.loadConfig()
    
        elif action == 'save':
            # save settings ini file
            self.saveConfig()
        
    def applySettings(self):
        '''update Settings class'''
        self.settings = Settings(
                grid=self.ui.spinGrid.value(),
                mixing=self.ui.spinMixing.value(),
                tolerence=self.ui.spinTolerence.value(),
                ech=self.ui.spinEch.value(),
                relic=self.ui.spinRelic.value(),
                atorb_exchange=self.ui.spinExchange.value(),
                relativistic_atorb=self.ui.checkRelativistic.isChecked(),
                write_files=self.ui.checkWriteFiles.isChecked(),
                delete_files_after_use=self.ui.checkDelete.isChecked(),
                filepath=str(self.ui.lineWritePath.text()),
                mufpot_exchange=self.ui.spinAlpha.value(),
                mufpot_nh=self.ui.spinMuffinTinZero.value(),
                phsh_method=self.phaseMethod,
                phsh_format=str(self.ui.comboFormat.currentText()),
                logfile=str(self.ui.lineLogFile.text()),
                logging=self.ui.checkLog.isChecked())

    def getLogPath(self):
        '''use file dialog to set path to log file'''
        logfile = str(QtGui.QFileDialog.getSaveFileName(self, 
                                    "Save Log File...", "log.txt", 
                                    'Text document (*.txt)'))

        if logfile != '':  # check for empty path
            if os.path.isfile(logfile):
                reply = QtGui.QMessageBox.warning(self, 'Warning!',
                            "File already exists! Do you wish to overwrite?", 
                            QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.No:
                    # get current time-stamp
                    timestamp = strftime("%Y-%m-%d_%H%M", gmtime()) 
                    logfile, extension = os.path.splitext(logfile)
                    # append time-stamp to filename
                    logfile = logfile + '_' + timestamp + '.' + extension
                    
            elif not os.path.exists(os.path.dirname(logfile)):  
                # no parent directory
                reply = QtGui.QMessageBox.warning(self, 'Warning!',
                            "'%s' directory does not exist! "
                            "Do you wish to create it?" 
                            % os.path.dirname(logfile), 
                            QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:  
                    # create parent directory
                    os.makedirs(os.path.dirname(logfile))  
                else:
                    return  # do not modify log file
            
            self.ui.lineLogFile.setText(logfile)
    
    def loadConfig(self):
        '''Attempts to loads the settings from valid ini file'''
        config = configparser.ConfigParser()
        default = os.path.join(str(QtGui.QDesktopServices.storageLocation(
                                QtGui.QDesktopServices.DataLocation)), 
                               '.phaseshifts', 'config.ini')
        filepath = str(QtGui.QFileDialog.getOpenFileName(self, 
                                'Open ini file...', default))
        
        if filepath == '':
            return  # user abort
        
        try:
            config.read(filepath)
        except IOError:
            QtGui.QMessageBox.critical(self, 'Error!', 
                       "Unable to read ini file '%s'" % filepath)
            return
        
        general_dict = {
                'write_files': 'self.ui.checkWriteFiles.setChecked(bool(%s))',
                'delete_files': 'self.ui.checkDelete.setChecked(bool(%s))',
                'path': 'self.ui.lineWritePath.setText(\"%s\")',
                'logging': 'self.ui.checkLog.setChecked(bool(%s))',
                'logfile': 'self.ui.lineLogFile.setText(\"%s\")'
                }
        
        atorb_dict = {
            'grid': 'self.ui.spinGrid.setValue(int(%s))',
            'mixing': 'self.ui.spinMixing.setValue(float(%s))',
            'tolerence': 'self.ui.spinTolerence.setValue(float(%s))',
            'ech': 'self.ui.spinEch.setValue(float(%s))',
            'relic': 'self.ui.spinRelic.setValue(float(%s))',
            'exchange': 'self.ui.spinExchange.setValue(float(%s))',
            'relativistic': 'self.ui.checkRelativistic.setChecked(bool(%s))',
            }
        
        phsh_dict = {
                     'alpha_exchange': 'self.ui.spinAlpha.setValue(float(%s))',
                     'phase_method': 'self.updatePhaseMethod(method=\"%s\")',
                     'nh': 'self.ui.spinMuffinTinZero.setValue(int(%s))',
                     'phase_format': 'self.updateFormat(formatting=\"%s\")'
                     }
        
        for option in general_dict:  # load each setting from general section
            try:
                value = config.get('GENERAL', option)
                line = general_dict.get(option)
                eval(line % value)
            except Exception as e:
                sys.stderr.write("Error loading %s in config file - %s\n" 
                                    % (option, e))
                sys.stderr.flush()
        
        for option in atorb_dict:  # load each setting from atorb section
            try:
                value = config.get('ATORB', option)
                line = atorb_dict.get(option)
                eval(line % value)
            except Exception as e:
                sys.stderr.write("Error loading %s in config file - %s\n" 
                                    % (option, e))
                sys.stderr.flush()
        
        for option in phsh_dict:  # load each setting from phaseshifts section
            try:
                value = config.get('PHASESHIFTS', option)
                line = phsh_dict.get(option)
                eval(line % value)
            except Exception as e:
                sys.stderr.write("Error loading %s in config file - %s\n" 
                                    % (option, e))
                sys.stderr.flush()
                
    def saveConfig(self):
        config = configparser.ConfigParser()
        default = os.path.join(str(QtGui.QDesktopServices.storageLocation(
                                QtGui.QDesktopServices.DataLocation)), 
                               '.phaseshifts', 'config.ini')
        filepath = str(QtGui.QFileDialog.getSaveFileName(self, 
                                    'Save ini file...', default))
        
        if filepath == '':
            return  # user abort
        
        # try to create parent directory if it doesn't exist
        if not os.path.exists(os.path.dirname(filepath)):
            try:
                os.makedirs(filepath)
            except IOError:
                QtGui.QMessageBox.critical(self, 'Error!', 
                            "Unable to create directory '%s'" % filepath)
                return
            
        config.add_section('GENERAL')
        config.add_section('ATORB')
        config.add_section('PHASESHIFTS')
        config.add_section('HISTORY')
        
        atorb_dict = {
            'grid': self.ui.spinGrid.value(),
            'mixing': self.ui.spinMixing.value(),
            'tolerence': self.ui.spinTolerence.value(),
            'ech': self.ui.spinEch.value(),
            'relic': self.ui.spinRelic.value(),
            'method': self.ui.comboMethod.currentText(),
            'exchange': self.ui.spinExchange.value(),
            'relativistic': self.ui.checkRelativistic.isChecked(),
            }
        
        general_dict = {
                'write_files': self.ui.checkWriteFiles.isChecked(),
                'delete_files': self.ui.checkDelete.isChecked(),
                'path': self.ui.lineWritePath.text(),
                'logging': self.ui.checkLog.isChecked(),
                'logfile': self.ui.lineLogFile.text()
                }
        
        phsh_dict = {
                     'alpha_exchange': self.ui.spinAlpha.value(),
                     'phase_method': self.phaseMethod,
                     'nh': self.ui.spinMuffinTinZero.value(),
                     'phase_format': self.ui.comboFormat.currentText()
                     }
        
        for option in atorb_dict:
            value = atorb_dict.get(option)
            config.set('ATORB', option, str(value))
            
        for option in general_dict:
            value = general_dict.get(option)
            config.set('GENERAL', option, str(value))

        for option in phsh_dict:
            value = phsh_dict.get(option)
            config.set('PHASESHIFTS', option, str(value))
            
        try:
            with open(filepath, 'w') as configfile:    # save
                config.write(configfile)
        except IOError:
            QtGui.QMessageBox.critical(self, 'Error!', 
                        "Unable to write settings file '%s'" % filepath)
    
    def getWritePath(self):
        '''use file dialog to set path to generated files'''
        path = str(QtGui.QFileDialog.getExistingDirectory(self, 
                    "Specify File Directory", self.ui.lineWritePath.text()))
        
        if path != '':  # check for empty path
            if not os.path.exists(path):
                reply = QtGui.QMessageBox.warning(self, 'Warning!',
                                "'%s' directory does not exist! "
                                "Do you wish to create it?" % os.pardir(path), 
                                QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:
                    try:
                        os.makedirs(path)
                    except IOError:
                        QtGui.QMessageBox.critical(self, 'Error!', 
                                    "Unable to create directory '%s'" % path)
                        return
                else:
                    return  # do not create directory, nor store new path
                
            self.ui.lineWritePath.setText(path)

    def resetSettings(self):
        # reset to original values
        self.ui.spinGrid.setValue(self.settings.grid)
        self.ui.spinMixing.setValue(self.settings.mixing)
        self.ui.spinTolerence.setValue(self.settings.tolerence)
        self.ui.spinEch.setValue(self.settings.ech)
        self.ui.spinRelic.setValue(self.settings.relic)
        self.ui.checkRelativistic.setChecked(
                                    self.settings.relativistic_atorb)
        
        # debugging
        self.ui.checkWriteFiles.setChecked(self.settings.write_files)
        self.ui.checkDelete.setChecked(
                                    self.settings.delete_files_after_use)
        self.ui.lineWritePath.setText(self.settings.filepath)
        self.ui.lineLogFile.setText(self.settings.logfile)
        self.ui.checkLog.setChecked(self.settings.logging)
        
        # muffin-tin potential and phase shifts 
        self.updatePhaseMethod(method=self.settings.phsh_method)
        self.updateFormat(formatting=self.settings.phsh_format)
        self.ui.spinAlpha.setValue(self.settings.mufpot_exchange)
        self.ui.spinMuffinTinZero.setValue(self.settings.mufpot_nh)


    def updateAlpha(self):
        '''update Hartree exchange term for muffin tin potential'''
        self.alpha = self.sender().value()

    def updateDeleteFiles(self):
        '''update whether to delete files when done'''
        self.deleteFiles = self.sender().isChecked()

    def updateEch(self):
        '''update atorb ech value'''
        self.ech = self.sender().value()

    def updateExchange(self):
        '''update exchange correlation value'''
        self.exchange = self.sender().value()
        
        if self.sender().value() == 0.0:
            self.ui.comboMethod.setCurrentIndex(0)
        elif self.sender().value() >= 0.0:
            self.ui.comboMethod.setCurrentIndex(1)
        elif self.sender().value() <= 0.0:
            self.ui.comboMethod.setCurrentIndex(2)
        
    def updateFormat(self, formatting=None):
        '''update formatting of phase shift files'''
        if isinstance(self.sender(), QtGui.QComboBox):
            self.format = self.sender().currentText()
            return

        if isinstance(formatting, str):
            for i in range(self.ui.comboFormat.count()):
                if formatting == self.ui.comboFormat.itemText(i):
                    self.ui.comboFormat.setCurrentIndex(i)
                    return
            
            sys.stderr.write("Invalid phaseshift format '%s'\n" % formatting)
            sys.stderr.flush()
        
    def updateGrid(self):
        '''update atorb grid value'''
        self.grid = self.sender().value()
        
    def updateLogging(self):
        '''updates whether to store logs to file or to stdout'''
        self.logging = self.sender().isChecked()
        self.ui.lineLogFile.setEnabled(self.logging)
        self.ui.toolLogFile.setEnabled(self.logging)

    def updateLogPath(self):
        '''update path to log file'''
        self.logfile = self.sender().text()
    
    def updateMethod(self):
        '''update calculation method (LDA, xalpha...)'''
        if self.sender().currentText() == 'Hartree-Fock':
            self.ui.spinExchange.setValue(0.0)
            
        elif self.sender().currentText() == 'LDA':
            self.ui.spinExchange.setValue(1.0)
            
        elif self.sender().currentText() == 'xalpha':
            self.ui.spinExchange.setValue(-1.0)
    
    def updateMixing(self):
        '''update atorb mixing value'''
        self.mixing = self.sender().value()
   
    def updateNh(self):
        '''update term for estimating muffin tin zero'''
        self.nh = self.sender().value()
        
    def updatePhaseMethod(self, method=None):
        '''update method for calculating phase shifts'''
        if isinstance(self.sender(), QtGui.QRadioButton):
            self.phaseMethod = str(self.sender().text())
            return
        
        # manual update using str
        if isinstance(method, str):
            method = str(method).capitalize()
            if method == str(self.ui.radioCavendish.text()):
                self.phaseMethod = method
                self.ui.radioCavendish.click()
            elif method == str(self.ui.radioWilliams.text()):
                self.phaseMethod = method
                self.ui.radioWilliams.click()
            elif method == str(self.ui.radioRelativistic.text()):
                self.phaseMethod = method
                self.ui.radioRelativistic.click()
        
    def updateRelativistic(self):
        '''update whether atorb calculations are relativistic'''
        self.relativistic = self.sender().isChecked()
        
    def updateRelic(self):
        '''update Relic level value'''
        self.relic = self.sender().value()
   
    def updateTolerence(self):
        '''update atorb tolerence value'''
        self.tolerence = self.sender().value()
       
    def updateWriteFiles(self):
        '''update whether to permenantly store intermediate files'''
        self.writeFiles = self.sender().isChecked()
        self.ui.checkDelete.setEnabled(self.writeFiles)
        self.ui.comboPath.setEnabled(self.writeFiles)
        self.ui.toolWritePath.setEnabled(self.writeFiles)
        self.ui.lineWritePath.setEnabled(self.writeFiles)

    def updateWritePath(self):
        '''update directory where intermediate files will be stored.
        This should be useful for debugging purposes'''
        path = self.ui.lineWritePath.text()
        if isinstance(self.sender(), QtGui.QComboBox):
            text = str(self.sender().currentText()).upper()
            app = 'phaseshifts'
            if text == 'DOCUMENTS':
                path = os.path.join(str(QtGui.QDesktopServices.storageLocation(
                            QtGui.QDesktopServices.DocumentsLocation)), app)
            elif text == 'APPDATA':
                path = os.path.join(str(QtGui.QDesktopServices.storageLocation(
                                QtGui.QDesktopServices.DataLocation)), app)
            elif text == 'TEMP':
                path = os.path.join(str(QtGui.QDesktopServices.storageLocation(
                                QtGui.QDesktopServices.TempLocation)), app)
            elif text == 'HOME':
                app = '.' + str(app)
                path = os.path.join(str(QtGui.QDesktopServices.storageLocation(
                                QtGui.QDesktopServices.HomeLocation)), app)
            else:  # enable custom path
                path = os.path.join('.', app)
                
        path = os.path.abspath(os.path.expandvars(path))
        
        # try to write path if it doesn't already exist
        if not os.path.exists(path):
            try:
                os.makedirs(path)
            except IOError:
                QtGui.QMessageBox.critical(self, 'Error!', 
                            "Unable to create directory '%s'" % path)
                return
            
        self.ui.lineWritePath.setText(path)
        self.writePath = path
