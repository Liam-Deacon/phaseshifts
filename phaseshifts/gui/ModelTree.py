#!/usr/bin/python
# -*- coding: utf-8 -*
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.m.deacon@gmail.com                                           #
#                                                                            #
# Copyright: Copyright (C) 2014-2015 Liam Deacon                             #
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
'''

'''
from __future__ import absolute_import, division, with_statement
from __future__ import print_function, unicode_literals

import os
import re
import sys

from qtsix import QtCore, QtGui
from qtsix.QtCore import QItemSelectionModel
from qtsix.QtGui import QIcon, QDesktopServices
from qtsix.QtWidgets import (QAction, QApplication, QCheckBox,
                             QFileDialog, QInputDialog, QLineEdit,
                             QMenu, QMessageBox, QTreeWidget, QTreeWidgetItem)


try:
    import res_rc
except ValueError:
    try:
        sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
        import res_rc
    except ImportError:
        import phaseshifts.gui.res_rc as res_rc

try:
    from ..model import Model, Atom
except ValueError:
    try:
        sys.path.insert(
            0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
        from model import Model, Atom
    except (ValueError, ImportError):
        from phaseshifts.model import Model, Atom


class ProjectTreeWidget(QTreeWidget):
    last_project_dir = os.path.expanduser('~')  # Project.default_dir
    default_dir = os.path.join(
        os.path.expanduser('~'), 'phaseshifts', 'models')

    def __init__(self, parent=None):
        super(ProjectTreeWidget, self).__init__(parent)

        self.setColumnCount(1)
        self.setHeaderLabel("Projects")

        # explorer actions
        self.renameAction = QAction(QIcon(":/tag_stroke.svg"),
                                    "&Rename", self,
                                    triggered=self.rename)
        self.renameAction.setToolTip("Rename project...")

        self.refreshAction = QAction(QIcon(":/spin.svg"),
                                     "Refresh", self,
                                     triggered=self.refresh,
                                     shortcut="F5")
        self.refreshAction.setToolTip("Refresh")

        self.newProjectAction = QAction(
            QIcon(":/document_alt_stroke.svg"),
            "&New Project", self,
            triggered=self.newProject)
        self.newProjectAction.setToolTip("Create new project...")

        self.importProjectAction = QAction(QIcon(":/import.svg"),
                                           "&Import Project", self,
                                           triggered=self.importProject)
        self.importProjectAction.setToolTip("Import existing project...")

        self.newModelAction = QAction(QIcon(":/atoms.svg"),
                                      "&New Model", self,
                                      triggered=self.newModel)
        self.newModelAction.setToolTip("Create new model...")

        self.importModelAction = QAction(QIcon(":/import.svg"),
                                         "&Import Model", self,
                                         triggered=self.importModel)
        self.importModelAction.setToolTip("Import existing model...")

        self.removeProjectAction = QAction(QIcon(":/x.svg"),
                                           "&Remove Project", self,
                                           triggered=self.removeProject,
                                           shortcut='Del')
        self.newProjectAction.setToolTip("Remove project")

        self.openProjectAction = QAction(
            QIcon(":/folder_fill.svg"),
            "Open Project &Location", self,
            triggered=self.openProjectLocation)
        self.newProjectAction.setToolTip(
            "Opens project location in file explorer")

        # explorer menus
        self.explorerDefaultMenu = QMenu()
        self.explorerDefaultMenu.addAction(self.newProjectAction)
        self.explorerDefaultMenu.addAction(self.importProjectAction)
        self.explorerDefaultMenu.addSeparator()
        # self.explorerDefaultMenu.addAction(self.copyAction)
        # self.explorerDefaultMenu.addAction(self.cutAction)
        # self.explorerDefaultMenu.addAction(self.pasteAction)
        self.explorerDefaultMenu.addAction(self.renameAction)
        self.explorerDefaultMenu.addSeparator()
        self.explorerDefaultMenu.addAction(self.refreshAction)

        self.explorerProjectMenu = QMenu()
        self.explorerProjectMenu.addAction(self.newModelAction)
        self.explorerProjectMenu.addAction(self.importModelAction)
        # self.explorerProjectMenu.addSeparator()
        # self.explorerProjectMenu.addAction(self.copyAction)
        # self.explorerProjectMenu.addAction(self.cutAction)
        # self.explorerProjectMenu.addAction(self.pasteAction)
        self.explorerProjectMenu.addAction(self.renameAction)
        self.explorerProjectMenu.addAction(self.openProjectAction)
        self.explorerProjectMenu.addSeparator()
        self.explorerProjectMenu.addAction(self.removeProjectAction)
        self.explorerProjectMenu.addSeparator()
        self.explorerProjectMenu.addAction(self.refreshAction)

        self.explorerFileMenu = QMenu()
        # self.explorerFileMenu.addAction(self.newAction)
        self.explorerFileMenu.addAction(self.refreshAction)

        # setup signals and slots
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.explorerPopupMenu)

        # recent projects
        self.recent_projects = []

    def openProjectLocation(self):
        ''' Opens the currently selected project's location '''
        import webbrowser
        filepath = self.currentProject()[
            'item'].project_dir or self.default_dir
        webbrowser.open(filepath)

    def getChildren(self, index):
        if not index.isValid():
            return []

        children = []
        for i in range(index.model().rowCount(index)):
            children += self.getChildren(index.child(i, 0))

        return children

    def expandChildren(self, index):
        ''' Recursely expands all children for the given index node'''
        if not index.isValid():
            return

        child_count = index.model().rowCount(index)
        for i in range(child_count):
            _child = index.child(i, 0)
            # Recursively call the function for each child node.
            self.expandChildren(_child)

        if not self.view.expanded(index):
            self.view.expand(index)

    def explorerPopupMenu(self, point):
        ''' Handles popup menu for explorer widget '''
        index = self.indexAt(point)
        if index.isValid():
            location = self.viewport().mapToGlobal(point)
            # show custom menu for file type held at given index
            item = self.itemFromIndex(index)
            if self.indexOfTopLevelItem(item) > -1:
                # then its a top-level item
                self.selectionModel().setCurrentIndex(index,
                                                      QItemSelectionModel.NoUpdate)
                self.explorerProjectMenu.popup(location)
            else:
                try:
                    self.currentItem().contextMenu.popup(location)
                    print('Handled right click of {} "{}"'
                          ''.format(self.currentItem(),
                                    self.currentItem().text(0)))
                except AttributeError:
                    print('Not handled right click of {} "{}"'
                          ''.format(self.currentItem(),
                                    self.currentItem().text(0)))
        else:
            # provide default menu
            self.explorerDefaultMenu.popup(
                self.viewport().mapToGlobal(point))

    def currentProject(self):
        '''returns the currently selected project'''
        item = self.currentItem()

        # get root item
        while self.indexOfTopLevelItem(item) < 0:
            item = item.parent()

        project = {'name': item.text(self.currentColumn()),
                   'item': item,
                   'index': self.currentIndex()}
        return project

    def newProject(self, projectName=None):
        if not projectName:
            projectName = "Untitled_Project"

        # get storage location for project
        homePath = QDesktopServices.storageLocation(
            QDesktopServices.HomeLocation)
        projectDir = os.path.join(homePath, "CLEED", "models")
        if not os.path.exists(projectDir):
            projectDir = self.examples_dir
        dlg = QFileDialog
        folder = QFileDialog.getExistingDirectory(parent=self,
                                                  caption="Select Project Base Directory",
                                                  directory=projectDir,
                                                  options=dlg.ShowDirsOnly |
                                                  dlg.DontResolveSymlinks)
        if folder:
            # do stuff
            items = [self.parent().topLevelItem(i).Path for i
                     in range(self.parent().topLevelItemCount())]
            if folder not in items:
                proj = ProjectItem(self.ui.parent(), path=folder)
            else:
                pass
                self.parent().setCurrentIndex(0, items.index(folder, ))

    def newModel(self, project, modelName=None):
        if not modelName:
            text, ok = QInputDialog.getText(self, 'Input Dialog',
                                            'Enter model name:')
            if not ok:
                return

            modelName = text

        try:
            index = self.selectedIndexes()[0]
            parent = self.itemFromIndex(index)
            path = os.path.join(parent.project_path, modelName)

            if not modelName:
                modelName = "New_Model"
                i = 1
                path = os.path.join(parent.Path, modelName)
                while os.path.isdir(modelName):
                    modelName = "New_Model%i" % i
                    path = os.path.join(parent.Path, modelName)
                    i += 1

            # a = parent.addChild(model)
            if not os.path.exists(path):
                os.makedirs(path, 755)
                # add new input files

            else:
                pass

        except IndexError:
            # no index selected (or created?)
            pass

    def importModel(self):
        ''' Import model from text file '''
        project = self.currentProject()['item']
        model = QFileDialog.getOpenFileName(parent=self,
                                            caption="Select CLEED project directory...",
                                            directory=project.project_dir,
                                            filter="*")
        if (os.path.exists(model) and
                model not in self.projects[project].models):
            self.projects[project].models.append(ModelItem.load(model))

    def importProject(self):
        ''' Import a project '''
        project = QFileDialog.getExistingDirectory(parent=self,
                                                   caption="Select CLEED "
                                                   "project directory...")
        if os.path.isdir(project) and project not in self.projects:
            self.projects.append(project)

    def removeProject(self):
        ''' Removes current project '''
        project = self.currentProject()

        msg_box = QMessageBox()
        msg_box.delete = QCheckBox("Delete project from filesystem",
                                   parent=msg_box)

        reply = msg_box.question(self, "Remove Project",
                                 "Are you sure you want to remove "
                                 "'{}'?".format(project['name']),
                                 QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            item = project['item']
            children = item.takeChildren()
            for child in children:
                del(child)
            self.takeTopLevelItem(self.indexOfTopLevelItem(item))
            del(item)

    def rename(self):
        ''' Renames current project '''
        project = self.currentProject()
        old_name = project['name']
        new_name, ok = QInputDialog.getText(self,
                                            self.tr("Rename Project"),
                                            self.tr("New name:"),
                                            QLineEdit.Normal,
                                            old_name)
        if ok and new_name is not old_name:
            item = project['item']
            item.setProject(dirname=item.project_dir, name=new_name,
                            overwrite=None, rename=True)

    def refresh(self):
        raise NotImplementedError('todo')

    def getChildItemsDict(self, obj):
        try:
            if isinstance(obj, QTreeWidget):
                root = obj.invisibleRootItem()
            elif isinstance(obj, QTreeWidgetItem):
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
        if isinstance(obj, QTreeWidget):
            root = obj.invisibleRootItem()
        elif isinstance(obj, QTreeWidgetItem):
            root = obj

        if isinstance(name, int):
            return root.child(name)
        elif isinstance(name, str):
            for i in range(root.childCount()):
                item = root.child(i)
                if str(item.text(0)) == name:
                    return item


class BaseItem(QTreeWidgetItem):

    def __init__(self, parent=None):
        super(BaseItem, self).__init__(parent)

        self.editAction = QAction(QIcon(":/document_edit_24x32.png"),
                                  "&Edit", None)
        self.refreshAction = QAction(QIcon(":/spin.svg"),
                                     "&Refresh", None,
                                     triggered=self.refresh,
                                     shortcut="F5")
        self.refreshAction.setToolTip("Refresh")

        self.contextMenu = QMenu()
        self.contextMenu.addAction(self.editAction)
        self.contextMenu.addAction(self.refreshAction)

    @classmethod
    def getChildren(cls, parent, recursive=True):
        ''' Get either immediate or all children of parent node '''
        children = []
        for i in range(parent.childCount()):
            child = parent.child(i)
            children += [child]
            if recursive:
                children += BaseItem.getChildren(child, recursive)
        return children

    def doubleClicked(self, index):
        """ Triggers doubleClick event """
        old_value = ''
        new_value, ok = QInputDialog.getText(self,
                                             self.tr("Rename"),
                                             self.tr("New value:"),
                                             QLineEdit.Normal,
                                             old_value)
        if ok and new_value is not old_value:
            item = self.parent.selectedIndexes()[0].model().itemFromIndex(
                index)
            item.setText(self.currentColumn(), new_value)

    def edit(self):
        ''' Edits the current item bringing up dialog if needed '''
        print("Edit of {}".format(self))

    def refresh(self):
        ''' Refreshes item and its children based on contained data '''
        print("Refresh of {}".format(self))


class ProjectItem(BaseItem):
    __new_name = 'NewProject'
    projects = []

    '''class for project items'''

    def __init__(self, parent=None, path=None, name=None, overwrite=None):
        super(ProjectItem, self).__init__(parent)
        self.setIcon(0, QIcon(":/folder_fill.svg"))
        # self.setFlags(self.flags() | QtCore.Qt.ItemIsEditable)
        self.setToolTip(0, "LEED-IV Project")
        self.models = []
        self.setProject(path, name, overwrite=overwrite)

        # add children
        self._init_children()

        ProjectItem.projects.append(self)

    def _init_children(self):
        pass

    def __del__(self):
        try:
            ProjectItem.projects.remove(self)
        except:
            pass

    @classmethod
    def load(cls, directory):
        pass

    @property
    def project_dir(self):
        if not hasattr(self, "_dir"):
            self._dir = None
        return self._dir or ProjectTreeWidget.default_dir

    @project_dir.setter
    def project_dir(self, dirname):
        self._dir = re.sub("[?<>\[\]{}&*^%$£@?|`~_]+",
                           "_", dirname or '') or ProjectTreeWidget.default_dir

    @property
    def project_path(self):
        return os.path.join(self.project_dir, self.name)

    def setProject(self, dirname, name, overwrite=None, rename=False):
        """ Sets the project directory and name """
        # old_name = self.name
        # old_dir = self.project_dir
        old_path = self.project_path

        self.name = name
        self.project_dir = dirname
        filepath = self.project_path

        # ask user to overwrite if not default specified
        if not os.path.isdir(self.project_dir):
            reply = QMessageBox.question(None, "Create missing directories",
                                         "Do you sure you want to create "
                                         "'{}'?".format(self.project_dir),
                                         QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.Yes:
                os.makedirs(self.project_dir, 755)

        # ask user to overwrite if not default specified
        if overwrite is None and os.path.exists(filepath):
            reply = QMessageBox.question(None, "{} exists".format(self.name),
                                         "Do you sure you want to overwrite '{}'?".format(
                                             filepath),
                                         QMessageBox.Yes | QMessageBox.No)
            overwrite = True if reply == QMessageBox.Yes else False

        # split number from end and enumerate as needed
        i = re.sub(".*?([0-9]+)$", "\\1", self.name)
        base_path = re.split("" + i + "$",  filepath)[0]
        i = int(i) if i.isdigit() else 0

        while os.path.exists(filepath) and not overwrite:
            filepath = "{}{:03d}".format(base_path, i)
            i += 1

        self.name = os.path.basename(filepath)
        self.setText(0, self.name)
        self.setToolTip(0, 'LEED-IV Project: "{}"'.format(self.project_path))
        self.setData(0, QtCore.Qt.UserRole, {'project': self.project_path,
                                             'models': self.models})

        if os.path.exists(old_path) and rename:
            from shutil import move
            try:
                move(old_path, self.project_path)
            except (IOError, OSError) as err:
                sys.stderr.write(err.message + '\n')

    @classmethod
    def newProject(cls, parent=None,
                   path=ProjectTreeWidget.default_dir,
                   name=__new_name, overwrite=None):
        """ Creates a new ProjectItem instance """
        return cls(parent, path, name, overwrite)

    @property
    def name(self):
        if not hasattr(self, "_name"):
            self._name = None
        return self._name or self.__new_name

    @name.setter
    def name(self, name):
        # ensure name is valid
        self._name = re.sub("[?<>\[\]{}&*^%$£@?|`~_]+",
                            "_", name or '') or self.__new_name


class ModelItem(BaseItem):
    ''' common class for both bulk and surface model items '''
    MODEL_CLASS = Model
    __tooltip__ = "Model"
    __description__ = ""

    def __init__(self, parent=None, model=None):
        super(BaseItem, self).__init__(parent)
        # self.setFlags(self.flags() | QtCore.Qt.ItemIsEditable)

        self.model = model

    @QtCore.Slot()
    def modelChanged(self, model):
        """ Actions to undertake once model has changed """
        print(self.__tooltip__ + " changed")

    def setModel(self, model):
        ''' sets the model '''
        if isinstance(model, basestring):
            # attempt to load model assuming filepath
            model = self.MODEL_CLASS().load(model)
        elif model == None:
            model = self.MODEL_CLASS()
        else:
            try:
                self.MODEL_CLASS(**model)
            except:
                raise ValueError(
                    "{} is not a supported model type".format(model))
        self._model = model

    def getModel(self):
        ''' gets the model '''
        return self._model

    def refresh(self):
        ''' Updates QTreeWidgetItem '''
        try:
            filename = self.model
            self.setToolTip(0, self.__tooltip__ +
                            ' ({})\n\n'.format(filename) + self.__description__)
        except:
            self.setToolTip(0, self.__tooltip__ + '\n' + self.__description__)

        print("TODO: clear all children")
        if self.model:
            atoms = [AtomItem(atom=a) for a in self.model.atoms]

    model = property(fset=setModel, fget=getModel)


class SettingsItem(BaseItem):
    '''class for local settings'''

    def __init__(self, parent=None, path=None):
        super(self.__class__, self).__init__(parent)
        self.setIcon(0, QIcon(":/wrench.svg"))
        self.setText(0, "Settings")
        self.setToolTip(0, "Defines settings for project")


class AtomItem(BaseItem):
    ''' Class for handling atoms '''

    def __init__(self, parent=None, atom=None):
        super(self.__class__, self).__init__(parent)

        self.atom = atom or self.newAtom()
        self.refresh()

    @QtCore.Slot()
    def atomChanged(self, atom):
        pass

    def newAtom(self, element='C', **kwargs):
        return Atom(element, **kwargs)

    def refresh(self):
        self.setIcon(0, QIcon(":/atom.svg"))
        self.setText(0, self.atom.symbol)
        self.setToolTip(
            0, "{} atom\n{}".format(self.atom.name, str(self.atom)))

        for i, j in enumerate(["x", "y", "z"]):
            eval("self.{} = QTreeWidgetItem(self)".format(j))
            eval('self.{}.setText(0, "{} = {}")'.format(
                j, j, self.atom.coordinates[i]))
            eval(
                'self.{}.setToolTip(0, "{}-coordinate for atom"'.format(j, j.upper()))

        # valence
        self.valence = QTreeWidgetItem(self)
        self.valence.setText(0, "valence = {}".format(self.atom.valence))
        self.valence.setToopTip(
            0, "Specifies the valency (charge) of the atom (or ion) e.g. +2")

        # fractional occupancy
        self.occupancy = QTreeWidgetItem(self)
        self.occupancy.setText(0, "occupancy = {}".format(self.atom.occupancy))
        self.occupancy.setToolTip(0, "Specifies the fractional occupancy "
                                  "of the atomic site.\nThis is useful if "
                                  "sites (e.g. on the surface) are known to \n"
                                  "contain vacancies or the structure is an "
                                  "alloy")
        # minimum radius
        self.radius = QTreeWidgetItem(self)
        self.radius.setText(0, "r_{min} = {}".format(self.atom.radius))
        self.radius.setToolTip(0, "Specifies the muffin-tin radius of "
                               "the atom or ion.\nThis radius is the smallest "
                               "interatomic distance before \na penalty is "
                               "invoked during geometry optimisation")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    explorer = ProjectTreeWidget()

    project = ProjectItem.newProject(overwrite=True)

    explorer.addTopLevelItem(project)

    project.setExpanded(True)
    for child in project.getChildren(project):
        child.setExpanded(True)

    explorer.show()
    app.exec_()
