#!/usr/bin/env python3
#
# The DREAM Theater GUI
#

from . import resolvedreampaths

from PyQt5 import QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox
from .ui import DREAMTheater_design

from DREAM import DREAMSettings, DREAMOutput
import dreampyface
from dreampyface import Simulation
import sys


class DREAMTheater(QtWidgets.QMainWindow):
    

    def __init__(self, settings):
        """
        Constructor.
        """
        QtWidgets.QMainWindow.__init__(self)

        self.ui = DREAMTheater_design.Ui_DREAMTheater()
        self.ui.setupUi(self)

        # Construct a simulation object
        self.simulation = Simulation(settings)

        self.treeViewModel = QtGui.QStandardItemModel()
        self.ui.treeView.setModel(self.treeViewModel)

        self.loadTreeView()

        self.bindEvents()


    def bindEvents(self):
        """
        Bind functions to control events.
        """
        self.ui.actionExit.triggered.connect(self.exit)
        pass


    def closeEvent(self, event):
        self.exit()


    def exit(self):
        self.close()


    def loadTreeView(self):
        """
        Load the list of unknowns in the simulation.
        """
        self.treeViewModel.clear()

        self.unknowns = self.simulation.unknowns.getInfo()
        unk = QtGui.QStandardItem('Unknown quantities')

        for uname in self.unknowns.keys():
            unk.appendRow(QtGui.QStandardItem(uname))

        self.treeViewModel.invisibleRootItem().appendRow(unk)



if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)

    win = DREAMTheater('../../examples/runaway/dream_settings.h5')
    win.show()
    sys.exit(app.exec_())


