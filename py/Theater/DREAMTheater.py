#!/usr/bin/env python3
#
# The DREAM Theater GUI
#

from . import resolvedreampaths

from PyQt5 import QtGui, QtWidgets, Qt
from PyQt5.QtWidgets import QMessageBox, QFileDialog
from PyQt5.QtGui import QCursor
from PyQt5.QtCore import Qt, QPoint
from .ui import DREAMTheater_design

from DREAM import DREAMSettings, DREAMOutput
from . controls.CodeEditor import CodeEditor
from . DataProviderOutput import DataProviderOutput
from . DataProviderSimulation import DataProviderSimulation
from . import evaluateExpression
from . PlotConfiguration import PlotConfiguration
from . SimulationThread import SimulationThread
import dreampyface
from dreampyface import Simulation
import h5py
import matplotlib as mpl
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import traceback
import random
import sys

from . DialogFluidQuantity import DialogFluidQuantity

try:
    import webbrowser
except Exception as ex:
    print("Unable to load the 'webbrowser' module. Will not be able to direct you to excellent music.")
    print(ex)

# Disable use of TeX for rendering text
mpl.rc('text', usetex=False)


class DREAMTheater(QtWidgets.QMainWindow):
    
    
    simulationThread = None
    scriptFilename = None
    output = None


    def __init__(self, data):
        """
        Constructor.

        :param data: May be either a string specifying a file name, a Python dict containing simulation settings, a DREAMSettings object or a DREAMOutput object.
        """
        QtWidgets.QMainWindow.__init__(self)

        self.ui = DREAMTheater_design.Ui_DREAMTheater()
        self.ui.setupUi(self)

        # Construct a simulation object
        #self.simulation = Simulation(settings)
        if type(data) == str:
            # Determine if the named file is an input or output file...
            isInput = None
            with h5py.File(data, 'r') as f:
                isInput = 'init' in f

            if isInput:
                data = DREAMSettings(data)
            else:
                data = DREAMOutput(data)

        if type(data) == DREAMOutput:
            self.setOutput(data)
            self.data = DataProviderOutput(data)
        else:
            self.data = DataProviderSimulation(Simulation(data))
            self.ui.pbSimulation.setValue(0)

        self.enableUI()

        self.treeViewModel = QtGui.QStandardItemModel()
        self.ui.treeView.setModel(self.treeViewModel)

        self.loadTreeView()

        # Hide unused labels
        self.ui.lblEquation.setText('')

        # Add code editor
        self.editor = CodeEditor(self.ui.tabEditor, self)
        self.ui.verticalLayout_6.addWidget(self.editor)
        self.editor.executeCommand.connect(self.runScript)

        self.bindEvents()
        self.setupPlots()


    def bindEvents(self):
        """
        Bind functions to control events.
        """
        self.ui.actionExit.triggered.connect(self.exit)
        self.ui.treeView.customContextMenuRequested.connect(self.plottableContextMenu)
        self.ui.treeView.selectionModel().selectionChanged.connect(self.plottableSelected)
        self.ui.btnRunSimulation.clicked.connect(self.runSimulation)

        self.ui.btnRunScript.clicked.connect(self.runScript)
        self.ui.btnLoadScript.clicked.connect(self.loadScript)
        self.ui.btnSaveScript.clicked.connect(self.saveScript)

        self.ui.actionDocumentation.triggered.connect(self.openDocumentation)
        self.ui.actionDREAM.triggered.connect(self.openDREAM)
        self.ui.actionThe_DREAM_Theater.triggered.connect(self.openDREAMTheater)
        self.ui.actionDream_Theater_2.triggered.connect(self.openMusic)

        self.editor.document().modificationChanged.connect(self.onScriptModified)


    def closeEvent(self, event):
        self.exit()

    
    def enableUI(self):
        """
        Enable and disable appropriate parts of the UI, depending on
        if a DREAMOutput object is loaded or a Simulation object is used.
        """
        out = not self.isSimulation()

        simOnly = not out
        outOnly = out

        self.ui.tabSimulation.setEnabled(simOnly)


    def enableEditor(self, enable=True):
        """
        Enable/disable script editor controls.
        """
        self.ui.btnRunScript.setEnabled(enable)
        self.ui.btnLoadScript.setEnabled(enable)
        self.ui.btnSaveScript.setEnabled(enable)
        self.editor.setEnabled(enable)


    def enableSimulationControls(self, enable=True):
        """
        Enable/disable simulation controls.
        """
        self.ui.btnRunSimulation.setEnabled(enable)


    def exit(self):
        self.close()


    def isDREAMOutput(self):
        """
        Returns ``True`` if the current simulation uses a
        DREAMOutput object, rather than a libdreampy Simulation
        object.
        """
        return (type(self.data) == DataProviderOutput)


    def isSimulation(self):
        """
        Returns ``True`` if the current simulation uses a
        libdreampy Simulation object, rather than just a DREAMOutput
        object.
        """
        return (type(self.data) == DataProviderSimulation)


    def loadTreeView(self):
        """
        Load the list of unknowns in the simulation.
        """
        self.treeViewModel.clear()

        # UNKNOWNS
        #self.unknowns = self.simulation.unknowns.getInfo()
        self.unknowns = self.data.getUnknownInfo()
        unk = QtGui.QStandardItem('Unknown quantities')

        for uname in self.unknowns.keys():
            itm = QtGui.QStandardItem(uname)
            itm.setData('unknown')
            unk.appendRow(itm)

        self.treeViewModel.invisibleRootItem().appendRow(unk)

        # OTHER QUANTITIES
        self.others = self.data.getOtherInfo()
        oth = QtGui.QStandardItem('Other quantities')

        def addItem(path, name, grps, treeItem):
            s, _, n = name.partition('/')
            fullname = path+'/'+s
            
            if n == "":
                itm = QtGui.QStandardItem(s)
                itm.setData('other')
                treeItem.appendRow(itm)
            else:
                if fullname not in grps:
                    t = QtGui.QStandardItem(s)
                    grps[fullname] = t
                    treeItem.appendRow(t)
                else:
                    t = grps[fullname]

                addItem(fullname, n, grps, t)

            return grps

        groups = {}
        for oname in self.others.keys():
            groups = addItem('/', oname, groups, oth)

        self.treeViewModel.invisibleRootItem().appendRow(oth)


    def loadScript(self):
        """
        Load a script from a file.
        """
        if self.scriptHasBeenModified():
            resp = QMessageBox.question(self, "Script contains unsaved changes",
                "The current script contains unsaved changes. Would you like to save it before opening a new script?",
                QMessageBox.Yes | QMessageBox.Discard | QMessageBox.Cancel)
            
            if resp == QMessageBox.Yes:
                self.saveScript()
            elif resp == QMessageBox.Cancel:
                return

        filename, _ = QFileDialog.getOpenFileName(parent=self, caption="Load Python script", filter="Python script (*.py);;All files (*.*)")

        if filename:
            self.scriptFilename = filename

            with open(filename, 'r', encoding='utf8') as f:
                text = f.read()

                self.editor.setPlainText(text.replace("\t", "    "))
                self.setScriptModified(False)
                self.editor.setFocus()


    def getPlottablePath(self, item):
        """
        Returns the treeview "path" to the given item.
        """
        if item == None:
            return ''
        else:
            return self.getPlottablePath(item.parent()) + '/' + item.text()


    def getSelectedPlottable(self):
        """
        Returns the plottable item currently selected in the TreeView.
        """
        return self.treeViewModel.itemFromIndex(self.ui.treeView.selectionModel().selectedIndexes()[0])


    def onFigureClicked(self, event):
        """
        Triggered when the figure is clicked.
        """
        if event.guiEvent.button() != Qt.RightButton:
            return
        if event.inaxes is None:
            return

        menu = QtWidgets.QMenu()

        # Monitor
        actionRemove = menu.addAction("Remove")
        
        #print(self.ui.framePlot.mapToGlobal(QPoint(event.x, event.y)))
        res = menu.exec_(self.ui.framePlot.mapToGlobal(QPoint(event.x, event.y)))
        if res == actionRemove:
            self.plotconfig.removeByAxes(event.inaxes)
            self.plotconfig.render(self.output, clearAxes=True)


    def openDocumentation(self):
        webbrowser.open_new_tab('https://ft.nephy.chalmers.se/dream')


    def openDREAM(self):
        webbrowser.open_new_tab('https://github.com/chalmersplasmatheory/DREAM')


    def openDREAMTheater(self):
        webbrowser.open_new_tab('https://ft.nephy.chalmers.se/dream/frontend/index.html#the-dream-theater-gui')


    def openMusic(self):
        """
        Direct the user to some excellent music.
        """
        selection = [
            # Pull Me Under
            'https://www.youtube.com/watch?v=mipc-JxrhRk',
            # Octavarium
            'https://www.youtube.com/watch?v=XYV8Zt2k0RQ',
            # The Enemy Inside
            'https://www.youtube.com/watch?v=RoVAUUFjl0I',
            # On The Backs of Angels
            'https://www.youtube.com/watch?v=28MmnThlYOo',
            # The Dance of Eternity
            'https://www.youtube.com/watch?v=eYCYGpu0OxM',
            # Along For The Ride
            'https://www.youtube.com/watch?v=bp85GgcyESs',
            # Another Day
            'https://www.youtube.com/watch?v=LYtiDCXLAcQ'
        ]

        webbrowser.open_new_tab(selection[random.randint(0, len(selection))])


    def plottableSelected(self, newIndex, oldIndex=None):
        """
        Event triggered when a plottable is selected in the TreeView.
        """
        item = self.treeViewModel.itemFromIndex(newIndex.indexes()[0])

        try:
            if item.data() == 'unknown':
                info = self.data.getUnknownInfo(item.text())
            elif item.data() == 'other':
                path = self.getPlottablePath(item)[1:]
                fullname = path.partition('/')[-1]
                info = self.data.getOtherInfo(fullname)
            else:
                # Ignore
                return
        except RuntimeError:
            QMessageBox.critical(self,
                "Unable to display information about quantity",
                "DREAM Theater is currently unable to display information about the quantity. Try launching the simulation.")
            return

        self.ui.lblDescription.setText(info['description'])
        self.ui.lblNElements.setText('{:d}'.format(info['nelements']))
        self.ui.lblNMultiples.setText('{:d}'.format(info['nmultiples']))

        if 'equation' in info:
            self.ui.lblEquation.setText(info['equation'])
        else:
            self.ui.lblEquation.setText('')


    def plottableContextMenu(self, position):
        """
        Event triggered when a quantity in the TreeView is right-clicked.
        """
        item = self.getSelectedPlottable()
        # Only show context menu for items which don't have children
        # (i.e. only those items which represent actual plottables)
        if item.hasChildren():
            return

        menu = QtWidgets.QMenu()

        # Monitor
        actionMonitor = menu.addAction("Monitor")
        #actionMonitor.setCheckable(True)
        #if self.plotconfig.monitorsQuantity(item.text()):
        #    actionMonitor.setChecked(True)

        # Separator
        menu.addSeparator()

        # Plot
        actionPlot = menu.addAction("Plot")

        action = menu.exec_(self.ui.treeView.viewport().mapToGlobal(position))

        if action == actionPlot:
            QMessageBox.critical(self,
                "Plotting not implemented yet",
                "Plotting of quantities has not been implemented yet")
        elif action == actionMonitor:
            # TODO determine if this is a kinetic, fluid or scalar quantity...
            diag = DialogFluidQuantity(r=self.output.grid.r, parent=self)
            if not diag.exec_():
                return
            
            if diag.getPlotType() == diag.PLOT_TYPE_RADIAL_PROFILE:
                x = 'r'
                y = lambda data : data[-1,:]
            else:
                ridx = diag.getRadialIndex()
                x = 't'
                y = lambda data : data[:,ridx]

            self.plotconfig.addQuantity(item.text(), x=x, y=y, xlabel=x)
            self.plotconfig.render(self.data)

    
    def runScript(self):
        """
        Run the current Python script.
        """
        self.enableEditor(False)
        self.editor.viewport().setCursor(QCursor(Qt.WaitCursor))

        # Check if a subset of the code is selected...
        code = self.editor.textCursor().selection().toPlainText()

        # If no selection was found, execute all code...
        if not code:
            code = self.editor.toPlainText()

        try:
            # TODO obtain output from running Simulation object...
            evaluateExpression.evaluate(code, output=self.output)
        except Exception as ex:
            QMessageBox.critical(self,
                'Exception was raised',
                '{}: The Python code raised the following exception:\n\n{}'.format(type(ex).__name__, traceback.format_exc()))

        self.editor.viewport().setCursor(QCursor(Qt.IBeamCursor))
        self.enableEditor(True)
        self.editor.setFocus()


    def runSimulation(self):
        """
        Run a DREAM simulation.
        """
        # Prevent multiple simulations from running simultaneously
        # in the same window...
        if self.simulationThread != None:
            return

        if not self.isSimulation():
            print('ERROR: No simulation loaded to run.')
            return

        self.enableSimulationControls(False)
        self.ui.pbSimulation.setValue(0)

        self.simulationThread = SimulationThread(self.data.simulation)
        self.simulationThread.timestepTaken.connect(self.timestepTaken)
        self.simulationThread.finished.connect(self.simulationFinished)
        self.simulationThread.start()


    def saveScript(self):
        """
        Save the current Python script.
        """
        filename, _ = QFileDialog.getSaveFileName(self, "Save as...", self.scriptFilename, "Python script (*.py);;All files (*.*)")

        if filename:
            with open(filename, 'w', encoding='utf8') as f:
                f.write(self.editor.toPlainText())

            self.setScriptModified(False)


    def scriptHasBeenModified(self):
        """
        Returns ``TRUE`` if the current script was modified.
        """
        return self.editor.document().isModified()


    def onScriptModified(self, mod):
        """
        Event triggered when the Python script has been modified.
        """
        self.ui.btnSaveScript.setEnabled(mod)


    def setOutput(self, output):
        """
        Set the DREAMOutput object to analyse.
        """
        self.output = output


    def setScriptModified(self, mod):
        """
        Set that the script has been modified (or not).
        """
        self.editor.document().setModified(mod)


    def setupPlots(self):
        """
        Initialization of the 'Plots' section on the window.
        """
        self.figure = Figure(tight_layout=True)
        self.canvas = FigureCanvas(self.figure)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        self.ui.framePlot.setLayout(layout)

        self.plotconfig = PlotConfiguration(self.figure)
        # Bind mouse click event
        self.plotconfig.fig.canvas.mpl_connect('button_press_event', self.onFigureClicked)


    def simulationFinished(self):
        """
        Method called when the simulation finishes.
        """
        self.setOutput(self.simulationThread.getOutput())
        self.simulationThread = None
        self.ui.pbSimulation.setValue(100)

        self.enableSimulationControls(True)


    def timestepTaken(self, simulation):
        """
        Method called whenever the currently running simulation has
        completed another time step.
        """
        maxt = simulation.getMaxTime()
        curt = simulation.getCurrentTime()

        self.ui.pbSimulation.setValue(int(curt/maxt*100))


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)

    win = DREAMTheater('../../examples/runaway/dream_settings.h5')
    win.show()
    sys.exit(app.exec_())


