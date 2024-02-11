#!/usr/bin/env python3

import h5py
import numpy as np
import sys
import traceback

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from ui import EqView_design
from PlotShapingWindow import PlotShapingWindow
from DREAM import DREAMIO
from DREAM.Settings.LUKEMagneticField import LUKEMagneticField

import EqFile
from EqViewConfig import EqViewConfig
from EQDSK import EQDSK

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


EQTYPE_NONE = 0
EQTYPE_LUKE = 1
EQTYPE_EQDSK = 2


class EqView(QtWidgets.QMainWindow):
    

    def __init__(self, argv):
        """
        Constructor.
        """
        QtWidgets.QMainWindow.__init__(self)

        self.ui = EqView_design.Ui_EqView()
        self.ui.setupUi(self)

        self.equil = None
        self.equil_type = 0
        self.frameAx = None

        self.canvas = FigureCanvas(Figure())
        self.frameToolbar = NavigationToolbar(self.canvas, self)

        self.frameLayout = QtWidgets.QVBoxLayout(self.ui.frame)
        self.frameLayout.addWidget(self.frameToolbar)
        self.frameLayout.addWidget(self.canvas)

        self.windows = {}

        self.bindEvents()


    def bindEvents(self):
        """
        Bind control events to methods.
        """
        self.ui.actionExit.triggered.connect(self.exit)
        self.ui.actionOpen.triggered.connect(self.open)
        self.ui.actionSaveAs.triggered.connect(self.save)
        self.ui.actionShapingProfiles.triggered.connect(self.exportShaping)


    def closeEvent(self, event):
        self.exit()


    def exit(self):
        """
        Close any child windows before exiting.
        """
        for _, w in self.windows.items():
            w.close()

        self.close()


    def open(self):
        """
        Open an equilibrium file.
        """
        filename, _ = QFileDialog.getOpenFileName(self, caption="Open equilibrium file", filter="All supported equilibria (*.eqdsk *.geq *.geqdsk *.h5 *.mat);;LUKE Equilibrium (*.h5 *.mat);;EQDSK file (*.eqdsk *.geq *.geqdsk);;All files (*.*)")

        if not filename:
            return

        retry = self.load(filename)
        oldparams = dict(cocos=1, override_psilim=False)
        while retry:
            retry, oldparams = self.special_load(filename, oldparams=oldparams)


    def load(self, filename, params={}):
        """
        Load the named equilibrium file.
        """
        if h5py.is_hdf5(filename):
            try:
                self.load_LUKE(filename)
            except Exception as ex:
                QMessageBox.critical(self, 'Error loading equilibrium', f'An error was encountered when loading the equilibrium.\n\n{ex}')
                return False
        else:
            try:
                self.load_EQDSK(filename, **params)
            except Exception as ex:
                rst = QMessageBox.critical(self, 'Error loading equilibrium', f'An error was encountered when loading the equilibrium.\n\n"{ex}"\n\nWould you like to try again with adjusted configuration parameters?', QMessageBox.Yes | QMessageBox.No)
                
                if rst == QMessageBox.Yes:
                    return True
                else: return False

        self.setupPlot()
        self.plot()

        return False


    def load_LUKE(self, filename):
        """
        Load a LUKE equilibrium.
        """
        self.equil = LUKEMagneticField(filename)
        self.equil_type = EQTYPE_LUKE

        self.ui.lblMagneticAxis.setText(f'({self.equil.Rp[0]:.3f}, {self.equil.Zp[0]:.3f})')
        self.ui.lblMinorRadius.setText(f'{self.equil.getMinorRadius():.3f} m')
        self.ui.lblB0.setText(f'{np.mean(self.equil.ptBPHI[:,0]):.3f} T')

        self.ui.actionSaveAs.setEnabled(False)
        self.ui.actionShapingProfiles.setEnabled(False)


    def load_EQDSK(self, filename, **params):
        """
        Load an EQDSK equilibrium.
        """
        self.equil = EQDSK(filename, plot_on_error=False, **params)
        self.equil_params = params
        self.equil_type = EQTYPE_EQDSK

        self.ui.lblMagneticAxis.setText(f'({self.equil.opoint[0]:.3f}, {self.equil.opoint[1]:.3f})')
        self.ui.lblMinorRadius.setText(f'{np.amax(self.equil.lcfs_R)-self.equil.opoint[0]:.3f} m')
        self.ui.lblB0.setText(f'{self.equil.f_psi(0)/self.equil.R0:.3f} T')

        self.ui.actionSaveAs.setEnabled(True)
        self.ui.actionShapingProfiles.setEnabled(True)


    def special_load(self, filename, oldparams):
        """
        Load equilibrium file with a special configuration.
        """
        params = EqViewConfig.configureEquil(oldparams)

        if params is None:
            return False, params

        return self.load(filename, params), params


    def setupPlot(self):
        """
        Set up the main plot.
        """
        if self.frameAx is not None:
            self.canvas.figure.clear()

        self.frameAx = self.canvas.figure.subplots()
        #self.frameAx.set_xlabel(r'$R$ (m)')
        #self.frameAx.set_ylabel(r'$Z$ (m)')
        #self.frameAx.figure.tight_layout()


    def plot(self):
        """
        Plot the equilibrium.
        """
        if self.equil_type == EQTYPE_LUKE:
            self.equil.visualize(ax=self.frameAx)
        elif self.equil_type == EQTYPE_EQDSK:
            self.equil.plot_psi(ax=self.frameAx)
            self.frameAx.figure.tight_layout()

        self.canvas.draw()


    def exportShaping(self):
        """
        Export shaping profiles from the current equilibrium.
        """
        if self.equil_type != EQTYPE_EQDSK:
            QMessageBox.critical(self, "Cannot export equilibrium parameters", "Shaping parameters can only be exported from EQDSK equilibria.")
            return
            
        filename, _ = QFileDialog.getSaveFileName(
            parent=self, caption="Export shaping parameters",
            filter="HDF5 file (*.h5);;All files (*.*)"
        )

        if filename:
            try:
                self.equil.save_eq_parameters(filename, nr=100)
                QMessageBox.information(self, "Successfully saved shaping parameters", f"Successfully saved the shaping parameters to '{filename}'.")
            except Exception as ex:
                QMessageBox.critical(self, "Failed to export shaping parameters", f"An error occurred while saving the equilibrium shaping parameters.\n\n{ex}")


    def save(self):
        """
        Save the equilibrium file to a LUKE format file.
        """
        if self.equil_type != EQTYPE_EQDSK:
            QMessageBox.critical(self, "Cannot save equilibrium", "Only EQDSK equilibria can be exported to LUKE HDF5 files.")
            return

        filename, _ = QFileDialog.getSaveFileName(
            parent=self, caption="Save equilibrium file",
            filter="HDF5 file (*.h5);;All files (*.*)"
        )

        if filename:
            try:
                # TODO allow selection of resolution
                self.equil.save_LUKE(filename, npsi=120, ntheta=150)
                QMessageBox.information(self, "Successfully saved equilibrium", f"Successfully saved equilibrium file to '{filename}'.")
            except Exception as ex:
                QMessageBox.critical(self, "Failed to save equilibrium", f"An error occurred while saving the equilibrium file.\n\n{ex}")


def main():
    app = QtWidgets.QApplication(sys.argv)

    win = EqView(sys.argv[1:])
    win.show()
    return app.exec_()


if __name__ == '__main__':
    sys.exit(main())


