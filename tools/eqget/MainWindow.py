
import numpy as np
import traceback

import sys
sys.path.append('../../py')

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from ui import MainWindow_design
from PlotWindow import PlotWindow
from PlotShapingWindow import PlotShapingWindow
from DREAM import DREAMIO
from pathlib import Path

import AUG
import EqFile

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MainWindow(QtWidgets.QMainWindow):
    

    def __init__(self, argv):
        """
        Constructor.
        """
        QtWidgets.QMainWindow.__init__(self)

        self.ui = MainWindow_design.Ui_EqGet()
        self.ui.setupUi(self)

        self.equil = None

        # Set up flux surface figure
        self.canvas = FigureCanvas(Figure())
        self.fluxSurfaceLayout = QtWidgets.QVBoxLayout(self.ui.frameFluxSurfaces)
        self.fluxSurfaceLayout.addWidget(self.canvas)
        self.setupPlot()

        # List of open windows
        self.windows = {}

        if AUG.isAvailable():
            self.ui.cbTokamak.addItem('ASDEX Upgrade', AUG)

        self.ui.cbTokamak.addItem('File', EqFile)
        self.tokamakChanged()

        self.toggleEnabled(False)
        self.bindEvents()

        if len(argv) == 1:
            self.parsearg(argv[0])


    def bindEvents(self):
        """
        Bind control events to methods.
        """
        self.ui.actionExit.triggered.connect(self.exit)
        self.ui.btnLoad.clicked.connect(self.load)

        self.ui.btnPlotPsi.clicked.connect(self.plotPsi)
        self.ui.btnPlotB.clicked.connect(self.plotB)
        self.ui.btnPlotBpol.clicked.connect(self.plotBpol)
        self.ui.btnPlotBr.clicked.connect(self.plotBr)
        self.ui.btnPlotBz.clicked.connect(self.plotBz)
        self.ui.btnPlotBphi.clicked.connect(self.plotBphi)

        self.ui.btnShaping.clicked.connect(self.calculateShaping)
        self.ui.btnSave.clicked.connect(self.save)

        self.ui.cbTokamak.currentTextChanged.connect(self.tokamakChanged)


    def closeEvent(self, event):
        self.exit()


    def exit(self):
        """
        Close any child windows before exiting.
        """
        for _, w in self.windows.items():
            w.close()

        self.close()


    def toggleEnabled(self, enabled=True):
        """
        Toggle the enabled state of controls which require data to be
        available.
        """
        self.ui.btnPlotPsi.setEnabled(enabled)
        self.ui.btnPlotB.setEnabled(enabled)
        self.ui.btnPlotBpol.setEnabled(enabled)
        self.ui.btnPlotBr.setEnabled(enabled)
        self.ui.btnPlotBz.setEnabled(enabled)
        self.ui.btnPlotBphi.setEnabled(enabled)
        self.ui.btnShaping.setEnabled(enabled)
        self.ui.btnSave.setEnabled(enabled)


    def tokamakChanged(self, e=None):
        """
        Event fired when the selected tokamak handler changes.
        """
        if self.ui.cbTokamak.currentText() == "File":
            self.ui.lblDischarge.setText("Equilibrium file")
            self.ui.btnLoad.setText("Open...")
            
            self.ui.lblTime.setText('Psi override')
            self.ui.tbTime.setText('')
        else:
            self.ui.lblDischarge.setText("Discharge")
            self.ui.btnLoad.setText("Load")

            self.ui.lblTime.setText('Time')
            self.ui.tbTime.setText('1.0')


    def getShot(self):
        shot = self.ui.tbShot.text()

        # Try to convert to integer. If that fails, the user may
        # have provided a file name instead...
        try: shot = int(shot)
        except: pass

        return shot

    
    def load(self):
        """
        Load data using the selected module.
        """
        shot = self.getShot()

        if self.ui.cbTokamak.currentText() == "File" and not shot:
            shot, _ = QFileDialog.getOpenFileName(self, caption="Open equilibrium file", filter="All supported equilibria (*.eqdsk *.geqdsk *.h5 *.mat);;LUKE equilibrium (*.h5 *.mat);;EQDSK file (*.eqdsk *.geqdsk);;All files (*.*)")
            if not shot:
                return

            self.ui.tbShot.setText(shot)

        self._load_internal(shot)


    def _load_internal(self, data):
        try:
            mod = self.ui.cbTokamak.currentData()
            time = self.ui.tbTime.text()
            if not time: time = False
            else: time = float(time)

            self.equil = mod.getLUKE(data, time)
            print("Loaded '{}'...".format(data))

            self.plotFluxSurfaces()
            self.toggleEnabled(True)
        except Exception as ex:
            QMessageBox.critical(self, 'Error loading data', f"The specified data file could not be loaded:\n\n{ex}")


    def calculateShaping(self):
        """
        Calculate shaping parameters for the loaded magnetic equilibrium.
        """
        pass
        try:
            mod = self.ui.cbTokamak.currentData()
            if not hasattr(mod, 'getShaping'):
                raise Exception("The selected equilibrium handler does not support calculating shaping parameters.")

            params = mod.getShaping(self.getShot(), equil=self.equil)
            self.plotShaping(params)
        except Exception as ex:
            QMessageBox.critical(self, 'Error loading shot', f"The specified shot file could not be loaded:\n\n{ex}\n\n{traceback.format_exc()}")


    def parsearg(self, arg):
        """
        Parse an input argument.
        """
        if Path(arg).is_file():
            self.ui.cbTokamak.setCurrentText('File')
            self._load_internal(arg)


    def plotFluxSurfaces(self):
        """
        Plot flux surfaces from loaded equilibrium data.
        """
        ax = self.fluxSurfaceAx
        ptx = self.equil['ptx']
        pty = self.equil['pty']
        Rp  = self.equil['Rp']
        Zp  = self.equil['Zp']

        ax.plot(ptx[:,:-1]+Rp, pty[:,:-1]+Zp, linewidth=0.7, color=(0.5, 0.5, 0.5))
        ax.plot(ptx[:,-1]+Rp, pty[:,-1]+Zp, linewidth=2, color='r')
        ax.plot(Rp, Zp, 's', color='r')
        ax.axis('equal')

        self.canvas.draw()


    def plotPsi(self):
        """
        Plot poloidal flux as function of minor radius.
        """
        if 'psi' in self.windows:
            self.windows['psi'].close()

        w = PlotWindow(600, 400)

        r        = self.equil['ptx'][0,:]
        psi_apRp = self.equil['psi_apRp']
        Rp       = self.equil['Rp']
        ap       = r[-1]

        psi = psi_apRp * (Rp/ap)

        w.ax.plot(r, psi)
        w.ax.set_xlim([0, ap])
        w.ax.set_xlabel(r'$r$ (m)')
        w.ax.set_ylabel(r'Poloidal flux $\Psi$ (Wb)')

        w.show()

        self.windows['psi'] = w


    def plot2D(self, name, data):
        """
        Plot the given magnetic field.
        """
        if name in self.windows:
            self.windows[name].close()

        w = PlotWindow()

        Rp = self.equil['Rp']
        Zp = self.equil['Zp']
        R  = self.equil['ptx'] + Rp
        Z  = self.equil['pty'] + Zp

        cnt = w.ax.contourf(R, Z, data, cmap='GeriMap', levels=40)
        cbar = w.figure.colorbar(cnt)
        w.ax.set_xlabel('$R$ (m)')
        w.ax.set_ylabel('$Z$ (m)')
        w.ax.axis('equal')

        cbar.set_label('{} (T)'.format(name))

        w.show()

        self.windows[name] = w


    def plotB(self):
        """
        Plot the magnetic field strength in (R, Z).
        """
        Br = self.equil['ptBx']
        Bz = self.equil['ptBy']
        Bp = self.equil['ptBPHI']

        self.plot2D('$|B|$', np.sqrt(Br**2 + Bz**2 + Bp**2))


    def plotBpol(self):
        """
        Plot the poloidal magnetic field.
        """
        Br = self.equil['ptBx']
        Bz = self.equil['ptBy']

        self.plot2D(r'$B_{\rm pol}$', np.sqrt(Br**2 + Bz**2))


    def plotBr(self):
        """
        Plot the radial magnetic field component.
        """
        self.plot2D(r'$B_r$', self.equil['ptBx'])


    def plotBz(self):
        """
        Plot the radial magnetic field component.
        """
        self.plot2D(r'$B_z$', self.equil['ptBy'])


    def plotBphi(self):
        """
        Plot the toroidal magnetic field component.
        """
        self.plot2D(r'$B_\varphi$', self.equil['ptBPHI'])


    def plotShaping(self, params):
        """
        Plot the DREAM shaping parameters.
        """
        if 'shaping' in self.windows:
            self.windows['shaping'].close()

        w = PlotShapingWindow(params)
        w.show()
        self.windows['shaping'] = w


    def save(self):
        """
        Save the loaded equilibrium to file.
        """
        filename, _ = QFileDialog.getSaveFileName(self, caption="Save LUKE equilibrium file", filter='HDF5 file (*.h5)')

        if filename:
            DREAMIO.SaveDictAsHDF5(filename, {'equil': self.equil})
            QMessageBox.information(self, "Equilibrium file saved", "The magnetic equilibrium data was saved to the file '{}'.".format(filename))


    def setupPlot(self):
        self.fluxSurfaceAx = self.canvas.figure.subplots()
        self.fluxSurfaceAx.set_xlabel(r'$R$ (m)')
        self.fluxSurfaceAx.set_ylabel(r'$Z$ (m)')
        self.fluxSurfaceAx.figure.tight_layout()


