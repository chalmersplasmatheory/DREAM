
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QFileDialog, QMessageBox

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from DREAM import DREAMIO


class PlotShapingWindow(QtWidgets.QFrame):
    def __init__(self, params, width=900, height=500, parent=None):
        super(PlotShapingWindow, self).__init__(parent)

        self.figure = Figure(tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.subplots(2,3)
        self.setWindowTitle('Shaping parameters')

        self.btnSave = QtWidgets.QPushButton(self)
        self.btnSave.setObjectName("btnSave")
        self.btnSave.setText("Save parameters")
        self.btnSave.clicked.connect(self.saveParameters)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.btnSave)
        self.setLayout(layout)
        self.resize(width,height)

        self.parameters = params

        self.plotShaping(params)


    def drawSafe(self):
        try:
            self.canvas.draw()
        except RuntimeError as e:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText(e.strerror)
            msg.setWindowTitle('Runtime Error')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()

    
    def plotShaping(self, params):
        self.plotParameter(self.ax[0,0], params['r'], params['psi'], r'$\psi$', 'Poloidal flux')
        self.plotParameter(self.ax[0,1], params['r'], params['kappa'], r'$\kappa$', 'Elongation')
        self.plotParameter(self.ax[0,2], params['r'], params['delta'], r'$\delta$', 'Triangularity')
        self.plotParameter(self.ax[1,0], params['r'], params['Delta'], r'$\Delta$', 'Shafranov shift')
        self.plotParameter(self.ax[1,1], params['r'], params['GOverR0'], r'$G/R_0$', 'Toroidal magnetic field function')
        self.figure.delaxes(self.ax[1,2])

        self.figure.tight_layout()


    def plotParameter(self, ax, r, p, name, title):
        """
        Plot a single shaping parameter on the specified axes.
        """
        ax.plot(r, p, 'r')
        ax.set_xlim([0, max(r)])

        if min(p) >= 0:
            ax.set_ylim([0, 1.1*max(p)])

        ax.set_xlabel(r'$r$')
        ax.set_ylabel(name)
        ax.set_title(title)


    def saveParameters(self):
        filename, _ = QFileDialog.getSaveFileName(self, caption="Save shaping parameters to HDF5", filter="HDF5 file (*.h5)")

        if filename:
            DREAMIO.SaveDictAsHDF5(filename, self.parameters)
            QMessageBox.information(self, "Parameters successfully saved", f"The equilibrium shaping parameters were successfully saved to '{filename}'.")


