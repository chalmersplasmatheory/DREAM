#!/usr/bin/env python3

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox
import ManualFit_design

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import pathlib
import sys

cwd = pathlib.Path(__file__).resolve()
ROOT = str(cwd.parent.parent.parent.resolve())
sys.path.append(ROOT)

import ADAS #import data, fit
ADAS.data._initADAS()


class ManualFit(QtWidgets.QMainWindow):
    

    def __init__(self):
        """
        Constructor.
        """
        QtWidgets.QMainWindow.__init__(self)

        self.ui = ManualFit_design.Ui_ManualFit()
        self.ui.setupUi(self)

        self.element_I, self.element_T = None, None

        self.figure = Figure(tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)

        self.fitLayout = QtWidgets.QVBoxLayout(self.ui.widget)
        self.fitLayout.addWidget(self.canvas)

        self.ui.cbMethod.addItem('Single cross-section')
        self.ui.cbMethod.addItem('Single cross-section, 3-parameter')
        self.ui.cbMethod.addItem('Double cross-sections')
        self.ui.cbMethod.setCurrentIndex(0)

        self.showParams({'C1': None, 'C2': None, 'DI1': None, 'DI2': None, 'betaStar': None, 'beta2': None})

        self.loadElements()

        self.bindEvents()

    
    def bindEvents(self):
        """
        Bind to control events.
        """
        self.ui.cbElements.currentIndexChanged.connect(self.elementSelected)
        self.ui.cbCS.currentIndexChanged.connect(self.chargeStateSelected)

        self.ui.hsTlower.valueChanged.connect(self.TlowerChanged)
        self.ui.hsTupper.valueChanged.connect(self.TupperChanged)

        self.ui.btnFit.clicked.connect(self.doFit)


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


    def loadElements(self):
        """
        Load the element named in the 'Element' text box.
        """
        for e in ADAS.data.ELEMENTS.keys():
            self.ui.cbElements.addItem(e)

        self.ui.cbElements.setCurrentIndex(0)
        self.elementSelected()


    def elementSelected(self):
        """
        A new element has been selected in the element combobox.
        """
        global ROOT

        el = self.ui.cbElements.currentText()
        self.element_I, self.element_Z, _, self.element_T, _ = ADAS.data.getIonizationData(el)

        self.ui.hsTlower.setMaximum(self.element_T.size)
        self.ui.hsTupper.setMaximum(self.element_T.size)

        self.TlowerChanged()
        self.TupperChanged()

        self.updateChargeStates()


    def updateChargeStates(self):
        """
        Update the list of available charge states.
        """
        Z = int(self.element_Z)
        
        self.ui.cbCS.clear()

        for i in range(Z):
            self.ui.cbCS.addItem(str(i))

        self.ui.cbCS.setCurrentIndex(0)
        self.chargeStateSelected()


    def chargeStateSelected(self):
        """
        A new charge state as been selected in the charge state combobox.
        """
        pass


    def getElement(self):
        """
        Returns the name of the currently selected element.
        """
        return self.ui.cbElements.currentText()


    def getTlower(self):
        """
        Returns the currently selected lower temperature cut-off.
        """
        return self.element_T[self.ui.hsTlower.value()]


    def getTupper(self):
        """
        Returns the currently selected upper temperature cut-off.
        """
        return self.element_T[self.element_T.size - self.ui.hsTupper.value() - 1]


    def getZ0(self):
        """
        Returns the currently selected charge state.
        """
        return int(self.ui.cbCS.currentText())


    def TlowerChanged(self):
        """
        Lower temperature bound changed.
        """
        Tlower = self.getTlower()
        self.ui.lblTlower.setText('{:.3f} eV'.format(Tlower))


    def TupperChanged(self):
        """
        Upper temperature bound changed.
        """
        Tupper = self.getTupper()
        self.ui.lblTupper.setText('{:.3f} eV'.format(Tupper))


    def doFit(self):
        """
        Fit the cross section to the selected element/charge state,
        using the specified model on the chosen temperature interval.
        """
        species = self.getElement()
        Z0 = self.getZ0()

        T_lower = self.getTlower()
        T_upper = self.getTupper()

        if T_lower >= T_upper:
            QMessageBox.critical(self,
                "Invalid temperature range selected",
                "The lower temperature cut-off must be strictly less than the upper temperature cut-off.")
            return
        
        idx = self.ui.cbMethod.currentIndex()
        if idx == 0:
            method = 'single'
        elif idx == 1:
            method = 'single3p'
        elif idx == 2:
            method = 'double'
        else:
            QMessageBox.critical(self,
                "Unrecognized fitting method",
                "Unrecognized fitting method selected: '{}'.".format(self.ui.cbMethod.currentText()))

        _, _, params = ADAS.fit.fitKineticIonizationForSpecies(species, Z0=Z0, fittype=method, T_lower=T_lower, T_upper=T_upper)

        self.drawFit(params)
        self.showParams(params)


    def drawFit(self, params):
        """
        Plot the curves resulting from a fit.
        """
        I_fit = ADAS.fit.evaluateAveragedCrossSection(T=self.element_T, **params)

        Z0 = self.getZ0()
        self.ax.clear()

        self.ax.loglog(self.element_T, self.element_I[Z0,:,0], 'k')
        self.ax.loglog(self.element_T, I_fit, 'r--')
        self.ax.loglog(self.element_T, I_fit, 'rx')

        Tmin, Tmax = self.element_T[0], self.element_T[-1]
        ymin, ymax = 1e-25, 1e-10

        Tl, Tu = self.getTlower(), self.getTupper()
        self.ax.loglog([Tl, Tl], [ymin, ymax], 'c--')
        self.ax.loglog([Tu, Tu], [ymin, ymax], 'c--')

        self.ax.legend(['ADAS', 'Fit'])
        self.ax.set_title('{}$^{{{}+}}$'.format(self.getElement(), Z0))

        self.ax.set_xlim([Tmin, Tmax])
        self.ax.set_ylim([ymin, ymax])

        self.drawSafe()


    def showParams(self, params):
        """
        Visualize the resulting fit parameters.
        """
        def show(val, lbl1, lbl2):
            s = val is not None
            if s:
                lbl2.setText('{:.12f}'.format(val))

            lbl1.setVisible(s)
            lbl2.setVisible(s)

        show(params['C1'], self.ui.lblC1l, self.ui.lblC1)
        show(params['C2'], self.ui.lblC2l, self.ui.lblC2)
        show(params['DI1'], self.ui.lblDI1l, self.ui.lblDI1)
        show(params['DI2'], self.ui.lblDI2l, self.ui.lblDI2)
        show(params['betaStar'], self.ui.lblBetaStarl, self.ui.lblBetaStar)
        show(params['beta2'], self.ui.lblBeta2l, self.ui.lblBeta2)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)

    win = ManualFit()
    win.show()
    sys.exit(app.exec_())


