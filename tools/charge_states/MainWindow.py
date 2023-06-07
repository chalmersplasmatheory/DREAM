
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

from ADAS import data, rates
import ionrate
from ui import MainWindow_design


class MainWindow(QtWidgets.QMainWindow):
    

    def __init__(self):
        """
        Constructor.
        """
        QtWidgets.QMainWindow.__init__(self)

        self.ui = MainWindow_design.Ui_MainWindow()
        self.ui.setupUi(self)

        data._initADAS()
        self.elements = data.ELEMENTS

        for e in self.elements.keys():
            self.ui.cbElements.addItem(e)

        self.figure = Figure(tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)

        self.canvasLayout = QtWidgets.QVBoxLayout(self.ui.widget)
        self.canvasLayout.addWidget(self.canvas)

        self.bindEvents()
        self.recalculate()


    def bindEvents(self):
        """
        Bind control events.
        """
        self.ui.cbElements.currentTextChanged.connect(self.recalculate)
        #self.ui.tbTe.textChanged.connect(self.recalculate)
        #self.ui.tbNe.textChanged.connect(self.recalculate)
        #self.ui.tbNi.textChanged.connect(self.recalculate)

        self.ui.btnCalc.clicked.connect(self.recalculate)


    def drawSafe(self):
        """
        Redraw canvas.
        """
        try:
            self.canvas.draw()
        except RuntimeError as e:
            pass


    def recalculate(self):
        """
        Recalculate the charge state distribution.
        """
        try:
            Te = float(self.ui.tbTe.text())
            ne = float(self.ui.tbNe.text())
            el = self.ui.cbElements.currentText()

            acd = rates.get_rate(el, 'acd')
            scd = rates.get_rate(el, 'scd')

            n = ionrate.solve_all(ni=1, ne=ne, Te=Te, acd=acd, scd=scd)
            Z = list(range(n.size))
            self.ax.clear()
            hbars = self.ax.bar(Z, n*100)
            lbls = [f'{z*100:.1f}%' if z>.001 else '' for z in n]
            self.ax.bar_label(hbars, labels=lbls)
            self.ax.set_xlabel('$Z$')
            self.ax.set_ylabel('Abundance (\%)')
            self.ax.set_ylim([0, 100])
            self.ax.set_xticks(Z)
            self.drawSafe()
        except ValueError: pass


