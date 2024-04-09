
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox
from SettingsDiff import SettingsDiff

from ui import MainWindow_design


class MainWindow(QtWidgets.QMainWindow):
    

    def __init__(self, file1, file2):
        """
        Constructor.
        """
        QtWidgets.QMainWindow.__init__(self)

        self.ui = MainWindow_design.Ui_MainWindow()
        self.ui.setupUi(self)

        self.filterBoxes = []
        self.makeDiff(file1, file2)

        self.bindEvents()


    def bindEvents(self):
        """
        Bind control events to methods.
        """
        self.ui.listWidget.currentItemChanged.connect(self.diffSelected)


    def diffSelected(self, newItem, prevItem):
        """
        An item in the diff list is selected.
        """
        if newItem is None:
            self.ui.tbValue1.setText('')
            self.ui.tbValue2.setText('')

            self.ui.lblReason.setText('n/a')
            self.ui.lblType1.setText('n/a')
            self.ui.lblType2.setText('n/a')
            self.ui.lblShape1.setText('n/a')
            self.ui.lblShape2.setText('n/a')
        else:
            diff = newItem.data(QtCore.Qt.UserRole)
            self.ui.tbValue1.setPlainText(diff.val1)
            self.ui.tbValue2.setPlainText(diff.val2)

            self.ui.lblReason.setText(diff.reason)

            if diff.type1 is not None:
                self.ui.lblType1.setText(diff.type1)
            else: self.ui.lblType1.setText('n/a')
            if diff.type2 is not None:
                self.ui.lblType2.setText(diff.type2)
            else: self.ui.lblType2.setText('n/a')

            if diff.shape1 is not None:
                self.ui.lblShape1.setText(str(diff.shape1))
            else: self.ui.lblShape1.setText('n/a')
            if diff.shape2 is not None:
                self.ui.lblShape2.setText(str(diff.shape2))
            else: self.ui.lblShape2.setText('n/a')


    def makeDiff(self, file1, file2):
        """
        Make a diff between the two files.
        """
        sd = SettingsDiff()
        difflist = sd.diff(file1, file2)

        self.ui.lblFile1.setText(sd.file1)
        self.ui.lblFile2.setText(sd.file2)

        for cb in self.filterBoxes:
            cb.deleteLater()
        self.filterBoxes = []

        reasons = []
        self.ui.listWidget.clear()
        for itm in difflist:
            ql = QtWidgets.QListWidgetItem(itm.path)
            ql.setData(QtCore.Qt.UserRole, itm)
            if 'Non-existant' in itm.reason:
                ql.setBackground(QtCore.Qt.lightGray)

            self.ui.listWidget.addItem(ql)
            
            if itm.reason not in reasons:
                reasons.append(itm.reason)

        for r in reasons:
            cb = QtWidgets.QCheckBox(r)
            cb.setChecked(True)
            cb.stateChanged.connect(self.filterChanged)
            self.filterBoxes.append(cb)
            self.ui.filterLayout.addWidget(cb)

        self.ui.lblStatus.setText(f'{len(difflist)} differences found ({len(difflist)} shown).')


    def filterChanged(self):
        """
        A filter checkbox was checked/unchecked.
        """
        reasons = []
        for cb in self.filterBoxes:
            if cb.isChecked():
                reasons.append(cb.text())

        items = 0
        for row in range(self.ui.listWidget.count()):
            itm = self.ui.listWidget.item(row)
            dif = itm.data(QtCore.Qt.UserRole)
            if dif.reason not in reasons:
                #self.ui.listWidget.setItemHidden(itm, True)
                itm.setHidden(True)
            else:
                #self.ui.listWidget.setItemHidden(itm, False)
                itm.setHidden(False)
                items += 1

        self.ui.lblStatus.setText(f'{self.ui.listWidget.count()} differences found ({items} shown).')


