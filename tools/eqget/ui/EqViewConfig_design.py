# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui/EqViewConfig.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DialogEqConfig(object):
    def setupUi(self, DialogEqConfig):
        DialogEqConfig.setObjectName("DialogEqConfig")
        DialogEqConfig.resize(400, 220)
        self.verticalLayout = QtWidgets.QVBoxLayout(DialogEqConfig)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(DialogEqConfig)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.tbCocos = QtWidgets.QLineEdit(DialogEqConfig)
        self.tbCocos.setObjectName("tbCocos")
        self.verticalLayout.addWidget(self.tbCocos)
        self.cbOverride = QtWidgets.QCheckBox(DialogEqConfig)
        self.cbOverride.setObjectName("cbOverride")
        self.verticalLayout.addWidget(self.cbOverride)
        self.tbOverride = QtWidgets.QLineEdit(DialogEqConfig)
        self.tbOverride.setObjectName("tbOverride")
        self.verticalLayout.addWidget(self.tbOverride)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.buttonBox = QtWidgets.QDialogButtonBox(DialogEqConfig)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(DialogEqConfig)
        self.buttonBox.accepted.connect(DialogEqConfig.accept) # type: ignore
        self.buttonBox.rejected.connect(DialogEqConfig.reject) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(DialogEqConfig)

    def retranslateUi(self, DialogEqConfig):
        _translate = QtCore.QCoreApplication.translate
        DialogEqConfig.setWindowTitle(_translate("DialogEqConfig", "Configure equilibrium"))
        self.label.setText(_translate("DialogEqConfig", "<html><head/><body><p>COCOS (<a href=\"https://crppwww.epfl.ch/~sauter/cocos/\"><span style=\" text-decoration: underline; color:#0000ff;\">documentation</span></a>):</p></body></html>"))
        self.cbOverride.setText(_translate("DialogEqConfig", "Override ψ at LCFS value:"))
