
from PyQt5 import QtCore, QtWidgets
from ui import EqViewConfig_design


class EqViewConfig(QtWidgets.QDialog):
    

    def __init__(self, parent=None):
        """
        Constructor.
        """
        super().__init__(parent=parent)

        self.ui = EqViewConfig_design.Ui_DialogEqConfig()
        self.ui.setupUi(self)

        self.toggleEnabled()
        self.bindEvents()


    def bindEvents(self):
        self.ui.cbOverride.toggled.connect(self.toggleEnabled)


    def toggleEnabled(self):
        """
        Toggle enabled state of relevant controls.
        """
        enbl = self.ui.cbOverride.isChecked()
        self.ui.tbOverride.setEnabled(enbl)


    def assemble(self):
        """
        Assemble the results of the configuration.
        """
        params = {}
        cocos = self.ui.tbCocos.text()
        try:
            cocos = int(cocos)
            if cocos < 1 or cocos > 18 or cocos in [9, 10]:
                raise ValueError("Invalid COCOS number.")

            params['cocos'] = cocos
        except ValueError as ex:
            QMessageBox.critical(self, "Invalid COCOS number", "An invalid COCOS number was specified. Must be 1-8 or 11-18.")

        if self.ui.cbOverride.isChecked():
            val = self.ui.tbOverride.text()
            if val:
                try:
                    params['override_psilim'] = float(val)
                except ValueError as ex:
                    params['override_psilim'] = True
            else:
                params['override_psilim'] = True

        return params


    @staticmethod
    def configureEquil(params):
        """
        Open the dialog to configure equilibrium parameters.
        """
        d = EqViewConfig()
        if 'cocos' in params:
            d.ui.tbCocos.setText(str(params['cocos']))
        if 'override_psilim' in params:
            if params['override_psilim'] != False:
                d.ui.cbOverride.setChecked(True)
                if type(params['override_psilim']) != bool:
                    d.ui.tbOverride.setText(str(params['override_psilim']))
        else:
            d.ui.cbOverride.setChecked(False)

        d.toggleEnabled()

        if d.exec():
            return d.assemble()
        else:
            return None

        
