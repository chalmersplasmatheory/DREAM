#!/usr/bin/env python3

from PyQt5 import QtWidgets
import Theater
import sys


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)

    win = Theater.DREAMTheater('dream_settings.h5')
    #win = Theater.DREAMTheater('../examples/theater/dream_settings.h5')
    #win = Theater.DREAMTheater('../examples/runaway/output.h5')
    win.show()
    sys.exit(app.exec_())


