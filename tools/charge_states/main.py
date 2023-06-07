#!/usr/bin/env python3

from PyQt5 import QtGui, QtWidgets
import sys

from pathlib import Path
p = str((Path(__file__).parent / '..').resolve().absolute())
sys.path.append(p)

from MainWindow import MainWindow

app = None


def show_main():
    global app
    window = MainWindow()
    window.show()
    return app.exec_()


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    sys.exit(show_main())


