#!/usr/bin/env python3

from PyQt5 import QtWidgets
import pathlib
import sys

from MainWindow import MainWindow

sys.path.append((pathlib.Path(__file__).parent.parent / "py").resolve())


def main(argv):
    app = QtWidgets.QApplication(sys.argv)

    win = MainWindow(argv)
    win.show()
    return app.exec_()


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


