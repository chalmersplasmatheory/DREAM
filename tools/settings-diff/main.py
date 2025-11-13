#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
p = str((Path(__file__).parent / '..').resolve().absolute())
sys.path.append(p)


def parse_args():
    parser = argparse.ArgumentParser("Compare two DREAM settings files")

    parser.add_argument("file1", help="Name of the first file to compare.")
    parser.add_argument("file2", help="Name of the second file to compare.")
    parser.add_argument("-c", "--cli", help="Run without the GUI", action="store_true")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    if args.cli:
        pass
    else:
        from PyQt5 import QtWidgets
        from MainWindow import MainWindow

        app = QtWidgets.QApplication(sys.argv)
        window = MainWindow(args.file1, args.file2)
        window.show()
        sys.exit(app.exec_())


