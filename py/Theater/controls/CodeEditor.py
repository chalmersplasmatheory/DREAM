# An extension of the Qt5 QPlainTextEdit control to make
# editing code more convenient.
#
# Obtained from https://gist.github.com/Axel-Erfurt/8c84b5e70a1faf894879cd2ab99118c2.
# Modified by Mathias Hoppe.

from PyQt5.QtWidgets import QPlainTextEdit, QCompleter
from PyQt5.QtCore import Qt, QStringListModel, pyqtSignal
from PyQt5.QtGui import QTextCursor
from . PythonHighlighter import PythonHighlighter


class CodeEditor(QPlainTextEdit):

    
    executeCommand = pyqtSignal()


    def __init__(self, parent=None, window=None):
        super().__init__(parent)

        self.installEventFilter(self)
        self._completer = None
        self._highlighter = PythonHighlighter(self.document())

        comp = QCompleter(window)
        comp.setModel(getCompletionWordlist(comp))
        comp.setModelSorting(QCompleter.CaseInsensitivelySortedModel)
        comp.setCaseSensitivity(Qt.CaseInsensitive)
        comp.setWrapAround(True)
        comp.setCompletionRole(Qt.EditRole)
        self.setCompleter(comp)

        self.setLineWrapMode(QPlainTextEdit.NoWrap)

        self.setTabStopWidth(40)
        self.setStyleSheet(editorStyle(self))


    def setCompleter(self, c):
        if self._completer is not None:
            self._completer.activated.disconnect()

        self._completer = c
#        c.popup().verticalScrollBar().hide()
        c.popup().setStyleSheet("background-color: #555753; color: #eeeeec; font-size: 8pt; selection-background-color: #4e9a06;")

        c.setWidget(self)
        c.setCompletionMode(QCompleter.PopupCompletion)
        c.activated.connect(self.insertCompletion)


    def completer(self):
        return self._completer


    def insertCompletion(self, completion):
        if self._completer.widget() is not self:
            return

        tc = self.textCursor()
        extra = len(completion) - len(self._completer.completionPrefix())
        tc.movePosition(QTextCursor.Left)
        tc.movePosition(QTextCursor.EndOfWord)
        tc.insertText(completion[-extra:])
        self.setTextCursor(tc)


    def textUnderCursor(self):
        tc = self.textCursor()
        tc.select(QTextCursor.WordUnderCursor)

        return tc.selectedText()


    def focusInEvent(self, e):
        if self._completer is not None:
            self._completer.setWidget(self)

        super().focusInEvent(e)


    def keyPressEvent(self, e):
        if (e.modifiers() == Qt.ControlModifier) and (e.key() in (Qt.Key_Enter, Qt.Key_Return)):
            self.executeCommand.emit()
            return
        if e.key() == Qt.Key_Tab:
            self.textCursor().insertText("    ")
            return
        if self._completer is not None and self._completer.popup().isVisible():
            # The following keys are forwarded by the completer to the widget.
            if e.key() in (Qt.Key_Enter, Qt.Key_Return):
                e.ignore()
                # Let the completer do default behavior.
                return

        isShortcut = ((e.modifiers() & Qt.ControlModifier) != 0 and e.key() == Qt.Key_Escape)
        if self._completer is None or not isShortcut:
            # Do not process the shortcut when we have a completer.
            super().keyPressEvent(e)

        ctrlOrShift = e.modifiers() & (Qt.ControlModifier | Qt.ShiftModifier)
        if self._completer is None or (ctrlOrShift and len(e.text()) == 0):
            return

        eow = "~!@#$%^&*()_+{}|:\"<>?,./;'[]\\-="
        hasModifier = (e.modifiers() != Qt.NoModifier) and not ctrlOrShift
        completionPrefix = self.textUnderCursor()

        if not isShortcut and (hasModifier or len(e.text()) == 0 or len(completionPrefix) < 2 or e.text()[-1] in eow):
            self._completer.popup().hide()
            return

        if completionPrefix != self._completer.completionPrefix():
            self._completer.setCompletionPrefix(completionPrefix)
            self._completer.popup().setCurrentIndex(
                    self._completer.completionModel().index(0, 0))

        cr = self.cursorRect()
        cr.setWidth(self._completer.popup().sizeHintForColumn(0) + self._completer.popup().verticalScrollBar().sizeHint().width())
        self._completer.complete(cr)


def editorStyle(self):
    return """
QPlainTextEdit {
    font-family: Courier New, monospace;
    font-size: 16px;
    background: #333333;
    color: #EEE;
    border: 1px solid #1EAE3D;
}
"""

def getCompletionWordlist(completer):
    """
    Returns the list of words to use for auto-completion.
    """
    return QStringListModel(['output', 'DREAMOutput', 'DREAMSettings'], completer)
