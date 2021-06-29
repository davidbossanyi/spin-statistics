from PyQt5 import QtWidgets
import sys
from spinstats_gui import Application

QtWidgets.QApplication.setStyle('Fusion')
app = QtWidgets.QApplication(sys.argv)
window = Application()
window.show()
sys.exit(app.exec_())