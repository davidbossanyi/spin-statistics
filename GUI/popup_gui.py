# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'popup.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_popup(object):
    def setupUi(self, popup):
        popup.setObjectName("popup")
        popup.resize(399, 388)
        self.picture = QtWidgets.QLabel(popup)
        self.picture.setGeometry(QtCore.QRect(10, 10, 379, 365))
        self.picture.setText("")
        self.picture.setObjectName("picture")

        self.retranslateUi(popup)
        QtCore.QMetaObject.connectSlotsByName(popup)

    def retranslateUi(self, popup):
        _translate = QtCore.QCoreApplication.translate
        popup.setWindowTitle(_translate("popup", "spinstats GUI"))

