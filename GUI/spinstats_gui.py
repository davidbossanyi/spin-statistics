# general
import os
import sys
import numpy as np
import pandas as pd
from scipy.constants import e, k

# gui
from PyQt5 import QtGui, QtWidgets
from main_gui import Ui_MainWindow as MainWindow
from popup_gui import Ui_popup as PopUpWindow

# spinstats classes
sys.path.insert(0, os.path.abspath(os.path.join('..')))
from spinstats.models import Model2 as Model
from spinstats import SpinHamiltonian

# graphs
import pyqtgraph as pg
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

if sys.platform == 'win32':
    import ctypes
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID('arbitrary_string')


class Popup(QtWidgets.QWidget):
    def __init__(self):
        super(Popup, self).__init__()
        self.ui = PopUpWindow()
        self.ui.setupUi(self)
        self.setWindowIcon(QtGui.QIcon('icon.png'))
        self.ui.picture.setPixmap(QtGui.QPixmap('icon.png'))


class Application(QtWidgets.QMainWindow):
    def __init__(self):
        super(Application, self).__init__()
        self.ui = MainWindow()
        self.ui.setupUi(self)
        self.setWindowIcon(QtGui.QIcon('icon.png'))
        self.file_directory = os.path.join(os.path.expanduser('~'), 'Documents')
        self.filetypes = 'SPINSTATS (*.spinstats)'
        self.m = Model()
        self.sh = SpinHamiltonian()
        self.xvars = pd.read_csv('xvars.csv', header=0, index_col=0)
        
        self.updates_on = False
        self.initialised = False
        
        self.ui.xVariable.addItems([
            'J',
            'D',
            'E',
            'X',
            'B',
            'theta',
            'phi',
            'alpha',
            'beta',
            'gamma',
            'kS',
            'kD',
            'kTTA',
            'kSF',
            'k_SF',
            'kTS',
            'kRISC',
            ])
        self.ui.yVariable.addItems([
            'PLQY',
            'UCQY',
            'eta',
            ])
        
        self.ui.valueJ.valueChanged.connect(self.update_J)
        self.ui.valueD.valueChanged.connect(self.update_D)
        self.ui.valueE.valueChanged.connect(self.update_E)
        self.ui.valueX.valueChanged.connect(self.update_X)
        self.ui.valueB.valueChanged.connect(self.update_B)
        self.ui.sliderTheta.valueChanged.connect(self.update_theta)
        self.ui.valueTheta.valueChanged.connect(self.update_theta_sb)
        self.ui.sliderPhi.valueChanged.connect(self.update_phi)
        self.ui.valuePhi.valueChanged.connect(self.update_phi_sb)
        self.ui.sliderAlpha.valueChanged.connect(self.update_alpha)
        self.ui.valueAlpha.valueChanged.connect(self.update_alpha_sb)
        self.ui.sliderBeta.valueChanged.connect(self.update_beta)
        self.ui.valueBeta.valueChanged.connect(self.update_beta_sb)
        self.ui.sliderGamma.valueChanged.connect(self.update_gamma)
        self.ui.valueGamma.valueChanged.connect(self.update_gamma_sb)
        self.ui.valueA.valueChanged.connect(self.update_A)
        self.ui.valueGGamma.valueChanged.connect(self.update_GGamma)
        self.ui.valueEvib.valueChanged.connect(self.update_Evib)
        self.ui.valueT.valueChanged.connect(self.update_kT)
        self.ui.valueET1.valueChanged.connect(self.update_ET1)
        self.ui.valueET2.valueChanged.connect(self.update_ET2)
        self.ui.valueKS.valueChanged.connect(self.update_kS)
        self.ui.valueKD.valueChanged.connect(self.update_kD)
        self.ui.valueKTTA.valueChanged.connect(self.update_kTTA)
        self.ui.valueKSF.valueChanged.connect(self.update_kSF)
        self.ui.valueK_SF.valueChanged.connect(self.update_k_SF)
        self.ui.valueKTS.valueChanged.connect(self.update_kTS)
        self.ui.valueKTF.valueChanged.connect(self.update_kTF)
        self.ui.valueKTFfactor.valueChanged.connect(self.update_kTFfactor)
        self.ui.valueKRISC.valueChanged.connect(self.update_kRISC)
        self.ui.valueKIC1.valueChanged.connect(self.update_kIC1)
        self.ui.valueKIC2.valueChanged.connect(self.update_kIC2)
        self.ui.valueKIC21.valueChanged.connect(self.update_kIC21)
        
        self.ui.checkBoxEnergyGapLaw.stateChanged.connect(self.update_use_energygap)
        self.ui.checkBoxKTFfactor.stateChanged.connect(self.update_use_kTFfactor)
        
        self.ui.valueXmin.valueChanged.connect(self.update_xmin)
        self.ui.valueXmax.valueChanged.connect(self.update_xmax)
        self.ui.checkBoxLogx.stateChanged.connect(self.update_logx)
        
        self.ui.saveParamsButton.clicked.connect(self.save_parameters_to_file)
        self.ui.loadParamsButton.clicked.connect(self.load_parameters_from_file)
        self.ui.saveDataButton.clicked.connect(self.save_plot_data)
        
        self.ui.showModelButton.clicked.connect(self.show_rate_model)
        
        self.use_energygap = True
        self.use_kTF_factor = True
        self.ET1 = 1
        self.ET2 = 2
        self.load_parameters('last_instance_parameters.spinstats')
        
        self.ui.checkBoxEnergyGapLaw.toggle()
        self.ui.checkBoxKTFfactor.toggle()
        self.ui.checkBoxEnergyGapLaw.toggle()
        self.ui.checkBoxKTFfactor.toggle()
        
        self.update_use_energygap()
        self.update_use_kTFfactor()
        
        self.ui.graph.plotItem.showAxis('top', show=True)
        self.ui.graph.plotItem.showAxis('right', show=True)
        
        self.update_logx()
        
        self.update_yvar()
        self.ui.xVariable.setCurrentIndex(int(self.params['xvar']))
        self.update_xvar()
        
        self.ui.xVariable.currentTextChanged.connect(self.update_xvar)
        self.ui.yVariable.currentTextChanged.connect(self.update_yvar)
        
        self.initialised = True
        self.updates_on = True
        
        self.display_status('application running', 'green')

    def display_status(self, message, colour, msecs=0):
        self.ui.statusbar.clearMessage()
        self.ui.statusbar.setStyleSheet('QStatusBar{color:'+colour+';}')
        self.ui.statusbar.showMessage(message, msecs=msecs)
        
    def show_rate_model(self):
        self.popup = Popup()
        self.popup.show()
        
    def load_parameters_from_file(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'load parameter file', self.file_directory, self.filetypes)[0]
        if fname != '':
            self.file_directory = os.path.dirname(fname)
            self.load_parameters(fname)
            self.display_status('parameters loaded from file', 'blue', msecs=5000)
            
    def save_parameters_to_file(self):
        fname = QtWidgets.QFileDialog().getSaveFileName(self, 'save parameter file as', self.file_directory, self.filetypes)[0]
        if fname != '':
            self.file_directory = os.path.dirname(fname)
            self.save_parameters(fname)
            self.display_status('parameters saved', 'blue', msecs=5000)
        
    def load_parameters(self, fpath):
        params = pd.read_csv(fpath, header=None, index_col=0, squeeze=True)
        self.ui.valueJ.setValue(params['J'])
        self.ui.valueD.setValue(params['D'])
        self.ui.valueE.setValue(params['E'])
        self.ui.valueX.setValue(params['X'])
        self.ui.valueB.setValue(params['B'])
        self.ui.sliderTheta.setValue(int(params['theta']))
        self.ui.sliderPhi.setValue(int(params['phi']))
        self.ui.sliderAlpha.setValue(int(params['alpha']))
        self.ui.sliderBeta.setValue(int(params['beta']))
        self.ui.sliderGamma.setValue(int(params['gamma']))
        self.ui.valueET1.setValue(params['ET1'])
        self.ui.valueET2.setValue(params['ET2'])
        self.ui.valueA.setValue(params['A'])
        self.ui.valueGGamma.setValue(params['ggamma'])
        self.ui.valueEvib.setValue(params['Evib'])
        self.ui.valueT.setValue(params['T'])
        self.ui.valueKS.setValue(params['kS'])
        self.ui.valueKD.setValue(params['kD'])
        self.ui.valueKTTA.setValue(params['kTTA'])
        self.ui.valueKSF.setValue(params['kSF'])
        self.ui.valueK_SF.setValue(params['k_SF'])
        self.ui.valueKTS.setValue(params['kTS'])
        self.ui.valueKTF.setValue(params['kTF'])
        self.ui.valueKTFfactor.setValue(params['kTFfactor'])
        self.ui.valueKRISC.setValue(params['kRISC'])
        self.ui.valueKIC1.setValue(params['kIC1'])
        self.ui.valueKIC2.setValue(params['kIC2'])
        self.ui.valueKIC21.setValue(params['kIC21'])
        self.ui.checkBoxKTFfactor.setChecked(int(params['use_kTFfactor']))
        self.ui.checkBoxEnergyGapLaw.setChecked(int(params['use_gaplaw']))
        self.ui.valueXmin.setValue(params['xmin'])
        self.ui.valueXmin.setValue(params['xmax'])
        self.ui.xVariable.setCurrentIndex(int(params['xvar']))
        self.ui.yVariable.setCurrentIndex(int(params['yvar']))
        self.ui.checkBoxLogx.setChecked(int(params['logx']))
        self.params = params
        
    def save_parameters(self, fpath):
        self.params['J'] = self.ui.valueJ.value()
        self.params['D'] = self.ui.valueD.value()
        self.params['E'] = self.ui.valueE.value()
        self.params['X'] = self.ui.valueX.value()
        self.params['B'] = self.ui.valueB.value()
        self.params['theta'] = self.ui.sliderTheta.value()
        self.params['phi'] = self.ui.sliderPhi.value()
        self.params['alpha'] = self.ui.sliderAlpha.value()
        self.params['beta'] = self.ui.sliderBeta.value()
        self.params['gamma'] = self.ui.sliderGamma.value()
        self.params['A'] = self.ui.valueA.value()
        self.params['ggamma'] = self.ui.valueGGamma.value()
        self.params['Evib'] = self.ui.valueEvib.value()
        self.params['T'] = self.ui.valueT.value()
        self.params['ET1'] = self.ui.valueET1.value()
        self.params['ET2'] = self.ui.valueET2.value()
        self.params['kS'] = self.ui.valueKS.value()
        self.params['kD'] = self.ui.valueKD.value()
        self.params['kTTA'] = self.ui.valueKTTA.value()
        self.params['kSF'] = self.ui.valueKSF.value()
        self.params['k_SF'] = self.ui.valueK_SF.value()
        self.params['kTS'] = self.ui.valueKTS.value()
        self.params['kTF'] = self.ui.valueKTF.value()
        self.params['kTFfactor'] = self.ui.valueKTFfactor.value()
        self.params['kRISC'] = self.ui.valueKRISC.value()
        self.params['kIC1'] = self.ui.valueKIC1.value()
        self.params['kIC2'] = self.ui.valueKIC2.value()
        self.params['kIC21'] = self.ui.valueKIC21.value()
        self.params['use_kTFfactor'] = int(self.ui.checkBoxKTFfactor.isChecked())
        self.params['use_gaplaw'] = int(self.ui.checkBoxEnergyGapLaw.isChecked())
        self.params['xmin'] = self.ui.valueXmin.value()
        self.params['xmax'] = self.ui.valueXmax.value()
        self.params['xvar'] = self.ui.xVariable.currentIndex()
        self.params['yvar'] = self.ui.yVariable.currentIndex()
        self.params['logx'] = int(self.ui.checkBoxLogx.isChecked())
        self.params.to_csv(fpath, header=False, index=True)
        
    def closeEvent(self, event):
        self.save_parameters('last_instance_parameters.spinstats')
        event.accept()
        
    def update_J(self):
        self.sh.J = self.ui.valueJ.value()*1e-6
        self.update()
        
    def update_D(self):
        self.sh.D = self.ui.valueD.value()*1e-6
        self.update()
        
    def update_E(self):
        self.sh.E = self.ui.valueE.value()*1e-6
        self.update()
        
    def update_X(self):
        self.sh.X = self.ui.valueX.value()*1e-9
        self.update()
        
    def update_B(self):
        self.sh.B = self.ui.valueB.value()*1e-3
        self.update()
        
    def update_theta(self):
        self.sh.theta = self.ui.sliderTheta.value()*np.pi/180
        self.ui.valueTheta.setValue(self.ui.sliderTheta.value())
        self.update()
        
    def update_theta_sb(self):
        self.sh.theta = self.ui.valueTheta.value()*np.pi/180
        self.ui.sliderTheta.setValue(self.ui.valueTheta.value())
        self.update()
        
    def update_phi(self):
        self.sh.phi = self.ui.sliderPhi.value()*np.pi/180
        self.ui.valuePhi.setValue(self.ui.sliderPhi.value())
        self.update()
        
    def update_phi_sb(self):
        self.sh.phi = self.ui.valuePhi.value()*np.pi/180
        self.ui.sliderPhi.setValue(self.ui.valuePhi.value())
        self.update()
        
    def update_alpha(self):
        self.sh.alpha = self.ui.sliderAlpha.value()*np.pi/180
        self.ui.valueAlpha.setValue(self.ui.sliderAlpha.value())
        self.update()
        
    def update_alpha_sb(self):
        self.sh.alpha = self.ui.valueAlpha.value()*np.pi/180
        self.ui.sliderAlpha.setValue(self.ui.valueAlpha.value())
        self.update()
        
    def update_beta(self):
        self.sh.beta = self.ui.sliderBeta.value()*np.pi/180
        self.ui.valueBeta.setValue(self.ui.sliderBeta.value())
        self.update()
        
    def update_beta_sb(self):
        self.sh.beta = self.ui.valueBeta.value()*np.pi/180
        self.ui.sliderBeta.setValue(self.ui.valueBeta.value())
        self.update()
        
    def update_gamma(self):
        self.sh.gamma = self.ui.sliderGamma.value()*np.pi/180
        self.ui.valueGamma.setValue(self.ui.sliderGamma.value())
        self.update()
        
    def update_gamma_sb(self):
        self.sh.gamma = self.ui.valueGamma.value()*np.pi/180
        self.ui.sliderGamma.setValue(self.ui.valueGamma.value())
        self.update()
        
    def update_kS(self):
        self.m.kS = self.ui.valueKS.value()
        self.update()
        
    def update_kD(self):
        self.m.kD = self.ui.valueKD.value()
        self.update()
        
    def update_kTTA(self):
        self.m.kTTA = self.ui.valueKTTA.value()
        self.update()
        
    def update_kSF(self):
        self.m.kSF = self.ui.valueKSF.value()
        self.update()
        
    def update_k_SF(self):
        self.m.k_SF = self.ui.valueK_SF.value()
        self.update()
        
    def update_kTS(self):
        self.m.kTS = self.ui.valueKTS.value()
        if self.use_kTF_factor:
            self.m.kTF = self.m.kTS/self.kTFfactor
            self.ui.valueKTF.setValue(self.m.kTF)
        self.update()
 
    def update_kTF(self):
        if not self.use_kTF_factor:
            self.m.kTF = self.ui.valueKTF.value()
            self.update()
        
    def update_kTFfactor(self):
        if self.use_kTF_factor:
            self.kTFfactor = self.ui.valueKTFfactor.value()
            self.m.kTF = self.m.kTS/self.kTFfactor
            self.ui.valueKTF.setValue(self.m.kTF)
            self.update()
        
    def update_use_kTFfactor(self):
        if self.ui.checkBoxKTFfactor.isChecked():
            self.use_kTF_factor = True
            self.ui.valueKTFfactor.setEnabled(True)
            self.ui.valueKTF.setEnabled(False)
            self.ui.xVariable.removeItem(self.ui.xVariable.findText('kTF'))
            self.update_kTFfactor()
        else:
            self.use_kTF_factor = False
            self.ui.valueKTFfactor.setEnabled(False)
            self.ui.valueKTF.setEnabled(True)
            self.ui.xVariable.addItem('kTF')
            self.update_kTF()
            
    def update_kRISC(self):
        self.m.kRISC = self.ui.valueKRISC.value()
        self.update()
        
    def update_use_energygap(self):
        if self.ui.checkBoxEnergyGapLaw.isChecked():
            self.use_energygap = True
            self.ui.valueA.setEnabled(True)
            self.ui.valueGGamma.setEnabled(True)
            self.ui.valueEvib.setEnabled(True)
            self.ui.valueT.setEnabled(True)
            self.ui.valueET1.setEnabled(True)
            self.ui.valueET2.setEnabled(True)
            self.ui.valueKIC1.setEnabled(False)
            self.ui.valueKIC2.setEnabled(False)
            self.ui.valueKIC21.setEnabled(False)
            self.ui.xVariable.removeItem(self.ui.xVariable.findText('kIC1'))
            self.ui.xVariable.removeItem(self.ui.xVariable.findText('kIC2'))
            self.ui.xVariable.removeItem(self.ui.xVariable.findText('kIC21'))
            self.ui.xVariable.addItems(['T', 'T1 energy', 'T2 energy'])
            self.m.gaplaw_A = self.ui.valueA.value()
            self.m.gaplaw_gamma = self.ui.valueGGamma.value()
            self.m.gaplaw_Evib = self.ui.valueEvib.value()
            self.m.gaplaw_kT = self.ui.valueT.value()*k/e
            self.ET1 = self.ui.valueET1.value()
            self.ET2 = self.ui.valueET2.value()
            self.calculate_ic_rates()
            self.update()
        else:
            self.use_energygap = False
            self.ui.valueA.setEnabled(False)
            self.ui.valueGGamma.setEnabled(False)
            self.ui.valueEvib.setEnabled(False)
            self.ui.valueT.setEnabled(False)
            self.ui.valueET1.setEnabled(False)
            self.ui.valueET2.setEnabled(False)
            self.ui.valueKIC1.setEnabled(True)
            self.ui.valueKIC2.setEnabled(True)
            self.ui.valueKIC21.setEnabled(True)
            self.ui.xVariable.removeItem(self.ui.xVariable.findText('T'))
            self.ui.xVariable.removeItem(self.ui.xVariable.findText('T1 energy'))
            self.ui.xVariable.removeItem(self.ui.xVariable.findText('T2 energy'))
            self.ui.xVariable.addItems(['kIC1', 'kIC2', 'kIC21'])
            self.m.kIC1 = self.ui.valueKIC1.value()
            self.m.kIC2 = self.ui.valueKIC2.value()
            self.m.kIC21 = self.ui.valueKIC21.value()
            self.update()
            
    def calculate_ic_rates(self):
        self.m.kIC1 = self.m.calculate_internal_conversion_rate(-self.ET1)
        self.m.kIC2 = self.m.calculate_internal_conversion_rate(self.ET2-2*self.ET1)
        self.m.kIC21 = self.m.calculate_internal_conversion_rate(self.ET1-self.ET2)
        self.ui.valueKIC1.setValue(self.m.kIC1)
        self.ui.valueKIC2.setValue(self.m.kIC2)
        self.ui.valueKIC21.setValue(self.m.kIC21)
        
    def update_kIC1(self):
        if not self.use_energygap:
           self.m.kIC1 = self.ui.valueKIC1.value()
           self.update()
           
    def update_kIC2(self):
        if not self.use_energygap:
           self.m.kIC2 = self.ui.valueKIC2.value()
           self.update()
           
    def update_kIC21(self):
        if not self.use_energygap:
           self.m.kIC21 = self.ui.valueKIC21.value()
           self.update()
           
    def update_A(self):
        if self.use_energygap:
            self.m.gaplaw_A = self.ui.valueA.value()
            self.calculate_ic_rates()
            self.update()
            
    def update_GGamma(self):
        if self.use_energygap:
            self.m.gaplaw_gamma = self.ui.valueGGamma.value()
            self.calculate_ic_rates()
            self.update()
            
    def update_Evib(self):
        if self.use_energygap:
            self.m.gaplaw_Evib = self.ui.valueEvib.value()
            self.calculate_ic_rates()
            self.update()
            
    def update_kT(self):
        if self.use_energygap:
            self.m.gaplaw_kT = self.ui.valueT.value()*k/e
            self.calculate_ic_rates()
            self.update()
            
    def update_ET1(self):
        if self.use_energygap:
            self.ET1 = self.ui.valueET1.value()
            self.calculate_ic_rates()
            self.update()
            
    def update_ET2(self):
        if self.use_energygap:
            self.ET2 = self.ui.valueET2.value()
            self.calculate_ic_rates()
            self.update()
            
    def update_logx(self):
        self.use_logx = self.ui.checkBoxLogx.isChecked()
        self.update()
        
    def update_xmin(self):
        self.xmin = self.ui.valueXmin.value()
        self.update()   
        
    def update_xmax(self):
        self.xmax = self.ui.valueXmax.value()
        self.update()
        
    def update_all_parameters(self):
        self.updates_on = False
        self.update_J()
        self.update_D()
        self.update_E()
        self.update_X()
        self.update_B()
        self.update_theta()
        self.update_phi()
        self.update_alpha()
        self.update_beta()
        self.update_gamma()
        self.update_kT()
        self.update_ET1()
        self.update_ET2()
        self.update_kS()
        self.update_kD()
        self.update_kTTA()
        self.update_kSF()
        self.update_k_SF()
        self.update_kTS()
        self.update_kTF()
        self.update_kRISC()
        self.update_kIC1()
        self.update_kIC2()
        self.update_kIC21()
        self.updates_on = True
        
    def update_xvar(self):
        self.updates_on = False
        if self.initialised:
            self.previous_xvar = self.current_xvar
            previous_xvar_ui_obj = self.xvars.loc[self.previous_xvar, 'ui']
            vars(self.ui)[previous_xvar_ui_obj].setEnabled(True)
        self.current_xvar = self.ui.xVariable.currentText()
        current_xvar_ui_obj = self.xvars.loc[self.current_xvar, 'ui']
        vars(self.ui)[current_xvar_ui_obj].setEnabled(False)
        self.ui.checkBoxEnergyGapLaw.setEnabled(True)
        self.ui.checkBoxKTFfactor.setEnabled(True)
        if self.current_xvar in ['T1 energy', 'T2 energy', 'T']:
            self.ui.checkBoxEnergyGapLaw.setEnabled(False)
        if self.current_xvar in ['kTF']:
            self.ui.checkBoxKTFfactor.setEnabled(False)
        self.ui.valueXmin.setValue(self.xvars.loc[self.current_xvar, 'xmin'])
        self.ui.valueXmax.setValue(self.xvars.loc[self.current_xvar, 'xmax'])
        self.update_xmin()
        self.update_xmax()
        self.ui.checkBoxLogx.setChecked(bool(self.xvars.loc[self.current_xvar, 'log']))
        self.update_logx()
        self.xvar_obj = self.xvars.loc[self.current_xvar, 'obj']
        self.xvar_str = self.xvars.loc[self.current_xvar, 'str']
        self.xvar_factor = self.xvars.loc[self.current_xvar, 'factor']
        self.update_all_parameters()
        self.updates_on = True
        self.update()
        
    def update_yvar(self):
        self.yvar_str = self.ui.yVariable.currentText()
        self.update()
        
    def calculate(self):
        self.sh.calculate_everything()
        self.m.set_spinhamiltonian_parameters(self.sh)
        self.m.calculate_eta()
        
    def calculate_plot_data(self):
        if self.use_logx and self.xmin > 0:
            xvarx = np.geomspace(self.xmin, self.xmax, 100)
        elif self.use_logx and self.xmin < 0:
            return False
        else:
            xvarx = np.linspace(self.xmin, self.xmax, 100)
        self.current_data = pd.DataFrame(index=xvarx, columns=['PLQY', 'UCQY', 'eta'])
        for value in xvarx:
            if self.xvar_obj == 'S':
                vars(self.sh)[self.xvar_str] = value*self.xvar_factor
            elif self.xvar_obj == 'M':
                vars(self.m)[self.xvar_str] = value*self.xvar_factor
            else:
                vars(self)[self.xvar_str] = value*self.xvar_factor
            if self.xvar_str in ['ET1', 'ET2', 'gaplaw_kT']:
                self.calculate_ic_rates()
            self.calculate()
            self.current_data.loc[value, ['PLQY', 'UCQY', 'eta']] = [self.m.plqy, self.m.ucqy, self.m.eta]
        self.current_plot_data = 100*self.current_data[self.yvar_str]
        return True
    
    def update_plot(self):
        self.ui.graph.plotItem.setLabels(left=self.yvar_str+' (%)', bottom=self.ui.xVariable.currentText())
        self.ui.graph.plotItem.plot(self.current_plot_data.index.values, self.current_plot_data.values, clear=True, pen='b')
        if self.use_logx:
            self.ui.graph.plotItem.ctrl.logXCheck.setChecked(True)
        else:
            self.ui.graph.plotItem.ctrl.logXCheck.setChecked(False)
        if self.current_plot_data.max()-self.current_plot_data.min() < 5:
            self.ui.graph.plotItem.setYRange(self.current_plot_data.min()-2.5, self.current_plot_data.max()+2.5)
        else:
            self.ui.graph.plotItem.enableAutoRange()
    
    def update(self):
        if self.updates_on:
            success = self.calculate_plot_data()
            if success:
                self.update_plot()
                self.display_status('application running', 'green')
            else:
                self.display_status('log x: x min cannot be zero', 'red', msecs=5000)
            
    def save_plot_data(self):
        fname = QtWidgets.QFileDialog().getSaveFileName(self, 'save data as', self.file_directory, 'CSV (*.csv)')[0]
        if fname != '':
            self.file_directory = os.path.dirname(fname)
            self.current_data.to_csv(fname, header=True, index=True)
            self.display_status('data saved', 'blue', msecs=5000)

 
if __name__ == '__main__':
    QtWidgets.QApplication.setStyle('Fusion')
    app = QtWidgets.QApplication(sys.argv)
    window = Application()
    window.show()
    sys.exit(app.exec_())
