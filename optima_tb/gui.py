#%% Imports
import logging
import logging.config

logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

from PyQt4 import QtGui
import sys

#%% GUI classes

class GUI(QtGui.QWidget):
    
    def __init__(self):
        super(GUI, self).__init__()
        self.initUI()        
        
    def initUI(self):
        
        self.setWindowTitle('GUI Selection Screen')   
        
        self.button = QtGui.QPushButton('Manual Calibration', self)
        self.button.clicked.connect(self.runGUICalibration)
        layout = QtGui.QVBoxLayout(self)
        layout.addWidget(self.button)
        self.setLayout(layout)
        
        screen = QtGui.QDesktopWidget().availableGeometry()
        self.setGeometry((screen.width()-self.width())/2, (screen.height()-self.height())/2, 
                         self.width(), self.height())   
    
        self.show()
        
    def runGUICalibration(self):
        self.sub_gui = GUICalibration()
        

class GUICalibration(QtGui.QWidget):
    
    def __init__(self):
        super(GUICalibration, self).__init__()
        self.initUI()
        
    def initUI(self):
        
        self.setWindowTitle('Manual Calibration')   
        
#        self.button = QtGui.QPushButton('Manual Calibration', self)
#        self.button.clicked.connect(self.runGUICalibration)
#        layout = QtGui.QVBoxLayout(self)
#        layout.addWidget(self.button)
#        self.setLayout(layout)
        
        screen = QtGui.QDesktopWidget().availableGeometry()
        self.setGeometry((screen.width()-self.width())/2, (screen.height()-self.height())/2, 
                         self.width(), self.height())   
    
        self.show()
    

#%% Functionality for opening GUIs

def runGUI():
    ''' Function that launches all available back-end GUIs as they are developed. '''
    
    app = QtGui.QApplication(sys.argv)
    app.setApplicationName('Optima GUI')
    gui = GUI()
    sys.exit(app.exec_())