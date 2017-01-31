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
        
        self.button_calibration = QtGui.QPushButton('Manual Calibration', self)
        self.button_calibration.clicked.connect(self.runGUICalibration)
        layout = QtGui.QVBoxLayout(self)
        layout.addWidget(self.button_calibration)
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
        
        self.label_cascade = QtGui.QLabel('Cascade Filename: ')
        self.edit_cascade = QtGui.QLineEdit()
        self.button_cascade = QtGui.QPushButton('Locate File', self)
        self.button_cascade.clicked.connect(self.makeSelectFile(display_field = self.edit_cascade))
        
        self.label_databook = QtGui.QLabel('Databook Filename: ')
        self.edit_databook = QtGui.QLineEdit()
        self.button_databook = QtGui.QPushButton('Locate File', self)
        self.button_databook.clicked.connect(self.makeSelectFile(display_field = self.edit_databook))
        
        self.label_parset = QtGui.QLabel('Parset: ')
        self.button_parset = QtGui.QToolButton(self)
        self.button_parset.setPopupMode(QtGui.QToolButton.MenuButtonPopup)
        self.button_parset.setMenu(QtGui.QMenu(self.button_parset))
        self.textbox_parset = QtGui.QTextBrowser(self)
        action_parset = QtGui.QWidgetAction(self.button_parset)
        action_parset.setDefaultWidget(self.textbox_parset)
        self.button_parset.menu().addAction(action_parset)
        
        grid = QtGui.QGridLayout()
        grid.setSpacing(10)          
        
        grid.addWidget(self.label_cascade, 0, 0)
        grid.addWidget(self.edit_cascade, 0, 1)
        grid.addWidget(self.button_cascade, 0, 2)
        grid.addWidget(self.label_databook, 1, 0)
        grid.addWidget(self.edit_databook, 1, 1)
        grid.addWidget(self.button_databook, 1, 2)
        grid.addWidget(self.label_parset, 2, 0)
        grid.addWidget(self.button_parset, 2, 1)
        
        self.setLayout(grid)
        
        screen = QtGui.QDesktopWidget().availableGeometry()
        self.resize(screen.width()*2.0/3.0, screen.height()*2.0/3.0)
        self.setGeometry((screen.width()-self.width())/2, (screen.height()-self.height())/2, 
                         self.width(), self.height())   
    
        self.show()
        
    def makeSelectFile(self, display_field):
        def selectFile(self):
            display_field.setText(QtGui.QFileDialog.getOpenFileName())
        return selectFile
    

#%% Functionality for opening GUIs

def runGUI():
    ''' Function that launches all available back-end GUIs as they are developed. '''
    
    app = QtGui.QApplication(sys.argv)
    app.setApplicationName('Optima GUI')
    gui = GUI()
    sys.exit(app.exec_())