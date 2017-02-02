#%% Imports
import logging
import logging.config

logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

from PyQt4 import QtGui
import sys

from optima_tb.project import Project

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
        
        self.project = None
        self.status = 'Status: No Project generated'
        self.flag_databook_visible = False
        self.flag_parset_visible = False
        self.flag_pars_visible = False
        
        self.setWindowTitle('Manual Calibration')
        
        self.label_project_name = QtGui.QLabel('Project Name: ')
        self.edit_project_name = QtGui.QLineEdit()
        self.edit_project_name.setText('New Project')
        
        self.label_cascade = QtGui.QLabel('Cascade Filename: ')
        self.edit_cascade = QtGui.QLineEdit()
        self.button_cascade = QtGui.QPushButton('Locate File', self)
        self.button_cascade.clicked.connect(self.makeSelectFile(display_field = self.edit_cascade))
        
        self.button_project_init = QtGui.QPushButton('Create Project', self)
        self.button_project_init.clicked.connect(self.createProject)
        
        self.label_databook = QtGui.QLabel('Databook Filename: ')
        self.edit_databook = QtGui.QLineEdit()
        self.button_databook = QtGui.QPushButton('Locate File', self)
        self.button_databook.clicked.connect(self.makeSelectFile(display_field = self.edit_databook))
        
        self.button_project_saturate = QtGui.QPushButton('Load Data', self)
        self.button_project_saturate.clicked.connect(self.loadData)
        
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
        
        grid.addWidget(self.label_project_name, 0, 0)
        grid.addWidget(self.edit_project_name, 0, 1)
        grid.addWidget(self.label_cascade, 1, 0)
        grid.addWidget(self.edit_cascade, 1, 1)
        grid.addWidget(self.button_cascade, 1, 2)
        grid.addWidget(self.button_project_init, 1, 3)
        grid.addWidget(self.label_databook, 2, 0)
        grid.addWidget(self.edit_databook, 2, 1)
        grid.addWidget(self.button_databook, 2, 2)
        grid.addWidget(self.button_project_saturate, 2, 3)
        grid.addWidget(self.label_parset, 3, 0)
        grid.addWidget(self.button_parset, 3, 1)
        
        self.status_bar = QtGui.QStatusBar()
        self.status_bar.showMessage(self.status)        
        
        total_layout = QtGui.QVBoxLayout()
        total_layout.addLayout(grid)
        total_layout.addStretch(1)
        total_layout.addWidget(self.status_bar)
        
        self.setLayout(total_layout)
        
        screen = QtGui.QDesktopWidget().availableGeometry()
        self.resize(screen.width()*2.0/3.0, screen.height()*2.0/3.0)
        self.setGeometry((screen.width()-self.width())/2, (screen.height()-self.height())/2, 
                         self.width(), self.height())   
    
        self.refreshVisibility()
        self.show()
        
    def refreshVisibility(self):
        self.status_bar.showMessage(self.status)
        
        self.label_databook.setVisible(self.flag_databook_visible)
        self.edit_databook.setVisible(self.flag_databook_visible)
        self.button_databook.setVisible(self.flag_databook_visible)
        self.button_project_saturate.setVisible(self.flag_databook_visible)
        
        self.label_parset.setVisible(self.flag_parset_visible)
        self.button_parset.setVisible(self.flag_parset_visible)
        
    def createProject(self):
        try:
            self.project = Project(name = self.edit_project_name.text(), cascade_path = self.edit_cascade.text(), validation_level = 'error')
            self.status = ('Status: Project "%s" generated, cascade settings loaded' % self.project.name)
            self.flag_databook_visible = True
            self.flag_parset_visible = False
            self.flag_pars_visible = False
        except:
            self.project = None
            self.status = ('Status: Attempt to generate Project failed')
            self.flag_databook_visible = False
            self.flag_parset_visible = False
            self.flag_pars_visible = False
        self.refreshVisibility()
        
    def loadData(self):
        try:
            self.project.loadSpreadsheet(databook_path = self.edit_databook.text())
            self.status = ('Status: Valid data loaded into Project "%s"' % self.project.name)
            self.flag_databook_visible = True
            self.flag_parset_visible = True
            self.flag_pars_visible = False
        except:
            self.project = None
            self.status = ('Status: Attempt to load data into Project failed, Project reset for safety')
            self.flag_databook_visible = False
            self.flag_parset_visible = False
            self.flag_pars_visible = False
        self.refreshVisibility()
        
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