#%% Imports
import logging
import logging.config

logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

from PyQt4 import QtGui, QtCore
import sys
import numpy as np

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
        self.selected_parset = None
        self.tvec = None        
        
#        self.calibration_items = []
        
        self.setWindowTitle('Manual Calibration')
        
        # Screen.
        screen = QtGui.QDesktopWidget().availableGeometry()
        self.resize(screen.width()*9.0/10.0, screen.height()*9.0/10.0)
        self.setGeometry((screen.width()-self.width())/2, (screen.height()-self.height())/2, 
                         self.width(), self.height())
        
        # Widgets.
        self.label_project_name = QtGui.QLabel('Project Name: ')
        self.edit_project_name = QtGui.QLineEdit()
        self.edit_project_name.setText('New Project')
        
        self.label_cascade = QtGui.QLabel('Cascade Filename: ')
        self.edit_cascade = QtGui.QLineEdit()
        self.button_cascade = QtGui.QPushButton('Locate File', self)
        self.button_cascade.clicked.connect(self.factorySelectFile(display_field = self.edit_cascade))
        
        self.button_project_init = QtGui.QPushButton('Create Project', self)
        self.button_project_init.clicked.connect(self.createProject)
        
        self.label_databook = QtGui.QLabel('Databook Filename: ')
        self.edit_databook = QtGui.QLineEdit()
        self.button_databook = QtGui.QPushButton('Locate File', self)
        self.button_databook.clicked.connect(self.factorySelectFile(display_field = self.edit_databook))
        
        self.button_project_saturate = QtGui.QPushButton('Load Data', self)
        self.button_project_saturate.clicked.connect(self.loadData)
        
        self.label_parset = QtGui.QLabel('Parset: ')
        self.combo_parset = QtGui.QComboBox(self)
        self.combo_parset.activated[str].connect(self.loadCalibration)
        
        self.status_bar = QtGui.QStatusBar()
        self.status_bar.showMessage(self.status)
        
        self.table_calibration = QtGui.QTableWidget()
        policy_min = QtGui.QSizePolicy.Minimum
        policy_exp = QtGui.QSizePolicy.Expanding
        self.parset_layout_stretch = QtGui.QSpacerItem(0, 0, policy_min, policy_exp)
        
        # Layout.   
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
        grid.addWidget(self.combo_parset, 3, 1)        
        
        parset_layout = QtGui.QVBoxLayout()
        parset_layout.addLayout(grid)
        parset_layout.addWidget(self.table_calibration)
        parset_layout.addItem(self.parset_layout_stretch)
        parset_layout.addWidget(self.status_bar)
        
        self.first_half = QtGui.QWidget()
        self.first_half.resize(self.width()/2.0, self.height())
        self.first_half.setLayout(parset_layout)
        self.second_half = QtGui.QWidget()
        self.second_half.resize(self.width()/2.0, self.height()/2.0)
        
        self.splitter_total = QtGui.QSplitter()
        self.splitter_total.addWidget(self.first_half)
        self.splitter_total.addWidget(self.second_half)
        total_layout = QtGui.QHBoxLayout()
        total_layout.addWidget(self.splitter_total)
        self.setLayout(total_layout)
    
        self.refreshVisibility()
        self.show()
        
    def refreshVisibility(self):
        self.status_bar.showMessage(self.status)
        
        self.label_databook.setVisible(self.flag_databook_visible)
        self.edit_databook.setVisible(self.flag_databook_visible)
        self.button_databook.setVisible(self.flag_databook_visible)
        self.button_project_saturate.setVisible(self.flag_databook_visible)
        
        self.label_parset.setVisible(self.flag_parset_visible)
        self.combo_parset.setVisible(self.flag_parset_visible)
        if self.flag_parset_visible:
            self.combo_parset.clear()
            for parset_name in self.project.parsets:
                self.combo_parset.addItem(parset_name)
        
        policy_min = QtGui.QSizePolicy.Minimum
        policy_exp = QtGui.QSizePolicy.Expanding
        if self.flag_pars_visible:
            self.parset_layout_stretch.changeSize(0, 0, policy_min, policy_min)
            self.makeParsetTable()
        else:
            self.parset_layout_stretch.changeSize(0, 0, policy_min, policy_exp)
        self.table_calibration.setVisible(self.flag_pars_visible)
        
    def createProject(self):
        try:
            self.project = Project(name = self.edit_project_name.text(), cascade_path = self.edit_cascade.text(), validation_level = 'error')
            self.tvec = np.arange(self.project.settings.tvec_start, self.project.settings.tvec_observed_end + 1.0/2)
            self.status = ('Status: Project "%s" generated, cascade settings loaded' % self.project.name)
            self.flag_databook_visible = True
            self.flag_parset_visible = False
            self.flag_pars_visible = False
        except:
            self.project = None
            self.tvec = None
            self.selected_parset = None
            self.status = ('Status: Attempt to generate Project failed')
            self.flag_databook_visible = False
            self.flag_parset_visible = False
            self.flag_pars_visible = False
        self.refreshVisibility()
        
    def loadData(self):
        try:
            self.project.loadSpreadsheet(databook_path = self.edit_databook.text())
            self.project.makeParset(name = 'default')
            self.selected_parset = 'default'
            self.status = ('Status: Valid data loaded into Project "%s", default parset generated' % self.project.name)
            self.flag_databook_visible = True
            self.flag_parset_visible = True
            self.flag_pars_visible = True
        except:
            self.project = None
            self.selected_parset = None
            self.status = ('Status: Attempt to load data into Project failed, Project reset for safety')
            self.flag_databook_visible = False
            self.flag_parset_visible = False
            self.flag_pars_visible = False
        self.refreshVisibility()
        
    def loadCalibration(self, parset_name):
        self.flag_pars_visible = True
        self.selected_parset = parset_name
        self.refreshVisibility()
        
    def factorySelectFile(self, display_field):
        def selectFile(self):
            display_field.setText(QtGui.QFileDialog.getOpenFileName())
        return selectFile
        
    def makeParsetTable(self):
        self.table_calibration.setVisible(False)    # Resizing columns requires table to be hidden first.

        parset = self.project.parsets[self.selected_parset]
        num_pops = len(parset.pop_labels)
        row_count = num_pops*(len(parset.pars['cascade'])-len(self.project.settings.par_funcs))
        self.table_calibration.setRowCount(row_count)
        self.table_calibration.setColumnCount(2+len(self.tvec))
        self.calibration_items = []
        
        k = 0
        par_labels = []
        for par in parset.pars['cascade']:
            if par.label not in self.project.settings.par_funcs:
                for pid in xrange(len(parset.pop_labels)):
                    pop_label = parset.pop_labels[pid]
                    par_labels.append(par.label+' ['+pop_label+']')
                    par_name = self.project.settings.linkpar_specs[par.label]['name']
                    temp = QtGui.QTableWidgetItem()
                    temp.setText(par_name)
                    temp.setFlags(QtCore.Qt.ItemIsEnabled or QtCore.Qt.ItemIsSelectable)
                    self.table_calibration.setItem(k*num_pops+pid, 0, temp)
                    temp = QtGui.QTableWidgetItem()
                    temp.setText(parset.pop_names[pid])
                    temp.setFlags(QtCore.Qt.ItemIsEnabled or QtCore.Qt.ItemIsSelectable)
                    self.table_calibration.setItem(k*num_pops+pid, 1, temp)
                    
                    for eid in xrange(len(par.t[pid])):
                        t = par.t[pid][eid]
                        y = par.y[pid][eid]
                        temp = QtGui.QTableWidgetItem()
                        temp.setText(str(y))
    #                    temp.setFlags(QtCore.Qt.ItemIsEnabled or QtCore.Qt.ItemIsEditable or QtCore.Qt.ItemIsSelectable)
                        self.table_calibration.setItem(k*num_pops+pid, 2+int(t)-self.tvec[0], temp)
                k += 1
        self.table_calibration.setVerticalHeaderLabels(par_labels)
        self.table_calibration.setHorizontalHeaderLabels(['Par. Name','Pop. Name']+[str(int(x)) for x in self.tvec])
        self.table_calibration.resizeColumnsToContents()



#%% Functionality for opening GUIs

def runGUI():
    ''' Function that launches all available back-end GUIs as they are developed. '''
    
    app = QtGui.QApplication(sys.argv)
    app.setApplicationName('Optima GUI')
    gui = GUI()
    sys.exit(app.exec_())