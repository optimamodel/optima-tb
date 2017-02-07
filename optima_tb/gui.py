#%% Imports
import logging
import logging.config

logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from PyQt4 import QtGui, QtCore
import sys
import numpy as np
import pylab as pl
from copy import deepcopy as dcp

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
        
    def resetAttributes(self):
        self.project = None     # This is the Project object that the GUI stores to calibrate and edit.
        self.parset = None      # This is the ParameterSet object that stores all edits made in the GUI.
                                # It is an (edited) copy, not a reference, to an existing Project ParameterSet.
        
        self.parset_source_name = None
        self.parset_comparison_name = None
        self.tvec = None
        
        self.combo_dict = {}    # Dictionary that maps Parset names to indices used in combo boxes.
        self.col_par_name = 0   # Column index for table calibration parameter names.
        self.col_pop_name = 1   # Column index for table calibration population names.
        
        self.guard_status = False   # Convenient flag that locks status updates if set to true.
        self.results_current = None
        self.results_comparison = None
        
    def initUI(self):
        
        self.resetAttributes()
        
        self.status = 'Status: No Project generated'    # A status message to show at the bottom of the screen.      
        
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
        
        self.label_parset = QtGui.QLabel('Parset To Edit: ')
        self.combo_parset = QtGui.QComboBox(self)
        self.combo_parset.activated[str].connect(self.loadCalibration)
        
        self.label_compare = QtGui.QLabel('Compare Edits With... ')
        self.combo_compare = QtGui.QComboBox(self)
        self.combo_compare.activated[str].connect(self.selectComparison)
        self.button_compare = QtGui.QPushButton('Compare Models', self)
        self.button_compare.clicked.connect(self.runComparison)
        
        self.label_overwrite = QtGui.QLabel('Save Edits To... ')
        self.edit_overwrite = QtGui.QLineEdit()
        self.button_overwrite = QtGui.QPushButton('Save Calibration', self)
        self.button_overwrite.clicked.connect(self.saveCalibration)
        
        self.status_bar = QtGui.QStatusBar()
        self.status_bar.showMessage(self.status)
        self.status_bar.setSizeGripEnabled(False)
        
        self.table_calibration = QtGui.QTableWidget()
        self.table_calibration.cellChanged.connect(self.updateParset)
        policy_min = QtGui.QSizePolicy.Minimum
        policy_exp = QtGui.QSizePolicy.Expanding
        self.parset_layout_stretch = QtGui.QSpacerItem(0, 0, policy_min, policy_exp)
        
        self.table_plotter = QtGui.QTableWidget()
        
        # Layout.   
        grid_upper = QtGui.QGridLayout()
        grid_upper.setSpacing(10)          
        
        grid_upper.addWidget(self.label_project_name, 0, 0)
        grid_upper.addWidget(self.edit_project_name, 0, 1)
        grid_upper.addWidget(self.label_cascade, 1, 0)
        grid_upper.addWidget(self.edit_cascade, 1, 1)
        grid_upper.addWidget(self.button_cascade, 1, 2)
        grid_upper.addWidget(self.button_project_init, 1, 3)
        grid_upper.addWidget(self.label_databook, 2, 0)
        grid_upper.addWidget(self.edit_databook, 2, 1)
        grid_upper.addWidget(self.button_databook, 2, 2)
        grid_upper.addWidget(self.button_project_saturate, 2, 3)
        grid_upper.addWidget(self.label_parset, 3, 0)
        grid_upper.addWidget(self.combo_parset, 3, 1)
        
        grid_lower = QtGui.QGridLayout()
        grid_lower.setSpacing(10)
        
        grid_lower.addWidget(self.label_compare, 0, 0)
        grid_lower.addWidget(self.combo_compare, 0, 1)
        grid_lower.addWidget(self.button_compare, 0, 2)
        grid_lower.addWidget(self.label_overwrite, 1, 0)
        grid_lower.addWidget(self.edit_overwrite, 1, 1)
        grid_lower.addWidget(self.button_overwrite, 1, 2)
        
        parset_layout = QtGui.QVBoxLayout()
        parset_layout.addLayout(grid_upper)
        parset_layout.addWidget(self.table_calibration)
        parset_layout.addLayout(grid_lower)
        parset_layout.addItem(self.parset_layout_stretch)
        parset_layout.addWidget(self.status_bar)
        
        self.plotter_layout = QtGui.QGridLayout()
        self.plotter_layout.setSpacing(5)
#        plotter_layout = QtGui.QVBoxLayout()
#        plotter_layout.addWidget(self.table_plotter)
        
        self.first_half = QtGui.QWidget()
        self.first_half.resize(self.width()/2.0, self.height())
        self.first_half.setLayout(parset_layout)
        self.second_half = QtGui.QWidget()
        self.second_half.resize(self.width()/2.0, self.height()/2.0)
        self.second_half.setLayout(self.plotter_layout)
        
        self.splitter_total = QtGui.QSplitter()
        self.splitter_total.addWidget(self.first_half)
        self.splitter_total.addWidget(self.second_half)
        total_layout = QtGui.QHBoxLayout()
        total_layout.addWidget(self.splitter_total)
        self.setLayout(total_layout)
    
        self.refreshVisibility()
        self.show()
        
    def refreshVisibility(self):
        self.refreshStatus()

        is_cascade_loaded = self.project is not None
        is_parset_loaded = is_cascade_loaded and len(self.project.parsets.keys()) > 0
        
        self.label_databook.setVisible(is_cascade_loaded)
        self.edit_databook.setVisible(is_cascade_loaded)
        self.button_databook.setVisible(is_cascade_loaded)
        self.button_project_saturate.setVisible(is_cascade_loaded)
        
        if is_parset_loaded:
            self.refreshComboBoxes()
        self.label_parset.setVisible(is_parset_loaded)
        self.combo_parset.setVisible(is_parset_loaded)
        
        policy_min = QtGui.QSizePolicy.Minimum
        policy_exp = QtGui.QSizePolicy.Expanding
        if is_parset_loaded:
            self.parset_layout_stretch.changeSize(0, 0, policy_min, policy_min)
            self.makeParsetTable()
        else:
            self.parset_layout_stretch.changeSize(0, 0, policy_min, policy_exp)
        self.table_calibration.setVisible(is_parset_loaded)
        
        self.label_compare.setVisible(is_parset_loaded)
        self.combo_compare.setVisible(is_parset_loaded)
        self.button_compare.setVisible(is_parset_loaded)
        self.label_overwrite.setVisible(is_parset_loaded)
        self.edit_overwrite.setVisible(is_parset_loaded)
        self.button_overwrite.setVisible(is_parset_loaded)
    
    # NOTE: This constant refreshing of comboboxes is inefficient and can be improved later.
    def refreshComboBoxes(self):
        combo_boxes = [self.combo_parset, self.combo_compare]
        combo_names = [self.parset_source_name, self.parset_comparison_name]
        self.combo_dict = {}
        for k in xrange(len(combo_boxes)):
            combo_box = combo_boxes[k]
            combo_name = combo_names[k]
            combo_box.clear()
            cid = 0
            for parset_name in self.project.parsets.keys():
                self.combo_dict[parset_name] = cid
                combo_box.addItem(parset_name)
                cid += 1
            combo_box.setCurrentIndex(self.combo_dict[combo_name])
            
    def refreshStatus(self):
        if not self.guard_status:
            self.status_bar.showMessage(self.status)
            
        
    def createProject(self):
        try:
            self.project = Project(name = self.edit_project_name.text(), cascade_path = self.edit_cascade.text(), validation_level = 'avert')
            self.tvec = np.arange(self.project.settings.tvec_start, self.project.settings.tvec_observed_end + 1.0/2)
            self.status = ('Status: Project "%s" generated, cascade settings loaded' % self.project.name)
        except:
            self.resetAttributes()
            self.status = ('Status: Attempt to generate Project failed')
        self.refreshVisibility()
        
    def loadData(self):
        try:
            self.project.loadSpreadsheet(databook_path = self.edit_databook.text())
            self.project.resetParsets()
            self.project.makeParset(name = 'default')
            self.loadCalibration(self.project.parsets[0].name, delay_refresh = True)
            self.selectComparison(self.project.parsets[0].name)
            self.status = ('Status: Valid data loaded into Project "%s", default Parset generated' % self.project.name)
        except:
            self.resetAttributes()
            self.status = ('Status: Attempt to load data into Project failed, Project reset for safety')
        self.refreshVisibility()
        
    def loadCalibration(self, parset_name, delay_refresh = False):
        self.parset_source_name = parset_name
        self.parset = dcp(self.project.parsets[parset_name])
        self.status = ('Status: Parset "%s" selected for editing' % parset_name)
        if not delay_refresh:
            self.refreshVisibility()
        
    def saveCalibration(self):
        parset_name = self.edit_overwrite.text()
        if parset_name == '':
            self.status = ('Status: Attempt to save Parset failed, no name provided')
        else:
            if parset_name in self.project.parsets.keys():
                self.status = ('Status: Parset "%s" successfully overwritten' % parset_name)
            else:
                self.status = ('Status: New Parset "%s" added to Project' % parset_name)
            self.project.parsets[parset_name] = dcp(self.parset)
        self.refreshVisibility()
        
    def selectComparison(self, parset_name):
        self.parset_comparison_name = parset_name
        
    def runComparison(self):
        self.status = ('Status: Running models for Parset comparison')
        self.refreshStatus()
        self.results_current = self.project.runSim(parset = self.parset)
        self.status = ('Status: Model successfully processed for edited Parset')
        self.refreshStatus()
        self.results_comparison = self.project.runSim(parset_name = self.parset_comparison_name)
        self.status = ('Status: Model successfully processed for Parset "%s"' % self.parset_comparison_name)
        self.refreshStatus()
        self.makePlotterTable()
        self.results_current = None
        self.results_comparison = None
        return
        
    def factorySelectFile(self, display_field):
        def selectFile(self):
            display_field.setText(QtGui.QFileDialog.getOpenFileName())
        return selectFile
        
    def makeParsetTable(self):
        self.table_calibration.setVisible(False)    # Resizing columns requires table to be hidden first.
        self.table_calibration.clear()
        
        # Disconnect the calibration table from cell change signals to avoid signal flooding during connection.
        try: self.table_calibration.cellChanged.disconnect()
        except: pass

#        parset = self.project.parsets[self.parset_source_name]
        parset = self.parset
        num_pops = len(parset.pop_labels)
        row_count = num_pops*(len(parset.pars['cascade'])-len(self.project.settings.par_funcs))
        self.table_calibration.setRowCount(row_count)
        self.table_calibration.setColumnCount(2+len(self.tvec))
        self.calibration_items = []
        
        k = 0
        par_labels = []
        for par_type in ['characs','cascade']:
            for par in parset.pars[par_type]:
                if ((par_type == 'cascade' and par.label not in self.project.settings.par_funcs.keys()) or (par_type == 'characs' and 'entry_point' in self.project.settings.charac_specs[par.label].keys())):
                    for pid in xrange(len(parset.pop_labels)):
                        pop_label = parset.pop_labels[pid]
                        par_labels.append(par.label+' ['+pop_label+']')
                        try:
                            par_name = self.project.settings.linkpar_specs[par.label]['name']
                        except:
                            par_name = self.project.settings.charac_specs[par.label]['name']
                        temp = QtGui.QTableWidgetItem()
                        temp.setText(par_name)
                        temp.setFlags(QtCore.Qt.ItemIsEnabled or QtCore.Qt.ItemIsSelectable)
                        self.table_calibration.setItem(k*num_pops+pid, self.col_par_name, temp)
                        temp = QtGui.QTableWidgetItem()
                        temp.setText(parset.pop_names[pid])
                        temp.setFlags(QtCore.Qt.ItemIsEnabled or QtCore.Qt.ItemIsSelectable)
                        self.table_calibration.setItem(k*num_pops+pid, self.col_pop_name, temp)
                        
                        for eid in xrange(len(par.t[pid])):
                            t = par.t[pid][eid]
                            y = par.y[pid][eid]
                            temp = QtGui.QTableWidgetItem()
                            temp.setText(str(y))
                            self.table_calibration.setItem(k*num_pops+pid, 2+int(t)-self.tvec[0], temp)
                    k += 1
        self.table_calibration.setVerticalHeaderLabels(par_labels)
        self.table_calibration.setHorizontalHeaderLabels(['Par. Name','Pop. Name']+[str(int(x)) for x in self.tvec])
        self.table_calibration.resizeColumnsToContents()
        
        self.table_calibration.cellChanged.connect(self.updateParset)
        
    def makePlotterTable(self):
        
        self.figure = pl.Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)

        self.axis = self.figure.add_subplot(111)
        self.axis.scatter([1,2,3],[4,5,6])        
        
        self.plotter_layout.addWidget(self.canvas, 0, 0)
        
        
#        self.table_calibration.setVisible(False)    # Resizing columns requires table to be hidden first.
#        self.table_calibration.clear()
        
        
        
#        # Disconnect the calibration table from cell change signals to avoid signal flooding during connection.
#        try: self.table_calibration.cellChanged.disconnect()
#        except: pass
#
##        parset = self.project.parsets[self.parset_source_name]
#        parset = self.parset
#        num_pops = len(parset.pop_labels)
#        row_count = num_pops*(len(parset.pars['cascade'])-len(self.project.settings.par_funcs))
#        self.table_calibration.setRowCount(row_count)
#        self.table_calibration.setColumnCount(2+len(self.tvec))
#        self.calibration_items = []
#        
#        k = 0
#        par_labels = []
#        for par_type in ['characs','cascade']:
#            for par in parset.pars[par_type]:
#                if ((par_type == 'cascade' and par.label not in self.project.settings.par_funcs.keys()) or (par_type == 'characs' and 'entry_point' in self.project.settings.charac_specs[par.label].keys())):
#                    for pid in xrange(len(parset.pop_labels)):
#                        pop_label = parset.pop_labels[pid]
#                        par_labels.append(par.label+' ['+pop_label+']')
#                        try:
#                            par_name = self.project.settings.linkpar_specs[par.label]['name']
#                        except:
#                            par_name = self.project.settings.charac_specs[par.label]['name']
#                        temp = QtGui.QTableWidgetItem()
#                        temp.setText(par_name)
#                        temp.setFlags(QtCore.Qt.ItemIsEnabled or QtCore.Qt.ItemIsSelectable)
#                        self.table_calibration.setItem(k*num_pops+pid, self.col_par_name, temp)
#                        temp = QtGui.QTableWidgetItem()
#                        temp.setText(parset.pop_names[pid])
#                        temp.setFlags(QtCore.Qt.ItemIsEnabled or QtCore.Qt.ItemIsSelectable)
#                        self.table_calibration.setItem(k*num_pops+pid, self.col_pop_name, temp)
#                        
#                        for eid in xrange(len(par.t[pid])):
#                            t = par.t[pid][eid]
#                            y = par.y[pid][eid]
#                            temp = QtGui.QTableWidgetItem()
#                            temp.setText(str(y))
#                            self.table_calibration.setItem(k*num_pops+pid, 2+int(t)-self.tvec[0], temp)
#                    k += 1
#        self.table_calibration.setVerticalHeaderLabels(par_labels)
#        self.table_calibration.setHorizontalHeaderLabels(['Par. Name','Pop. Name']+[str(int(x)) for x in self.tvec])
#        self.table_calibration.resizeColumnsToContents()
#        
#        self.table_calibration.cellChanged.connect(self.updateParset)
        
    def updateParset(self, row, col):
        new_val_str = str(self.table_calibration.item(row,col).text())
        year = float(str(self.table_calibration.horizontalHeaderItem(col).text()))
        
        par_name = str(self.table_calibration.item(row, self.col_par_name).text())
        pop_name = str(self.table_calibration.item(row, self.col_pop_name).text())
        try:
            par_label = self.project.settings.linkpar_name_labels[par_name]
        except:
            par_label = self.project.settings.charac_name_labels[par_name]
        pop_label = self.project.data['pops']['name_labels'][pop_name]
        
        par = self.parset.getPar(par_label)
        try:
            new_val = float(new_val_str)
        except:
            if new_val_str == '':
                remove_success = par.removeValueAt(t = year, pop_label = pop_label)
                if remove_success:
                    self.status = ('Status: Value successfully deleted from Parset')
                    self.refreshStatus()
                    return
                else: self.status = ('Status: Attempt to remove item in Parset failed, at least one value per row required')
            else: self.status = ('Status: Attempt to edit item in Parset failed, only numbers allowed')
            self.refreshStatus()
            self.guard_status = True
            self.table_calibration.item(row,col).setText(str(par.interpolate(tvec = [year], pop_label = pop_label)[0]))
            self.guard_status = False
            return
        par.insertValuePair(t = year, y = new_val, pop_label = pop_label)
        self.status = ('Status: Parset "%s" given value "%f" for parameter "%s", population "%s", year "%i"' % (self.parset.name, new_val, par_label, pop_label, year))
        self.refreshStatus()
            
#        if not new_val_str == '':
#            try:
#                new_val = float(new_val_str)
#            except:
#                self.table_calibration.item(row,col).setText('')
#                self.status = ('Status: Attempt to edit item in Parset failed, only numbers allowed')
#                self.refreshStatus()
#                return
#            par.insertValuePair(t = year, y = new_val, pop_label = pop_label)
#        else:
#            remove_success = par.removeValueAt(t = year, pop_label = pop_label)
#            if not remove_success:
#                self.table_calibration.item(row,col).setText(par.interpolate(tvec = [year], pop_label = pop_label)[0])
#                self.status = ('Status: Attempt to remove item in Parset failed, at least one value per row required')
#                self.refreshStatus()
#                return
##            except:
##                self.table_calibration.item(row,col).setText('')
##                self.status = ('Status: Attempt to edit item in Parset failed, Parset may be formatted incorrectly')
##                self.refreshStatus()
#        return



#%% Functionality for opening GUIs

def runGUI():
    ''' Function that launches all available back-end GUIs as they are developed. '''
    
    app = QtGui.QApplication(sys.argv)
    app.setApplicationName('Optima GUI')
    gui = GUI()
    sys.exit(app.exec_())