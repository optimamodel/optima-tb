#%% Imports
import logging
import logging.config

logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

from matplotlib import pyplot as pp
pp.ioff()   # Turn off interactive mode.

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
import sys
import numpy as np
from copy import deepcopy as dcp

from optima_tb.project import Project
from optima_tb.plotting import extractCharacteristic, _plotLine
from optima_tb.dataio import saveObject, loadObject

#%% GUI classes

class GUI(qtw.QWidget):
    
    def __init__(self):
        super(GUI, self).__init__()
        self.initUI()        
        
    def initUI(self):
        
        self.setWindowTitle('GUI Selection Screen')   
        
        self.button_calibration = qtw.QPushButton('Manual Calibration', self)
        self.button_calibration.clicked.connect(self.runGUICalibration)
        layout = qtw.QVBoxLayout(self)
        layout.addWidget(self.button_calibration)
        self.setLayout(layout)
        
        screen = qtw.QDesktopWidget().availableGeometry()
        self.setGeometry((screen.width()-self.width())/2, (screen.height()-self.height())/2, 
                         self.width(), self.height())   
    
        self.show()
        
    def runGUICalibration(self):
        self.sub_gui = GUICalibration()
        

class GUICalibration(qtw.QMainWindow):
    
    def __init__(self):
        super(GUICalibration, self).__init__()
        self.initUI()
        
    def resetAttributes(self):
        self.project = None     # This is the Project object that the GUI stores to calibrate and edit.
        self.parset = None      # This is the ParameterSet object that stores all edits made in the GUI.
                                # It is an (edited) copy, not a reference, to an existing Project ParameterSet.
        
        self.parset_source_name = None
        self.parset_comparison_name = None
        self.charac_plot_name = None
        self.pop_plot_name = None
        self.tvec = None
        
        self.combo_parset_dict = {}     # Dictionary that maps Parset names to indices used in combo boxes.
        self.combo_charac_dict = {}     # Dictionary that maps Characteristic names to indices used in combo boxes.
        self.combo_pop_dict = {}        # Dictionary that maps Population names to indices used in combo boxes.
        self.col_par_name = 0   # Column index for table calibration parameter names.
        self.col_pop_name = 1   # Column index for table calibration population names.
        
        self.guard_status = False   # Convenient flag that locks status updates if set to true.
        self.results_current = None
        self.results_comparison = None
        
    def initUI(self):
        
        self.resetAttributes()
        
        self.status = 'Status: No Project generated'    # A status message to show at the bottom of the screen.      
        
        self.setWindowTitle('Manual Calibration')
        
        # Screen.
        screen = qtw.QDesktopWidget().availableGeometry()
        self.resize(screen.width()*9.0/10.0, screen.height()*9.0/10.0)
        self.setGeometry((screen.width()-self.width())/2, (screen.height()-self.height())/2, 
                         self.width(), self.height())
        
        # Widgets.
        self.label_parset = qtw.QLabel('Parset To Edit: ')
        self.combo_parset = qtw.QComboBox(self)
        self.combo_parset.activated[str].connect(self.loadCalibration)
        
        self.label_compare = qtw.QLabel('Compare Edits With... ')
        self.combo_compare = qtw.QComboBox(self)
        self.combo_compare.activated[str].connect(self.selectComparison)
        self.button_compare = qtw.QPushButton('Compare Models', self)
        self.button_compare.clicked.connect(self.runComparison)
        
        self.label_overwrite = qtw.QLabel('Save Edits To... ')
        self.edit_overwrite = qtw.QLineEdit()
        self.button_overwrite = qtw.QPushButton('Save Calibration', self)
        self.button_overwrite.clicked.connect(self.saveCalibration)
        
        self.status_bar = qtw.QStatusBar()
        self.status_bar.showMessage(self.status)
        self.status_bar.setSizeGripEnabled(False)
        
        self.table_calibration = qtw.QTableWidget()
        self.table_calibration.cellChanged.connect(self.updateParset)
        policy_min = qtw.QSizePolicy.Minimum
        policy_exp = qtw.QSizePolicy.Expanding
        self.parset_layout_stretch = qtw.QSpacerItem(0, 0, policy_min, policy_exp)
        
#        self.table_plotter = qtw.QTableWidget()
        self.label_plotter_charac = qtw.QLabel('Select Characteristic: ')
        self.combo_plotter_charac = qtw.QComboBox(self)
        self.combo_plotter_charac.activated[str].connect(self.selectCharacteristic)
        self.label_plotter_pop = qtw.QLabel('Select Population: ')
        self.combo_plotter_pop = qtw.QComboBox(self)
        self.combo_plotter_pop.activated[str].connect(self.selectPopulation)
        self.button_plotter = qtw.QPushButton('Plot Figures', self)
        self.button_plotter.clicked.connect(self.showPlots)
        
        self.scroller_plotter = qtw.QScrollArea()
        self.scroller_plotter.setWidgetResizable(False)
        self.plotter = qtw.QWidget()
        
        # Layout.   
        grid_upper = qtw.QGridLayout()
        grid_upper.setSpacing(10)          
        
        grid_upper.addWidget(self.label_parset, 0, 0)
        grid_upper.addWidget(self.combo_parset, 0, 1)
        
        grid_lower = qtw.QGridLayout()
        grid_lower.setSpacing(10)
        
        grid_lower.addWidget(self.label_overwrite, 0, 0)
        grid_lower.addWidget(self.edit_overwrite, 0, 1)
        grid_lower.addWidget(self.button_overwrite, 0, 2)
        grid_lower.addWidget(self.label_compare, 1, 0)
        grid_lower.addWidget(self.combo_compare, 1, 1)
        grid_lower.addWidget(self.button_compare, 1, 2)
        grid_lower.addWidget(self.label_plotter_charac, 2, 0)
        grid_lower.addWidget(self.combo_plotter_charac, 2, 1)
        grid_lower.addWidget(self.label_plotter_pop, 3, 0)
        grid_lower.addWidget(self.combo_plotter_pop, 3, 1)
        grid_lower.addWidget(self.button_plotter, 3, 2)
        
        parset_layout = qtw.QVBoxLayout()
        parset_layout.addLayout(grid_upper)
        parset_layout.addWidget(self.table_calibration)
        parset_layout.addLayout(grid_lower)
        parset_layout.addItem(self.parset_layout_stretch)
#        parset_layout.addWidget(self.status_bar)
        
        self.plotter_layout = qtw.QGridLayout()#.QVBoxLayout()
        self.plotter_layout.setSpacing(5)
        self.scroller_layout = qtw.QVBoxLayout()
#        plotter_layout.addWidget(self.table_plotter)
        
        self.first_half = qtw.QWidget()
        self.first_half.resize(self.width()/2.0, self.height())
        self.first_half.setLayout(parset_layout)
        self.second_half = qtw.QWidget()
        self.second_half.resize(self.width()/2.0, self.height()/2.0)
        self.second_half.setLayout(self.plotter_layout)
#        self.second_half.setLayout(self.scroller_layout)
#        self.scroller_layout.addWidget(self.scroller_plotter)
#        self.scroller_plotter.setWidget(self.plotter)
#        self.plotter.setLayout(self.plotter_layout)
        
        self.splitter_total = qtw.QSplitter()
        self.splitter_total.addWidget(self.first_half)
        self.splitter_total.addWidget(self.second_half)
        
        self.setCentralWidget(self.splitter_total)
        
        # Menu.
        menu_bar = self.menuBar()
        self.file_menu = menu_bar.addMenu('Project')
        self.action_new_proj = qtw.QAction('Create', self)
        self.action_import_proj = qtw.QAction('Import', self)
        self.action_export_proj = qtw.QAction('Export', self)
        self.file_menu.addAction(self.action_new_proj)
        self.file_menu.addAction(self.action_import_proj)
        self.file_menu.addAction(self.action_export_proj)
        
        self.action_new_proj.triggered.connect(self.createProject)
        self.action_import_proj.triggered.connect(self.importProject)
        self.action_export_proj.triggered.connect(self.exportProject)
    
        self.refreshVisibility()
        self.show()
        
    def createProject(self):
        self.resetAttributes()
        project_name, project_name_chosen = qtw.QInputDialog.getText(self, 'Create Project', 'Enter project name:')
        if project_name_chosen:
            cascade_path = qtw.QFileDialog.getOpenFileName(self, 'Select Cascade File')[0]
            try:
                self.project = Project(name = project_name, cascade_path = cascade_path, validation_level = 'avert')
                self.tvec = np.arange(self.project.settings.tvec_start, self.project.settings.tvec_observed_end + 1.0/2)
                self.status = ('Status: Project "%s" generated, cascade settings loaded' % self.project.name)
            except:
                self.status = ('Status: Attempt to generate Project failed')
                self.refreshVisibility()
                return
                
            databook_path = qtw.QFileDialog.getOpenFileName(self, 'Select Databook File')[0]
            try:
                self.project.loadSpreadsheet(databook_path = databook_path)
                self.project.resetParsets()
                self.project.makeParset(name = 'default')
                self.loadCalibration(self.project.parsets[0].name, delay_refresh = True)
                self.selectComparison(self.project.parsets[0].name)
                self.status = ('Status: Valid data loaded into Project "%s", default Parset generated' % self.project.name)
            except:
                self.resetAttributes()
                self.status = ('Status: Attempt to load data into Project failed, Project reset for safety')
        self.refreshVisibility()
        
    def importProject(self):
        self.resetAttributes()
        import_path = qtw.QFileDialog.getOpenFileName(self, 'Import Project From File')[0]
        try:
            self.project = loadObject(filename=import_path)
            self.tvec = np.arange(self.project.settings.tvec_start, self.project.settings.tvec_observed_end + 1.0/2)
            if len(self.project.parsets) < 1:
                self.project.makeParset(name = 'default')
            self.loadCalibration(self.project.parsets[0].name, delay_refresh = True)
            self.selectComparison(self.project.parsets[0].name)
            self.status = ('Status: Project "%s" successfully imported' % self.project.name)
        except:
            self.status = ('Status: Attempt to import Project failed')
        self.refreshVisibility()
        
    def exportProject(self):
        export_path = qtw.QFileDialog.getSaveFileName(self, 'Export Project To File')[0]
        try:
            saveObject(filename=export_path, obj=self.project)
            self.status = ('Status: Project "%s" successfully exported' % self.project.name)
        except:
            self.status = ('Status: Attempt to export Project failed')
        self.refreshVisibility()

    def refreshStatus(self):
        if not self.guard_status:
            self.statusBar().showMessage(self.status)

    def refreshVisibility(self):
        self.refreshStatus()

        is_cascade_loaded = self.project is not None
        is_parset_loaded = is_cascade_loaded and len(self.project.parsets.keys()) > 0
        are_results_generated = (self.results_current is not None and self.results_comparison is not None)
        
        if is_parset_loaded:
            self.refreshParsetComboBoxes()
        self.label_parset.setVisible(is_parset_loaded)
        self.combo_parset.setVisible(is_parset_loaded)
        
        policy_min = qtw.QSizePolicy.Minimum
        policy_exp = qtw.QSizePolicy.Expanding
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
        
        if are_results_generated:
            self.refreshCharacComboBox()
            self.refreshPopComboBox()
        self.label_plotter_charac.setVisible(are_results_generated)
        self.combo_plotter_charac.setVisible(are_results_generated)
        self.label_plotter_pop.setVisible(are_results_generated)
        self.combo_plotter_pop.setVisible(are_results_generated)
        self.button_plotter.setVisible(are_results_generated)
            
    
    def refreshParsetComboBoxes(self):
        combo_boxes = [self.combo_parset, self.combo_compare]
        combo_names = [self.parset_source_name, self.parset_comparison_name]
        self.combo_parset_dict = {}
        for k in xrange(len(combo_boxes)):
            combo_box = combo_boxes[k]
            combo_name = combo_names[k]
            combo_box.clear()
            cid = 0
            for parset_name in self.project.parsets.keys():
                self.combo_parset_dict[parset_name] = cid
                combo_box.addItem(parset_name)
                cid += 1
            combo_box.setCurrentIndex(self.combo_parset_dict[combo_name])
            
    def refreshCharacComboBox(self):
        self.combo_charac_dict = {}
        self.combo_plotter_charac.clear()
        cid = 0
        for charac_label in self.project.settings.charac_specs.keys():
            charac_name = self.project.settings.charac_specs[charac_label]['name']
            self.combo_charac_dict[charac_name] = cid
            self.combo_plotter_charac.addItem(charac_name)
            cid += 1
        if self.charac_plot_name is None:
            self.charac_plot_name = self.combo_plotter_charac.itemText(self.combo_plotter_charac.currentIndex())
        self.combo_plotter_charac.setCurrentIndex(self.combo_charac_dict[self.charac_plot_name])
        
    def refreshPopComboBox(self):
        self.combo_pop_dict = {}
        self.combo_plotter_pop.clear()
        pid = 0
        for pop_name in self.parset.pop_names:
            self.combo_pop_dict[pop_name] = pid
            self.combo_plotter_pop.addItem(pop_name)
            pid += 1
        if self.pop_plot_name is None:
            self.pop_plot_name = self.combo_plotter_pop.itemText(self.combo_plotter_pop.currentIndex())
        self.combo_plotter_pop.setCurrentIndex(self.combo_pop_dict[self.pop_plot_name])
        
    def loadCalibration(self, parset_name, delay_refresh = False):
        self.parset_source_name = parset_name
        self.parset = dcp(self.project.parsets[parset_name])
        self.status = ('Status: Parset "%s" selected for editing' % parset_name)
        self.results_current = None
        self.results_comparison = None
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
        
    def selectCharacteristic(self, charac_name):
        self.charac_plot_name = charac_name
        
    def selectPopulation(self, pop_name):
        self.pop_plot_name = pop_name
        
    def runComparison(self):
        self.status = ('Status: Running models for Parset comparison')
        self.refreshStatus()
        self.results_current = self.project.runSim(parset = self.parset)
        self.status = ('Status: Model successfully processed for current Parset')
        self.refreshStatus()
        self.results_comparison = self.project.runSim(parset_name = self.parset_comparison_name)
        self.status = ('Status: Model successfully processed for Parset "%s"' % self.parset_comparison_name)
        self.refreshStatus()
        self.refreshVisibility()
        return
        
    def makeParsetTable(self):
        self.table_calibration.setVisible(False)    # Resizing columns requires table to be hidden first.
        self.table_calibration.clear()
        
        # Disconnect the calibration table from cell change signals to avoid signal flooding during connection.
        try: self.table_calibration.cellChanged.disconnect()
        except: pass

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
                        temp = qtw.QTableWidgetItem()
                        temp.setText(par_name)
                        temp.setFlags(qtc.Qt.ItemIsEnabled or qtc.Qt.ItemIsSelectable)
                        self.table_calibration.setItem(k*num_pops+pid, self.col_par_name, temp)
                        temp = qtw.QTableWidgetItem()
                        temp.setText(parset.pop_names[pid])
                        temp.setFlags(qtc.Qt.ItemIsEnabled or qtc.Qt.ItemIsSelectable)
                        self.table_calibration.setItem(k*num_pops+pid, self.col_pop_name, temp)
                        
                        for eid in xrange(len(par.t[pid])):
                            t = par.t[pid][eid]
                            y = par.y[pid][eid]
                            temp = qtw.QTableWidgetItem()
                            temp.setText(str(y))
                            self.table_calibration.setItem(k*num_pops+pid, 2+int(t)-self.tvec[0], temp)
                    k += 1
        self.table_calibration.setVerticalHeaderLabels(par_labels)
        self.table_calibration.setHorizontalHeaderLabels(['Par. Name','Pop. Name']+[str(int(x)) for x in self.tvec])
        self.table_calibration.resizeColumnsToContents()
        
        self.table_calibration.cellChanged.connect(self.updateParset)
        
    def showPlots(self):
        # Clear plots.
        for i in reversed(range(self.plotter_layout.count())): 
            self.plotter_layout.itemAt(i).widget().setParent(None)        
        
        if self.charac_plot_name is not None and self.pop_plot_name is not None:
            charac_plot_label = self.project.settings.charac_name_labels[self.charac_plot_name]
            pop_plot_label = self.project.data['pops']['name_labels'][self.pop_plot_name]
#            pid = self.combo_pop_dict[self.pop_plot_name]
            
            y_values_cur, t_values_cur = self.results_current.getValuesAt(label = charac_plot_label, year_init = self.tvec[0], year_end = self.tvec[-1], pop_labels = [pop_plot_label])
            y_values_com, t_values_com = self.results_comparison.getValuesAt(label = charac_plot_label, year_init = self.tvec[0], year_end = self.tvec[-1], pop_labels = [pop_plot_label])

#            y_values_cur, t_values_cur, final_dict_cur = extractCharacteristic(results=self.results_current, charac_label=charac_plot_label, charac_specs=self.project.settings.charac_specs, data=self.project.data)
#            y_values_com, t_values_com, final_dict_com = extractCharacteristic(results=self.results_comparison, charac_label=charac_plot_label, charac_specs=self.project.settings.charac_specs, data=self.project.data)

            print y_values_cur
            print t_values_cur
            print y_values_com
            print t_values_com

            try: y_data = self.project.data['characs'][charac_plot_label][pop_plot_label]['y']
            except: y_data = []
            try: t_data = self.project.data['characs'][charac_plot_label][pop_plot_label]['t']
            except: t_data = []

            print y_data
            print t_data

            figure = _plotLine(ys = [y_values_cur,y_values_com], ts = [t_values_cur,t_values_com], labels = ['Edited "%s"' % self.parset.name,'Unedited "%s"' % self.parset_comparison_name], y_hat=[y_data,y_data], t_hat=[t_data,t_data])
            
            canvas = FigureCanvasQTAgg(figure)

            self.plotter_layout.addWidget(canvas)        
        
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
#                        temp = qtw.QTableWidgetItem()
#                        temp.setText(par_name)
#                        temp.setFlags(qtc.Qt.ItemIsEnabled or qtc.Qt.ItemIsSelectable)
#                        self.table_calibration.setItem(k*num_pops+pid, self.col_par_name, temp)
#                        temp = qtw.QTableWidgetItem()
#                        temp.setText(parset.pop_names[pid])
#                        temp.setFlags(qtc.Qt.ItemIsEnabled or qtc.Qt.ItemIsSelectable)
#                        self.table_calibration.setItem(k*num_pops+pid, self.col_pop_name, temp)
#                        
#                        for eid in xrange(len(par.t[pid])):
#                            t = par.t[pid][eid]
#                            y = par.y[pid][eid]
#                            temp = qtw.QTableWidgetItem()
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
        self.status = ('Status: Current edited Parset uses value "%f" for parameter "%s", population "%s", year "%i"' % (new_val, par_label, pop_label, year))
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
    
    app = qtw.QApplication(sys.argv)
    app.setApplicationName('Optima GUI')
    gui = GUI()
    sys.exit(app.exec_())