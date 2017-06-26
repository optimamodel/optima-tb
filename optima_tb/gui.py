# %% Imports
import logging
import logging.config

logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

from matplotlib import pyplot as pp
pp.ioff()   # Turn off interactive mode.

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import sys
import numpy as np
from copy import deepcopy as dcp

from optima_tb.project import Project
from optima_tb.plotting import _plotLine
from optima_tb.dataio import saveObject, loadObject
from optima_tb.defaults import defaultOptimOptions

# %% PyQt imports

def importPyQt():
    ''' Try to import pyqt, either PyQt4 or PyQt5, but allow it to fail '''
    try:
        from PyQt5 import QtCore as qtc
        from PyQt5 import QtWidgets as qtw
    except:
        try:
            from PyQt4 import QtGui as qtw
            from PyQt4 import QtCore as qtc
        except Exception as E:
            errormsg = 'PyQt could not be imported: %s' % E.__repr__()
            raise Exception(errormsg)
    return  qtc, qtw

qtc, qtw = importPyQt()

# %% GUI classes

class GUI(qtw.QWidget):

    def __init__(self):
        super(GUI, self).__init__()
        self.initUI()

    def initUI(self):

        self.setWindowTitle('GUI Selection Screen')

        self.button_calibration = qtw.QPushButton('Manual Calibration', self)
        self.button_calibration.clicked.connect(self.runGUICalibration)
        self.button_scenario_parameter = qtw.QPushButton('Parameter Scenario', self)
        self.button_scenario_parameter.clicked.connect(self.runGUIScenarioParameter)
        self.button_scenario_budget = qtw.QPushButton('Budget Scenario', self)
        self.button_scenario_budget.clicked.connect(self.runGUIScenarioBudget)
        layout = qtw.QVBoxLayout(self)
        layout.addWidget(self.button_calibration)
        layout.addWidget(self.button_scenario_parameter)
        layout.addWidget(self.button_scenario_budget)
        self.setLayout(layout)

        screen = qtw.QDesktopWidget().availableGeometry()
        self.setGeometry((screen.width() - self.width()) / 2, (screen.height() - self.height()) / 2,
                         self.width(), self.height())

        self.show()

    def runGUICalibration(self):
        self.sub_gui = GUICalibration()

    def runGUIScenarioParameter(self):
        self.sub_gui = GUIResultPlotterIntermediate()

    def runGUIScenarioBudget(self):
        self.sub_gui = GUIBudgetScenario()


# A base GUI class that enables project creation, import and export.
class GUIProjectManagerBase(qtw.QMainWindow):

    def __init__(self):
        super(GUIProjectManagerBase, self).__init__()
        self.initUIProjectManager()


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to initialise all attributes used by this GUI.
    def resetAttributes(self):  self.resetAttributesProjectManager()
    def resetAttributesProjectManager(self):
        self.project = None     # This is the Project object that the GUI stores to calibrate and edit.
        self.tvec = None

        self.guard_status = False   # Convenient flag that locks status updates if set to true.
        self.status = 'Status: No Project generated'    # Initial status.


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to update widgets inside the GUI and what the user sees, as required.
    def refreshVisibility(self):  self.refreshVisibilityProjectManager()
    def refreshVisibilityProjectManager(self):
        self.refreshStatus()


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to acknowledge the existence of a project and run initial processes required by derived GUI.
    def acknowledgeProject(self): self.acknowledgeProjectProjectManager()
    def acknowledgeProjectProjectManager(self): pass


    def initUIProjectManager(self):
        self.resetAttributesProjectManager()    # This call is specific to base class attributes.

        self.setWindowTitle('Project Manager')

        # Screen.
        screen = qtw.QDesktopWidget().availableGeometry()
        self.resize(screen.width() * 9.0 / 10.0, screen.height() * 9.0 / 10.0)
        self.setGeometry((screen.width() - self.width()) / 2, (screen.height() - self.height()) / 2,
                         self.width(), self.height())

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

        self.refreshVisibilityProjectManager()      # This call is specific to base class components.
        self.show()

    def createProject(self):
        self.resetAttributes()      # This call is to a derived method, representing a full GUI reset.
        project_name, project_name_chosen = qtw.QInputDialog.getText(self, 'Create Project', 'Enter project name:')
        if project_name_chosen:
            try: cascade_path = str(qtw.QFileDialog.getOpenFileNameAndFilter(self, 'Select Cascade File')[0])
            except: cascade_path = str(qtw.QFileDialog.getOpenFileName(self, 'Select Cascade File')[0])
            try:
                self.project = Project(name=project_name, cascade_path=cascade_path, validation_level='avert')
                self.tvec = np.arange(self.project.settings.tvec_start, self.project.settings.tvec_observed_end + 1.0 / 2)
                self.status = ('Status: Project "%s" generated, cascade settings loaded' % self.project.name)
            except Exception as E:
                self.status = ('Status: Attempt to generate Project failed (%s)' % E.__repr__())
                self.refreshVisibility()
                return

            try: databook_path = str(qtw.QFileDialog.getOpenFileNameAndFilter(self, 'Select Databook File')[0])
            except: databook_path = str(qtw.QFileDialog.getOpenFileName(self, 'Select Databook File')[0])
            try:
                self.project.loadSpreadsheet(databook_path=databook_path)
                self.project.resetParsets()
                self.acknowledgeProject()
                self.status = ('Status: Valid data loaded into Project "%s", default Parset generated' % self.project.name)
            except Exception as E:
                self.resetAttributes()      # This call is to a derived method, representing a full GUI reset.
                self.status = ('Status: Attempt to load data into Project failed, Project reset for safety (%s)' % E.__repr__())
        self.refreshVisibility()

    def importProject(self):
        self.resetAttributes()      # This call is to a derived method, representing a full GUI reset.
        try: import_path = str(qtw.QFileDialog.getOpenFileNameAndFilter(self, 'Import Project From File')[0])
        except: import_path = str(qtw.QFileDialog.getOpenFileName(self, 'Import Project From File')[0])
        try:
            self.project = loadObject(filename=import_path)
            self.tvec = np.arange(self.project.settings.tvec_start, self.project.settings.tvec_observed_end + 1.0 / 2)
            self.acknowledgeProject()
            self.status = ('Status: Project "%s" successfully imported' % self.project.name)
        except:
            self.status = ('Status: Attempt to import Project failed')
        self.refreshVisibility()

    def exportProject(self):
        try: export_path = str(qtw.QFileDialog.getSaveFileNameAndFilter(self, 'Export Project To File')[0])
        except: export_path = str(qtw.QFileDialog.getSaveFileName(self, 'Export Project To File')[0])
        try:
            saveObject(filename=export_path, obj=self.project)
            self.status = ('Status: Project "%s" successfully exported' % self.project.name)
        except:
            self.status = ('Status: Attempt to export Project failed')
        self.refreshVisibility()

    def refreshStatus(self):
        if not self.guard_status:
            self.statusBar().showMessage(self.status)


# An intermediate GUI class that extends project management and allows for result comparison.
class GUIResultPlotterIntermediate(GUIProjectManagerBase):

    def __init__(self):
        super(GUIResultPlotterIntermediate, self).__init__()
        self.initUIResultPlotter()


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to initialise all attributes used by this GUI.
    def resetAttributes(self):
        self.resetAttributesProjectManager()
        self.resetAttributesResultPlotter()
    def resetAttributesResultPlotter(self):
        self.result_1_plot_name = None
        self.result_2_plot_name = None
        self.charac_plot_name = None
        self.pop_plot_name = None
        self.combo_result_dict = {}     # Dictionary that maps Result names to indices used in combo boxes.
        self.combo_charac_dict = {}     # Dictionary that maps Characteristic names to indices used in combo boxes.
        self.combo_pop_dict = {}        # Dictionary that maps Population names to indices used in combo boxes.


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to update widgets inside the GUI and what the user sees, as required.
    def refreshVisibility(self):
        self.refreshVisibilityProjectManager()
        self.refreshVisibilityResultPlotter()
    def refreshVisibilityResultPlotter(self):
        is_project_loaded = self.project is not None
        does_project_have_results = is_project_loaded and len(self.project.results) > 0

        if does_project_have_results:
            self.refreshResultComboBoxes()
            self.refreshCharacComboBox()
            self.refreshPopComboBox()
        self.label_plotter_result_1.setVisible(does_project_have_results)
        self.combo_plotter_result_1.setVisible(does_project_have_results)
        self.label_plotter_result_2.setVisible(does_project_have_results)
        self.combo_plotter_result_2.setVisible(does_project_have_results)
        self.label_plotter_charac.setVisible(does_project_have_results)
        self.combo_plotter_charac.setVisible(does_project_have_results)
        self.label_plotter_pop.setVisible(does_project_have_results)
        self.combo_plotter_pop.setVisible(does_project_have_results)
        self.button_plotter.setVisible(does_project_have_results)


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to acknowledge the existence of a project and run initial processes required by derived GUI.
    def acknowledgeProject(self):
        self.acknowledgeProjectProjectManager()
        self.acknowledgeProjectResultPlotter()
    def acknowledgeProjectResultPlotter(self):
        if self.project is not None:
            self.acknowledgeResults()


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to extend the user interface ahead of result plotting.
    def developLayout(self, layout): self.developLayoutResultPlotter()
    def developLayoutResultPlotter(self):
        pass


    def initUIResultPlotter(self):
        self.resetAttributesResultPlotter()    # This call is specific to base class attributes.

        self.setWindowTitle('Result Plotter')

        # Widgets.
        self.label_plotter_result_1 = qtw.QLabel('Select First Result: ')
        self.combo_plotter_result_1 = qtw.QComboBox(self)
        self.combo_plotter_result_1.activated[str].connect(self.selectResultOne)
        self.label_plotter_result_2 = qtw.QLabel('Select Second Result: ')
        self.combo_plotter_result_2 = qtw.QComboBox(self)
        self.combo_plotter_result_2.activated[str].connect(self.selectResultTwo)
        self.label_plotter_charac = qtw.QLabel('Select Characteristic: ')
        self.combo_plotter_charac = qtw.QComboBox(self)
        self.combo_plotter_charac.activated[str].connect(self.selectCharacteristic)
        self.label_plotter_pop = qtw.QLabel('Select Population: ')
        self.combo_plotter_pop = qtw.QComboBox(self)
        self.combo_plotter_pop.activated[str].connect(self.selectPopulation)
        self.button_plotter = qtw.QPushButton('Plot Figures', self)
        self.button_plotter.clicked.connect(self.showPlots)

        # Layout.
        grid_lower = qtw.QGridLayout()
        grid_lower.setSpacing(10)

        grid_lower.addWidget(self.label_plotter_result_1, 0, 0)
        grid_lower.addWidget(self.combo_plotter_result_1, 0, 1)
        grid_lower.addWidget(self.label_plotter_result_2, 1, 0)
        grid_lower.addWidget(self.combo_plotter_result_2, 1, 1)
        grid_lower.addWidget(self.label_plotter_charac, 2, 0)
        grid_lower.addWidget(self.combo_plotter_charac, 2, 1)
        grid_lower.addWidget(self.label_plotter_pop, 3, 0)
        grid_lower.addWidget(self.combo_plotter_pop, 3, 1)
        grid_lower.addWidget(self.button_plotter, 3, 2)

        policy_min = qtw.QSizePolicy.Minimum
        policy_exp = qtw.QSizePolicy.Expanding
        self.process_layout_stretch = qtw.QSpacerItem(0, 0, policy_min, policy_exp)

        process_layout = qtw.QVBoxLayout()
        self.developLayout(process_layout)
        process_layout.addLayout(grid_lower)
        process_layout.addItem(self.process_layout_stretch)

        self.plotter_layout = qtw.QGridLayout()
        self.plotter_layout.setSpacing(5)

        self.first_half = qtw.QWidget()
        self.first_half.resize(self.width() / 2.0, self.height())
        self.first_half.setLayout(process_layout)
        self.second_half = qtw.QWidget()
        self.second_half.resize(self.width() / 2.0, self.height() / 2.0)
        self.second_half.setLayout(self.plotter_layout)

        # Window splitter.
        self.splitter_total = qtw.QSplitter()
        self.splitter_total.addWidget(self.first_half)
        self.splitter_total.addWidget(self.second_half)
        self.setCentralWidget(self.splitter_total)

        self.refreshVisibilityResultPlotter()      # This call is specific to base class components.
        self.show()

    def refreshResultComboBoxes(self):
        combo_boxes = [self.combo_plotter_result_1, self.combo_plotter_result_2]
        combo_names = [self.result_1_plot_name, self.result_2_plot_name]
        self.combo_result_dict = {}
        for k in xrange(len(combo_boxes)):
            combo_box = combo_boxes[k]
            combo_name = combo_names[k]
            combo_box.clear()
            rid = 0
            for result_name in self.project.results.keys():
                self.combo_result_dict[result_name] = rid
                combo_box.addItem(result_name)
                rid += 1
            try: combo_box.setCurrentIndex(self.combo_result_dict[combo_name])
            except: pass    # Should be triggered if there are no results.

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
            self.charac_plot_name = str(self.combo_plotter_charac.itemText(self.combo_plotter_charac.currentIndex()))
        self.combo_plotter_charac.setCurrentIndex(self.combo_charac_dict[self.charac_plot_name])

    def refreshPopComboBox(self):
        self.combo_pop_dict = {}
        self.combo_plotter_pop.clear()
        pid = 0
        for pop_name in self.project.data['pops']['name_labels'].keys():
            self.combo_pop_dict[pop_name] = pid
            self.combo_plotter_pop.addItem(pop_name)
            pid += 1
        if self.pop_plot_name is None:
            self.pop_plot_name = str(self.combo_plotter_pop.itemText(self.combo_plotter_pop.currentIndex()))
        self.combo_plotter_pop.setCurrentIndex(self.combo_pop_dict[self.pop_plot_name])

    def selectResultOne(self, result_name):
        self.result_1_plot_name = str(result_name)

    def selectResultTwo(self, result_name):
        self.result_2_plot_name = str(result_name)

    def selectCharacteristic(self, charac_name):
        self.charac_plot_name = str(charac_name)

    def selectPopulation(self, pop_name):
        self.pop_plot_name = str(pop_name)
        
    # When a calibration or scenario is run, results should be acknowledged.
    # This means that result name references for the plotter should be linked to the first and last results in the project.
    def acknowledgeResults(self):
        if len(self.project.results) > 0:
                self.result_1_plot_name = self.project.results.keys()[0]
                self.result_2_plot_name = self.project.results.keys()[-1]

    def showPlots(self):
        # Clear plots.
        for i in reversed(range(self.plotter_layout.count())):
            self.plotter_layout.itemAt(i).widget().setParent(None)

        if self.result_1_plot_name is None or self.result_2_plot_name is None:
            self.status = ('Status: No results could be selected for plotting purposes')
            self.refreshVisibility()
            return

        if self.charac_plot_name is not None and self.pop_plot_name is not None:

            charac_plot_label = self.project.settings.charac_name_labels[self.charac_plot_name]
            pop_plot_label = self.project.data['pops']['name_labels'][self.pop_plot_name]

            y_values_cur, t_values_cur = self.project.results[self.result_1_plot_name].getValuesAt(label=charac_plot_label, year_init=self.tvec[0], year_end=self.tvec[-1], pop_labels=[pop_plot_label])
            y_values_com, t_values_com = self.project.results[self.result_2_plot_name].getValuesAt(label=charac_plot_label, year_init=self.tvec[0], year_end=self.tvec[-1], pop_labels=[pop_plot_label])

            try: y_data = self.project.data['characs'][charac_plot_label][pop_plot_label]['y']
            except: y_data = []
            try: t_data = self.project.data['characs'][charac_plot_label][pop_plot_label]['t']
            except: t_data = []

            figure = _plotLine(ys=[y_values_cur, y_values_com], ts=[t_values_cur, t_values_com], labels=['%s' % self.result_1_plot_name, '%s' % self.result_2_plot_name], y_hat=[y_data, y_data], t_hat=[t_data, t_data])

            canvas = FigureCanvasQTAgg(figure)

            self.plotter_layout.addWidget(canvas)




class GUICalibration(GUIResultPlotterIntermediate):

    def __init__(self):
        super(GUICalibration, self).__init__()
        self.initUICalibration()


    def resetAttributes(self):
        self.resetAttributesProjectManager()
        self.resetAttributesResultPlotter()

        self.parset = None      # This is the ParameterSet object that stores all edits made in the GUI.
                                # It should be an (edited) copy, not a reference, to an existing Project ParameterSet.
        self.parset_name = None
        self.combo_parset_dict = {}     # Dictionary that maps Parset names to indices used in combo boxes.
        self.col_par_name = 0   # Column index for table calibration parameter names.
        self.col_pop_name = 1   # Column index for table calibration population names.


    def refreshVisibility(self):
        self.refreshVisibilityProjectManager()
        self.refreshVisibilityResultPlotter()

        is_project_loaded = self.project is not None
        does_parset_exist = is_project_loaded and len(self.project.parsets.keys()) > 0

        if is_project_loaded:
            self.refreshParsetComboBox()
        self.label_parset.setVisible(is_project_loaded)
        self.combo_parset.setVisible(is_project_loaded)

        policy_min = qtw.QSizePolicy.Minimum
        policy_exp = qtw.QSizePolicy.Expanding
        if does_parset_exist:
            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_min)
            self.makeParsetTable()
        else:
            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_exp)
        self.table_calibration.setVisible(does_parset_exist)

        self.label_model_run.setVisible(does_parset_exist)
        self.edit_model_run.setVisible(does_parset_exist)
        self.button_model_run.setVisible(does_parset_exist)
        self.label_overwrite.setVisible(does_parset_exist)
        self.edit_overwrite.setVisible(does_parset_exist)
        self.button_overwrite.setVisible(does_parset_exist)


    def acknowledgeProject(self):
        self.acknowledgeProjectProjectManager()
        self.acknowledgeProjectResultPlotter()
        if not self.project is None:
            if len(self.project.parsets) < 1:
                self.project.makeParset(name='default')
            self.loadCalibration(self.project.parsets[0].name, delay_refresh=True)


    # While UI initialisation can extend the interface, this method is where widgets for the core process should be set up.
    def developLayout(self, layout):
        self.developLayoutResultPlotter()

        # Widgets.
        self.label_parset = qtw.QLabel('Parset To Edit: ')
        self.combo_parset = qtw.QComboBox(self)
        self.combo_parset.activated[str].connect(self.loadCalibration)

        self.table_calibration = qtw.QTableWidget()
        self.table_calibration.cellChanged.connect(self.updateParset)

        self.label_model_run = qtw.QLabel('Run & Save Calibration Results As... ')
        self.edit_model_run = qtw.QLineEdit()
        self.button_model_run = qtw.QPushButton('Generate Results', self)
        self.button_model_run.clicked.connect(self.runCalibration)

        self.label_overwrite = qtw.QLabel('Save Edits To... ')
        self.edit_overwrite = qtw.QLineEdit()
        self.button_overwrite = qtw.QPushButton('Save Calibration', self)
        self.button_overwrite.clicked.connect(self.saveCalibration)

        # Layout.
        grid_parset_load = qtw.QGridLayout()
        grid_parset_load.setSpacing(10)

        grid_parset_load.addWidget(self.label_parset, 0, 0)
        grid_parset_load.addWidget(self.combo_parset, 0, 1)

        grid_parset_save = qtw.QGridLayout()
        grid_parset_save.setSpacing(10)

        grid_parset_save.addWidget(self.label_overwrite, 0, 0)
        grid_parset_save.addWidget(self.edit_overwrite, 0, 1)
        grid_parset_save.addWidget(self.button_overwrite, 0, 2)
        grid_parset_save.addWidget(self.label_model_run, 1, 0)
        grid_parset_save.addWidget(self.edit_model_run, 1, 1)
        grid_parset_save.addWidget(self.button_model_run, 1, 2)

        layout.addLayout(grid_parset_load)
        layout.addWidget(self.table_calibration)
        layout.addLayout(grid_parset_save)


    def initUICalibration(self):
        self.resetAttributes()

        self.setWindowTitle('Manual Calibration')

        self.refreshVisibility()
        self.show()


    def refreshParsetComboBox(self):
        self.combo_parset_dict = {}
        self.combo_parset.clear()
        pid = 0
        for parset_name in self.project.parsets.keys():
            self.combo_parset_dict[parset_name] = pid
            self.combo_parset.addItem(parset_name)
            pid += 1
        try: self.combo_parset.setCurrentIndex(self.combo_parset_dict[self.parset_name])
        except: pass


    def loadCalibration(self, parset_name, delay_refresh=False):
        self.parset_name = str(parset_name)
        self.parset = dcp(self.project.parsets[self.parset_name])
        self.status = ('Status: Parset "%s" selected for editing' % self.parset_name)
        if not delay_refresh:
            self.refreshVisibility()

    def runCalibration(self):
        self.status = ('Status: Running model for Parset "%s"' % self.parset_name)
        self.refreshStatus()
        result_name = str(self.edit_model_run.text())
        if result_name == '':
            result_name = None
        self.project.runSim(parset=self.parset, store_results=True, result_type='calibration', result_name=result_name)
        self.acknowledgeResults()
        self.status = ('Status: Model successfully processed for Parset "%s"' % self.parset_name)
        self.refreshVisibility()

    def saveCalibration(self):
        parset_name = str(self.edit_overwrite.text())
        if parset_name == '':
            self.status = ('Status: Attempt to save Parset failed, no name provided')
        else:
            if parset_name in self.project.parsets.keys():
                self.status = ('Status: Parset "%s" successfully overwritten' % parset_name)
            else:
                self.status = ('Status: New Parset "%s" added to Project' % parset_name)
            self.parset.name = parset_name
            self.project.parsets[parset_name] = dcp(self.parset)
        self.refreshVisibility()


    def makeParsetTable(self):
        self.table_calibration.setVisible(False)    # Resizing columns requires table to be hidden first.
        self.table_calibration.clear()

        # Disconnect the calibration table from cell change signals to avoid signal flooding during connection.
        try: self.table_calibration.cellChanged.disconnect()
        except: pass

        parset = self.parset
        num_pops = len(parset.pop_labels)
        row_count = num_pops * (len(parset.pars['cascade']) - len(self.project.settings.par_funcs))
        self.table_calibration.setRowCount(row_count)
        self.table_calibration.setColumnCount(2 + len(self.tvec))
        self.calibration_items = []

        k = 0
        par_labels = []
        for par_type in ['characs', 'cascade']:
            for par in parset.pars[par_type]:
                if ((par_type == 'cascade' and par.label not in self.project.settings.par_funcs.keys()) or (par_type == 'characs' and 'entry_point' in self.project.settings.charac_specs[par.label].keys())):
                    for pid in xrange(len(parset.pop_labels)):
                        pop_label = parset.pop_labels[pid]
                        par_labels.append(par.label + ' [' + pop_label + ']')
                        try:
                            par_name = self.project.settings.linkpar_specs[par.label]['name']
                        except:
                            par_name = self.project.settings.charac_specs[par.label]['name']
                        temp = qtw.QTableWidgetItem()
                        temp.setText(par_name)
                        temp.setFlags(qtc.Qt.ItemIsEnabled or qtc.Qt.ItemIsSelectable)
                        self.table_calibration.setItem(k * num_pops + pid, self.col_par_name, temp)
                        temp = qtw.QTableWidgetItem()
                        temp.setText(parset.pop_names[pid])
                        temp.setFlags(qtc.Qt.ItemIsEnabled or qtc.Qt.ItemIsSelectable)
                        self.table_calibration.setItem(k * num_pops + pid, self.col_pop_name, temp)

                        for eid in xrange(len(par.t[pid])):
                            t = par.t[pid][eid]
                            y = par.y[pid][eid]
                            temp = qtw.QTableWidgetItem()
                            temp.setText(str(y))
                            self.table_calibration.setItem(k * num_pops + pid, 2 + int(t) - self.tvec[0], temp)
                    k += 1
        self.table_calibration.setVerticalHeaderLabels(par_labels)
        self.table_calibration.setHorizontalHeaderLabels(['Par. Name', 'Pop. Name'] + [str(int(x)) for x in self.tvec])
        self.table_calibration.resizeColumnsToContents()

        self.table_calibration.cellChanged.connect(self.updateParset)


    def updateParset(self, row, col):
        new_val_str = str(self.table_calibration.item(row, col).text())
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
                remove_success = par.removeValueAt(t=year, pop_label=pop_label)
                if remove_success:
                    self.status = ('Status: Value successfully deleted from Parset')
                    self.refreshStatus()
                    return
                else: self.status = ('Status: Attempt to remove item in Parset failed, at least one value per row required')
            else: self.status = ('Status: Attempt to edit item in Parset failed, only numbers allowed')
            self.refreshStatus()
            self.guard_status = True
            self.table_calibration.item(row, col).setText(str(par.interpolate(tvec=[year], pop_label=pop_label)[0]))
            self.guard_status = False
            return
        par.insertValuePair(t=year, y=new_val, pop_label=pop_label)
        self.status = ('Status: Current edited Parset uses value "%f" for parameter "%s", population "%s", year "%i"' % (new_val, par_label, pop_label, year))
        self.refreshStatus()



        
class GUIParameterScenario(GUIResultPlotterIntermediate):

    def __init__(self):
        super(GUIParameterScenario, self).__init__()
        self.initUIParameterScenario()


    def resetAttributes(self):
        self.resetAttributesProjectManager()
        self.resetAttributesResultPlotter()

        # Initialise attributes specific to your parameter scenario GUI.


    def refreshVisibility(self):
        self.refreshVisibilityProjectManager()
        self.refreshVisibilityResultPlotter()

        # Update the visibility of widgets depending on if they have anything to show.


    def acknowledgeProject(self):
        self.acknowledgeProjectProjectManager()
        self.acknowledgeProjectResultPlotter()
        
        # If a project is loaded, do whatever initial pre-processing is needed.


    # While UI initialisation can extend the interface, this method is where widgets for the core process should be set up.
    def developLayout(self, layout):
        self.developLayoutResultPlotter()

        # Initialise all widgets specific to the GUI process.
        # They will appear on the left side of the screen beneath the file menu and above the project plotter.


    def initUICalibration(self):
        self.resetAttributes()

        self.setWindowTitle('Parameter Scenario')

        self.refreshVisibility()
        self.show()


    # Whatever auxiliary functions you require to link to widgets and so on.
    # E.g...
    def translateToParameterScenario(self, param_diag, param_linkage, param_treatment_succ, param_treatment_fail):
        """
        
        using the returned scvalues, correct usage should be
        
        proj.createScenario(scvalues)
        resultsset = proj.runScenarios()
        
        """
        return None # scvalues odict
        
    
# DJK to CK: Do this one after the Parameter Scenario as I -might- get to it before you.
class GUIBudgetScenario(GUIResultPlotterIntermediate):

    def __init__(self):
        super(GUIBudgetScenario, self).__init__()
        self.initUIBudgetScenario()


    def resetAttributes(self):
        self.resetAttributesProjectManager()
        self.resetAttributesResultPlotter()
        
        self.parset_name = None
        self.progset_name = None
        self.combo_parset_dict = {}     # Dictionary that maps Parset names to indices used in combo boxes.
        self.combo_progset_dict = {}    # Dictionary that maps Progset names to indices used in combo boxes.
        
        self.options = None     # The options dictionary for running a budget scenario.

        # Initialise attributes specific to your budget scenario GUI.


    def refreshVisibility(self):
        self.refreshVisibilityProjectManager()
        self.refreshVisibilityResultPlotter()
        
        is_project_loaded = self.project is not None
        does_parset_exist = is_project_loaded and len(self.project.parsets.keys()) > 0
        does_progset_exist = is_project_loaded and len(self.project.progsets.keys()) > 0

        if is_project_loaded:
            self.refreshParsetComboBox()
            self.refreshProgsetComboBox()
        self.label_parset.setVisible(is_project_loaded)
        self.combo_parset.setVisible(is_project_loaded)
        self.label_progset.setVisible(is_project_loaded)
        self.combo_progset.setVisible(is_project_loaded)

#        policy_min = qtw.QSizePolicy.Minimum
#        policy_exp = qtw.QSizePolicy.Expanding
#        if does_parset_exist:
#            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_min)
#            self.makeParsetTable()
#        else:
#            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_exp)
#        self.table_calibration.setVisible(does_parset_exist)

        self.label_model_run.setVisible(does_parset_exist & does_progset_exist)
        self.edit_model_run.setVisible(does_parset_exist & does_progset_exist)
        self.button_model_run.setVisible(does_parset_exist & does_progset_exist)

        # Update the visibility of widgets depending on if they have anything to show.


    def acknowledgeProject(self):
        self.acknowledgeProjectProjectManager()
        self.acknowledgeProjectResultPlotter()
        if not self.project is None:
            if len(self.project.parsets) < 1:
                self.project.makeParset(name='default')
            self.loadCalibration(self.project.parsets[0].name, delay_refresh=True)
            if len(self.project.progsets) < 1:
                self.project.makeProgset(name='default')
            self.loadPrograms(self.project.progsets[0].name, delay_refresh=True)
        
        # If a project is loaded, do whatever initial pre-processing is needed.


    # While UI initialisation can extend the interface, this method is where widgets for the core process should be set up.
    def developLayout(self, layout):
        self.developLayoutResultPlotter()
        
        # Widgets.
        self.label_parset = qtw.QLabel('Parset To Use: ')
        self.combo_parset = qtw.QComboBox(self)
        self.combo_parset.activated[str].connect(self.loadCalibration)
        
        self.label_progset = qtw.QLabel('Progset To Use: ')
        self.combo_progset = qtw.QComboBox(self)
        self.combo_progset.activated[str].connect(self.loadPrograms)

        self.label_model_run = qtw.QLabel('Run & Save Budget Scenario Results As... ')
        self.edit_model_run = qtw.QLineEdit()
        self.button_model_run = qtw.QPushButton('Generate Results', self)
        self.button_model_run.clicked.connect(self.runBudgetScenario)

        # Layout.
        grid_progset_load = qtw.QGridLayout()
        grid_progset_load.setSpacing(10)

        grid_progset_load.addWidget(self.label_parset, 0, 0)
        grid_progset_load.addWidget(self.combo_parset, 0, 1)
        grid_progset_load.addWidget(self.label_progset, 1, 0)
        grid_progset_load.addWidget(self.combo_progset, 1, 1)

        grid_progset_save = qtw.QGridLayout()
        grid_progset_save.setSpacing(10)

        grid_progset_save.addWidget(self.label_model_run, 0, 0)
        grid_progset_save.addWidget(self.edit_model_run, 0, 1)
        grid_progset_save.addWidget(self.button_model_run, 0, 2)

        layout.addLayout(grid_progset_load)
#        layout.addWidget(self.table_calibration)
        layout.addLayout(grid_progset_save)


    def initUIBudgetScenario(self):
        self.resetAttributes()

        self.setWindowTitle('Budget Scenario')

        self.refreshVisibility()
        self.show()
        
        
    def refreshParsetComboBox(self):
        self.combo_parset_dict = {}
        self.combo_parset.clear()
        pid = 0
        for parset_name in self.project.parsets.keys():
            self.combo_parset_dict[parset_name] = pid
            self.combo_parset.addItem(parset_name)
            pid += 1
        try: self.combo_parset.setCurrentIndex(self.combo_parset_dict[self.parset_name])
        except: pass
    
    def refreshProgsetComboBox(self):
        self.combo_progset_dict = {}
        self.combo_progset.clear()
        pid = 0
        for progset_name in self.project.progsets.keys():
            self.combo_progset_dict[progset_name] = pid
            self.combo_progset.addItem(progset_name)
            pid += 1
        try: self.combo_progset.setCurrentIndex(self.combo_progset_dict[self.progset_name])
        except: pass

    def loadCalibration(self, parset_name, delay_refresh=False):
        self.parset_name = str(parset_name)
        self.parset = dcp(self.project.parsets[self.parset_name])
        self.status = ('Status: Parset "%s" selected for budget scenario' % self.parset_name)
        if not delay_refresh:
            self.refreshVisibility()
            
    def loadPrograms(self, progset_name, delay_refresh=False):
        self.progset_name = str(progset_name)
#        self.progset = dcp(self.project.progsets[progset_name])
        self.options = defaultOptimOptions(settings = self.project.settings, progset = self.project.progsets[self.progset_name])
        self.status = ('Status: Progset "%s" selected for budget scenario' % self.progset_name)
        if not delay_refresh:
            self.refreshVisibility()
            
    def runBudgetScenario(self):
        self.status = ('Status: Running model for Parset "%s" and Progset "%s"' % (self.parset_name, self.progset_name))
        self.refreshStatus()
        result_name = str(self.edit_model_run.text())
        if result_name == '':
            result_name = None
        self.project.runSim(parset_name=self.parset_name, progset_name=self.progset_name, options=self.options, store_results=True, result_type='scen_budget', result_name=result_name)
        self.acknowledgeResults()
        self.status = ('Status: Model successfully processed for Parset "%s" and Progset "%s"' % (self.parset_name, self.progset_name))
        self.refreshVisibility()



# %% Functionality for opening GUIs

def runGUI():
    ''' Function that launches all available back-end GUIs as they are developed. '''

    app = qtw.QApplication(sys.argv)
    app.setApplicationName('Optima GUI')
    gui = GUI()
    sys.exit(app.exec_())
