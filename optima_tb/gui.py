# %% Imports
import logging
import logging.config

logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

import sys
import numpy as np
from copy import deepcopy as dcp
from matplotlib import pyplot as pp
import pylab as pl
import itertools as it
pp.ioff()   # Turn off interactive mode.

from optima_tb.project import Project
from optima_tb.plotting import plotScenarios, plotCompareResults # _plotLine
from optima_tb.dataio import saveObject, loadObject
from optima_tb.defaults import defaultOptimOptions
from optima_tb.utils import odict
from optima_tb.settings import PlottingSettings, DO_NOT_SCALE, DEFAULT_YFACTOR

##### TODO remove hardcoded from Belarus:
colors = ['#8C8984',
          '#333399']
#### /hardcoded ref


# %% PyQt imports

def importPyQt():
    ''' Try to import pyqt, either PyQt4 or PyQt5, but allow it to fail '''
    try:
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
        from PyQt5 import QtCore as qtc
        from PyQt5 import QtWidgets as qtw
    except:
        try:
            from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
            from PyQt4 import QtGui as qtw
            from PyQt4 import QtCore as qtc
        except Exception as E:
            errormsg = 'PyQt could not be imported: %s' % E.__repr__()
            raise Exception(errormsg)
    return  qtc, qtw, FigureCanvasQTAgg

qtc, qtw, FigureCanvasQTAgg = importPyQt()

# %% Sanitized path import
def sanitizedFileDialog(instance=None, which=None, title=None):
    '''
    Get a sanitized path from a file dialog, either open or save. Examples:
    
    path = sanitizedFileDialog(self, 'open', 'Choose file to open')
    path = sanitizedFileDialog(self, 'save', 'Save file as')
    '''
    if which == 'open':
        try:    path = str(qtw.QFileDialog.getOpenFileNameAndFilter(instance, title)[0])
        except: path = str(qtw.QFileDialog.getOpenFileName(instance, title)[0])
    elif which == 'save':
        try:    path = str(qtw.QFileDialog.getSaveFileNameAndFilter(instance, title)[0])
        except: path = str(qtw.QFileDialog.getSaveFileName(instance, title)[0])
    else:
        raise Exception('The argument "which" must be either "open" or "save", not %s' % which)
    return path

# %% Wrapper GUI class

class GUI(qtw.QWidget):

    def __init__(self):
        super(GUI, self).__init__()
#         self.setStyleSheet("background-color:white;")
        self.initUI()

    def initUI(self):

        self.setWindowTitle('GUI selection screen')

        self.button_calibration = qtw.QPushButton('Edit parameters', self)
        self.button_calibration.clicked.connect(self.runGUICalibration)
        self.button_reconciliation = qtw.QPushButton('Edit programs', self)
        self.button_reconciliation.clicked.connect(self.runGUIReconciliation)
#        self.button_scenario_parameter = qtw.QPushButton('Parameter scenario', self)
#        self.button_scenario_parameter.clicked.connect(self.runGUIParameterScenario)
        self.button_scenario_budget = qtw.QPushButton('Budget scenario', self)
        self.button_scenario_budget.clicked.connect(self.runGUIBudgetScenario)
        layout = qtw.QVBoxLayout(self)
        layout.addWidget(self.button_calibration)
        layout.addWidget(self.button_reconciliation)
#        layout.addWidget(self.button_scenario_parameter)
        layout.addWidget(self.button_scenario_budget)
        self.setLayout(layout)

        screen = qtw.QDesktopWidget().availableGeometry()
        self.setGeometry((screen.width() - self.width()) / 2, (screen.height() - self.height()) / 2,
                         self.width(), self.height())

        self.show()

    def runGUICalibration(self):
        self.sub_gui = GUICalibration()
        
    def runGUIReconciliation(self):
        self.sub_gui = GUIReconciliation()

    def runGUIParameterScenario(self):
        self.sub_gui = GUIParameterScenario()

    def runGUIBudgetScenario(self):
        self.sub_gui = GUIBudgetScenario()

# %% GUI class for project management

# A base GUI class that enables project creation, import and export.
class GUIProjectManagerBase(qtw.QMainWindow):

    def __init__(self):
        super(GUIProjectManagerBase, self).__init__()
#         self.setStyleSheet("background-color:white;")
        self.initUIProjectManager()


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to initialise all attributes used by this GUI.
    def resetAttributes(self):  self.resetAttributesProjectManager()
    def resetAttributesProjectManager(self):
        self.project = None     # This is the Project object that the GUI stores to calibrate and edit.
        self.tvec = None

        self.guard_status = False   # Convenient flag that locks status updates if set to true.
        self.status = 'Status: No project generated'    # Initial status.


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to update widgets inside the GUI and what the user sees, as required.
    def refreshVisibility(self):
        self.refreshVisibilityProjectManager()
    def refreshVisibilityProjectManager(self):
        self.refreshStatus()


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to acknowledge the existence of a project and run initial processes required by derived GUI.
    def acknowledgeProject(self): self.acknowledgeProjectProjectManager()
    def acknowledgeProjectProjectManager(self): pass


    def initUIProjectManager(self):
        self.resetAttributesProjectManager()    # This call is specific to base class attributes.

        self.setWindowTitle('Project manager')

        # Screen.
        widget = qtw.QDesktopWidget()
#         p = widget.palette()
#         p.setColor(widget.backgroundRole(), Qt.red)
#         widget.setPalette(p)

        screen = widget.availableGeometry()
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
        project_name, project_name_chosen = qtw.QInputDialog.getText(self, 'Create project', 'Enter project name:')
        if project_name_chosen:
            cascade_path = sanitizedFileDialog(self, 'open', 'Select cascade file')
            try:
                self.project = Project(name=project_name, cascade_path=cascade_path, validation_level='avert', plotting_level='gui')
                self.tvec = np.arange(self.project.settings.tvec_start, self.project.settings.tvec_observed_end + 1.0 / 2)
                self.status = ('Status: Project "%s" generated, cascade settings loaded' % self.project.name)
            except Exception as E:
                self.status = ('Status: Attempt to generate project failed (%s)' % E.__repr__())
                self.refreshVisibility()
                return

            # set year ##### TODO remove hardcoded references
            plot_over = [2000, 2035]
            dt = 0.25
            self.project.setYear(plot_over, False)
            self.project.setYear(plot_over, True)
            self.project.settings.plot_settings['x_ticks'] = [np.arange(plot_over[0], plot_over[1] + dt, 5, dtype=int), np.arange(plot_over[0], plot_over[1] + dt, 5, dtype=int)]
            self.project.settings.plot_settings['xlim'] = (plot_over[0] - dt, plot_over[1] + dt)
            print self.project.settings.tvec_end
            #### /hardcoded references

            databook_path = sanitizedFileDialog(self, 'open', 'Select databook file')
            try:
                self.project.loadSpreadsheet(databook_path=databook_path)
                self.project.resetParsets()
                self.acknowledgeProject()
                self.status = ('Status: Valid data loaded into project "%s", default parameter set generated' % self.project.name)
            except Exception as E:
                self.resetAttributes()      # This call is to a derived method, representing a full GUI reset.
                self.status = ('Status: Attempt to load data into project failed, project reset for safety (%s)' % E.__repr__())
        self.refreshVisibility()

    def importProject(self):
        self.resetAttributes()      # This call is to a derived method, representing a full GUI reset.
        import_path = sanitizedFileDialog(self, 'open', 'Import project from file')
        try:
            self.project = loadObject(filename=import_path)
            self.tvec = np.arange(self.project.settings.tvec_start, self.project.settings.tvec_end + 1.0 / 2)
            PlottingSettings("gui") # updates rcParams for this instance
            self.acknowledgeProject()
            self.status = ('Status: Project "%s" successfully imported' % self.project.name)
        except:
            self.status = ('Status: Attempt to import project failed')
        self.refreshVisibility()

    def exportProject(self):
        export_path = sanitizedFileDialog(self, 'save', 'Export project to file')
        try:
            saveObject(filename=export_path, obj=self.project)
            self.status = ('Status: Project "%s" successfully exported' % self.project.name)
        except:
            self.status = ('Status: Attempt to export Project failed')
        self.refreshVisibility()

    def refreshStatus(self):
        if not self.guard_status:
            self.statusBar().showMessage(self.status)

# %% GUI class for result plotting

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

        self.plot_window = None         # A new window for displaying plots.


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

            # Clear all figures in the canvas.
            for i in reversed(range(self.plotter_layout.count())):
                self.plotter_layout.itemAt(i).widget().setParent(None)


    # The following wrapper function can be overloaded by derived classes.
    # Intended usage is to extend the user interface ahead of result plotting.
    def developLayout(self, layout): self.developLayoutResultPlotter()
    def developLayoutResultPlotter(self):
        pass


    def initUIResultPlotter(self):
        self.resetAttributesResultPlotter()    # This call is specific to base class attributes.

        # Widgets.
        self.label_plotter_result_1 = qtw.QLabel('Select first result: ')
        self.combo_plotter_result_1 = qtw.QComboBox(self)
        self.combo_plotter_result_1.activated[str].connect(self.selectResultOne)
        self.label_plotter_result_2 = qtw.QLabel('Select second result: ')
        self.combo_plotter_result_2 = qtw.QComboBox(self)
        self.combo_plotter_result_2.activated[str].connect(self.selectResultTwo)
        self.label_plotter_charac = qtw.QLabel('Select characteristic: ')
        self.combo_plotter_charac = qtw.QComboBox(self)
        self.combo_plotter_charac.activated[str].connect(self.selectCharacteristic)
        self.label_plotter_pop = qtw.QLabel('Select population: ')
        self.combo_plotter_pop = qtw.QComboBox(self)
        self.combo_plotter_pop.activated[str].connect(self.selectPopulation)
        self.button_plotter = qtw.QPushButton('Plot figures', self)
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
#        self.splitter_total = qtw.QSplitter()
#        self.splitter_total.addWidget(self.first_half)
#        self.splitter_total.addWidget(self.second_half)
        self.setCentralWidget(self.first_half)

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
            for result_name in self.project.results:
                self.combo_result_dict[result_name] = rid
                combo_box.addItem(result_name)
                rid += 1
            try: combo_box.setCurrentIndex(self.combo_result_dict[combo_name])
            except: pass    # Should be triggered if there are no results.

    def refreshCharacComboBox(self):
        self.combo_charac_dict = {}
        self.combo_plotter_charac.clear()
        cid = 0
        for charac_label in self.project.settings.charac_specs:
            charac_name = self.project.settings.charac_specs[charac_label]['name']
            self.combo_charac_dict[charac_name] = cid
            self.combo_plotter_charac.addItem(charac_name)
            cid += 1

        for flow_label in self.project.settings.linkpar_specs:
            
            if 'tag' in self.project.settings.linkpar_specs[flow_label]:
                flow_name = self.project.settings.linkpar_specs[flow_label]['name']
                self.combo_charac_dict[flow_name] = cid
                self.combo_plotter_charac.addItem(flow_name)
                cid += 1
            
            # Add flowrates that are publishable by 'output' tag
            elif 'output' in self.project.settings.linkpar_specs[flow_label] and self.project.settings.linkpar_specs[flow_label]['output'] == 'y':
                flow_name = self.project.settings.linkpar_specs[flow_label]['name']
                self.combo_charac_dict[flow_name] = cid
                self.combo_plotter_charac.addItem(flow_name)
                cid += 1


        self.combo_plotter_charac.setCurrentIndex(0)    # Should be triggered if there are no results.

        if self.charac_plot_name is None:
            self.charac_plot_name = str(self.combo_plotter_charac.itemText(self.combo_plotter_charac.currentIndex()))
        self.combo_plotter_charac.setCurrentIndex(self.combo_charac_dict[self.charac_plot_name])

    def refreshPopComboBox(self):
        self.combo_pop_dict = {}
        self.combo_plotter_pop.clear()
        pid = 0
        # <------------- Added total populations
        self.combo_pop_dict['Total populations'] = pid
        self.combo_plotter_pop.addItem('Total populations')
        pid += 1
        # <-------------
        for pop_name in self.project.data['pops']['name_labels']:
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

        print "------"
        print self.result_1_plot_name
        print self.result_2_plot_name
        print "------"

        if self.charac_plot_name is not None and self.pop_plot_name is not None:

            try:
                pop_plot_label = [self.project.data['pops']['name_labels'][self.pop_plot_name]]
                plot_observed = True
            except:
                logger.info("Could not identify a population. Setting to total populations")
                pop_plot_label = None
                plot_observed = False

            result_set = odict()
            result_set['%s' % self.result_1_plot_name] = self.project.results[self.result_1_plot_name]
            result_set['%s' % self.result_2_plot_name] = self.project.results[self.result_2_plot_name]

            try:
                plot_label = self.project.settings.charac_name_labels[self.charac_plot_name]
                logger.info("GUI plot characteristic: %s" % plot_label)
            except:
                try:
                    plot_label = self.project.settings.linkpar_name_labels[self.charac_plot_name]
                    try:
                        # Convert the parameter label to a transition tag if possible, so as to derive effective flow rate.
                        # TODO: Sort out getValuesAt() structure at some point to avoid potential label/tag confusion.
                        plot_label = self.project.settings.linkpar_specs[plot_label]['tag']
                        logger.info("GUI plot effective flow: %s" % plot_label)
                    except:
                        logger.info("GUI plot 'output' parameter: %s" % plot_label)
                except:
                    logger.info('Unable to plot for "%s"' % self.charac_plot_name)

            figure = plotCompareResults(self.project,
                                       result_set,
                                       output_labels=[plot_label],
                                       pop_labels=pop_plot_label,
                                       plot_observed_data=plot_observed,
                                       plot_total=True,
                                       colors=colors,
                                       save_fig=False)

            canvas = FigureCanvasQTAgg(figure)

            # Create new plotting window
            self.plot_window = qtw.QMainWindow(self)
            self.plot_window.setWindowTitle('%s - %s' % (self.charac_plot_name, self.pop_plot_name))
            widget = qtw.QDesktopWidget()
            screen = widget.availableGeometry()
            self.plot_window.resize(screen.width() * 4.0 / 10.0, screen.height() * 4.0 / 10.0)
            self.plot_window.setGeometry((screen.width() - self.plot_window.width()) / 2, (screen.height() - self.plot_window.height()) / 2,
                                         self.plot_window.width(), self.plot_window.height())
            self.plot_window.setCentralWidget(canvas)
            self.plot_window.show()

# %% GUI class for parameter calibration

class GUICalibration(GUIResultPlotterIntermediate):

    def __init__(self):
        super(GUICalibration, self).__init__()
        self.initUICalibration()


    def initUICalibration(self):
        self.setWindowTitle('Parameter calibration')
        self.resetAttributes()
        self.refreshVisibility()
        self.show()


    def resetAttributes(self):
        self.resetAttributesProjectManager()
        self.resetAttributesResultPlotter()

        self.parset = None      # This is the ParameterSet object that stores all edits made in the GUI.
                                # It should be an (edited) copy, not a reference, to an existing Project ParameterSet.
        self.parset_name = None
        
        self.combo_parset_dict = {}     # Dictionary that maps Parset names to indices used in combo boxes.
        self.calibration_id_dict = {}   # Dictionary that maps custom calibration-widget label to parameter and population labels.
        self.par_rows_dict = {}         # Dictionary that maps parameter and characteristic label to the rows of the parset table dedicated to them.
        self.fitted_characs_dict = {}   # Dictionary with characteristic-population label pairs as keys, denoting fitting metric inclusions for autocalibration.

        self.check_option = 'one'   # String that denotes how parset checkboxes will be ticked.
                                    # Can be 'one', 'par' or 'all', depending whether selection is standard, grouped-across-populations or all-or-nothing.


    def refreshVisibility(self):
        self.refreshVisibilityProjectManager()
        self.refreshVisibilityResultPlotter()

        is_project_loaded = self.project is not None
        does_parset_exist = is_project_loaded and len(self.project.parsets) > 0

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
        self.label_check_options.setVisible(does_parset_exist)
        self.radio_check_normal.setVisible(does_parset_exist)
        self.radio_check_par.setVisible(does_parset_exist)
        self.radio_check_all.setVisible(does_parset_exist)
        self.table_calibration.setVisible(does_parset_exist)

        self.label_autocalibrate.setVisible(does_parset_exist)
        self.edit_autocalibrate.setVisible(does_parset_exist)
        self.button_autocalibrate.setVisible(does_parset_exist)
        self.label_overwrite.setVisible(does_parset_exist)
        self.edit_overwrite.setVisible(does_parset_exist)
        self.button_overwrite.setVisible(does_parset_exist)
        self.label_model_run.setVisible(does_parset_exist)
        self.edit_model_run.setVisible(does_parset_exist)
        self.button_model_run.setVisible(does_parset_exist)


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
        self.label_parset = qtw.QLabel('Parameter set to edit: ')
        self.combo_parset = qtw.QComboBox(self)
        self.combo_parset.activated[str].connect(self.loadCalibration)
        
        self.label_check_options = qtw.QLabel('Toggle Checkbox Selection: ')
        self.radio_check_normal = qtw.QRadioButton('Normally')
        self.radio_check_normal.setChecked(True)
        self.radio_check_normal.toggled.connect(lambda:self.checkOptionState(self.radio_check_normal))
        self.radio_check_par = qtw.QRadioButton('Across populations')
        self.radio_check_par.toggled.connect(lambda:self.checkOptionState(self.radio_check_par))
        self.radio_check_all = qtw.QRadioButton('By all possible parameters')
        self.radio_check_all.toggled.connect(lambda:self.checkOptionState(self.radio_check_all))

        self.table_calibration = qtw.QTableWidget()
        self.table_calibration.cellChanged.connect(self.updateParset)

        self.label_autocalibrate = qtw.QLabel('Autocalibration time in seconds... ')
        self.edit_autocalibrate = qtw.QLineEdit()
        self.edit_autocalibrate.setText(str(10))
        self.button_autocalibrate = qtw.QPushButton('Autocalibrate parameters', self)
        self.button_autocalibrate.clicked.connect(self.autocalibrate)

        self.label_overwrite = qtw.QLabel('Save edits to... ')
        self.edit_overwrite = qtw.QLineEdit()
        self.button_overwrite = qtw.QPushButton('Save calibration', self)
        self.button_overwrite.clicked.connect(self.saveCalibration)

        self.label_model_run = qtw.QLabel('Run & save calibration results as... ')
        self.edit_model_run = qtw.QLineEdit()
        self.button_model_run = qtw.QPushButton('Generate results', self)
        self.button_model_run.clicked.connect(self.runCalibration)

        # Layout.
        grid_parset_load = qtw.QGridLayout()
        grid_parset_load.setSpacing(10)

        grid_parset_load.addWidget(self.label_parset, 0, 0)
        grid_parset_load.addWidget(self.combo_parset, 0, 1)
        
        grid_check_option = qtw.QGridLayout()
        grid_check_option.setSpacing(10)
        
        grid_check_option.addWidget(self.label_check_options, 0, 0)
        grid_check_option.addWidget(self.radio_check_normal, 0, 1)
        grid_check_option.addWidget(self.radio_check_par, 1, 1)
        grid_check_option.addWidget(self.radio_check_all, 2, 1)

        grid_parset_save = qtw.QGridLayout()
        grid_parset_save.setSpacing(10)

        grid_parset_save.addWidget(self.label_autocalibrate, 0, 0)
        grid_parset_save.addWidget(self.edit_autocalibrate, 0, 1)
        grid_parset_save.addWidget(self.button_autocalibrate, 0, 2)
        grid_parset_save.addWidget(self.label_overwrite, 1, 0)
        grid_parset_save.addWidget(self.edit_overwrite, 1, 1)
        grid_parset_save.addWidget(self.button_overwrite, 1, 2)
        grid_parset_save.addWidget(self.label_model_run, 2, 0)
        grid_parset_save.addWidget(self.edit_model_run, 2, 1)
        grid_parset_save.addWidget(self.button_model_run, 2, 2)

        layout.addLayout(grid_parset_load)
        layout.addLayout(grid_check_option)
        layout.addWidget(self.table_calibration)
        layout.addLayout(grid_parset_save)

    def checkOptionState(self, button):
        if button.text() == 'Normally':
            if button.isChecked() == True:
                self.status = ('Status: Parameters in the table will be ticked and unticked individually')
                self.check_option = 'one'
        elif button.text() == 'Across populations':
            if button.isChecked() == True:
                self.status = ('Status: Parameters in the table will be ticked and unticked as a group across populations')
                self.check_option = 'par'
        elif button.text() == 'By all possible parameters':
            if button.isChecked() == True:
                self.status = ('Status: Parameters in the table will be ticked and unticked across the entire parset')
                self.check_option = 'all'
        self.refreshStatus()    

    def refreshParsetComboBox(self):
        self.combo_parset_dict = {}
        self.combo_parset.clear()
        pid = 0
#        print self.project.parsets.keys()
        for parset_name in self.project.parsets:
            self.combo_parset_dict[parset_name] = pid
            self.combo_parset.addItem(parset_name)
            pid += 1
        try: self.combo_parset.setCurrentIndex(self.combo_parset_dict[self.parset_name])
        except: pass

    def loadCalibration(self, parset_name, delay_refresh=False):
        self.parset_name = str(parset_name)
        self.parset = dcp(self.project.parsets[self.parset_name])
        self.fitted_characs_dict = {pair:True for pair in it.product(self.parset.par_ids['characs'].keys(), self.parset.pop_labels)}
        self.status = ('Status: Parameter set "%s" selected for editing' % self.parset_name)
        if not delay_refresh:
            self.refreshVisibility()

    def saveCalibration(self):
        parset_name = str(self.edit_overwrite.text())
        if parset_name == '':
            self.status = ('Status: Attempt to save parameter set failed, no name provided')
        else:
            if parset_name in self.project.parsets:
                self.status = ('Status: Parameter set "%s" successfully overwritten' % parset_name)
            else:
                self.status = ('Status: New parameter set "%s" added to project' % parset_name)
            self.parset.name = parset_name
            self.project.parsets[parset_name] = dcp(self.parset)
        self.refreshVisibility()

    def runCalibration(self):
        self.status = ('Status: Running model for parameter set "%s"' % self.parset_name)
        self.refreshStatus()
        result_name = str(self.edit_model_run.text())
        if result_name == '':
            result_name = None
        self.project.runSim(parset=self.parset, store_results=True, result_type='calibration', result_name=result_name)
        self.acknowledgeResults()
        self.status = ('Status: Model successfully processed for parameter set "%s"' % self.parset_name)
        self.refreshVisibility()

    def autocalibrate(self):
        try:
            calibration_time = float(str(self.edit_autocalibrate.text()))
            if calibration_time < 0:
                raise Exception('autocalibration time cannot be negative')
        except Exception as E:
            self.status = ('Status: Autocalibration aborted because "%s"' % E.message)
            self.refreshStatus()
            return
        self.status = ('Status: Autocalibrating checked selection of parameter set "%s" for %s seconds' % (self.parset.name, str(calibration_time)))
        self.refreshStatus()
        try:
            self.parset = self.project.runAutofitCalibration(parset=self.parset, new_parset_name=self.parset.name, target_characs=self.fitted_characs_dict.keys(), max_time=calibration_time, save_parset=False)
            self.status = ('Status: Autocalibration complete (but unsaved) for parameter set "%s"' % self.parset.name)
        except Exception as E:
            self.status = ('Status: Autocalibration was unsuccessful, perhaps because no parameters were selected to calibrate or no parameter-associated data was chosen to fit against')
        self.refreshVisibility()

    def makeParsetTable(self):
        self.table_calibration.setVisible(False)    # Resizing columns requires table to be hidden first.
        self.table_calibration.clear()

        # Disconnect the calibration table from cell change signals to avoid signal flooding during connection.
        try: self.table_calibration.cellChanged.disconnect()
        except: pass

        parset = self.parset
        num_pops = len(parset.pop_labels)
        row_count = num_pops * (len(parset.pars['cascade']) + len(parset.pars['characs']) - len(self.project.settings.par_funcs))
        self.table_calibration.setRowCount(row_count)
        self.table_calibration.setColumnCount(len(self.tvec))

        k = 0
        custom_ids = []
        for par_type in ['characs', 'cascade']:
            for par in parset.pars[par_type]:
                if ((par_type == 'cascade' and par.label not in self.project.settings.par_funcs) or 
                    (par_type == 'characs')): #and 'entry_point' in self.project.settings.charac_specs[par.label])):
                    for pid in xrange(len(parset.pop_labels)):
                        pop_label = parset.pop_labels[pid]
                        try:
                            par_name = self.project.settings.linkpar_specs[par.label]['name']
                        except:
                            par_name = self.project.settings.charac_specs[par.label]['name']
                        custom_id = par_name + '\n' + parset.pop_names[pid]
                        custom_ids.append(custom_id)
                        self.calibration_id_dict[custom_id] = {'par_label':par.label, 'pop_label':pop_label}
                        if par.label not in self.par_rows_dict:
                            self.par_rows_dict[par.label] = {}
                        self.par_rows_dict[par.label][k * num_pops + pid] = True
                        
                        # Create autocalibration checkbox column.
                        temp = qtw.QTableWidgetItem()
                        if par_type == 'characs' and 'entry_point' not in self.project.settings.charac_specs[par.label]:
                            temp.setFlags(qtc.Qt.NoItemFlags)
                        else:
                            temp.setText(str(par.y_factor[pid]))
                            temp.setTextAlignment(qtc.Qt.AlignCenter)
                            temp.setFlags(qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsUserCheckable)
                            if par.autocalibrate[pid]:
                                temp.setCheckState(qtc.Qt.Checked)
                            else:
                                temp.setCheckState(qtc.Qt.Unchecked)
                        self.table_calibration.setItem(k * num_pops + pid, 0, temp)
                        
                        # Create data fitting checkbox column.
                        temp = qtw.QTableWidgetItem()
                        if par_type == 'characs':
                            temp.setText(str(1))    # Until calibration weighting is implemented, the default uneditable value will be 1.
                            temp.setTextAlignment(qtc.Qt.AlignCenter)
                            temp.setFlags(qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsUserCheckable)
                            if (par.label,pop_label) in self.fitted_characs_dict:
                                temp.setCheckState(qtc.Qt.Checked)
                            else:
                                temp.setCheckState(qtc.Qt.Unchecked)
                        else:
                            temp.setFlags(qtc.Qt.NoItemFlags)
                        self.table_calibration.setItem(k * num_pops + pid, 1, temp)

                        # Insert the actual values.
                        for eid in xrange(len(par.t[pid])):
                            t = par.t[pop_label][eid]
                            y = par.y[pop_label][eid]
                            temp = qtw.QTableWidgetItem()
                            temp.setText(str(y))
                            temp.setTextAlignment(qtc.Qt.AlignCenter)
                            self.table_calibration.setItem(k * num_pops + pid, 2 + int(t) - self.tvec[0], temp)
                    k += 1
        self.table_calibration.setVerticalHeaderLabels(custom_ids)
        self.table_calibration.setHorizontalHeaderLabels(['Parameter Scaling Factor\n(Tick: Autocalibrate Parameter)','Data Weighting Factor\n(Tick: Use in Fitting Metric)'] + [str(int(x)) for x in self.tvec])
        self.table_calibration.resizeColumnsToContents()
        self.table_calibration.resizeRowsToContents()

        self.table_calibration.cellChanged.connect(self.updateParset)


    def updateParset(self, row, col):
        custom_id = str(self.table_calibration.verticalHeaderItem(row).text())
        par_label = self.calibration_id_dict[custom_id]['par_label']
        pop_label = self.calibration_id_dict[custom_id]['pop_label']
        par = self.parset.getPar(par_label)

        new_val_str = str(self.table_calibration.item(row, col).text())
        if col == 0:
            calib_prefix = ''
            if self.table_calibration.item(row, col).checkState() == qtc.Qt.Checked:
                if not self.guard_status:   # Used here, this is a guard against recursion.
                    self.guard_status = True
                    if self.check_option == 'all':
                        for other_row in xrange(self.table_calibration.rowCount()):
                            if not self.table_calibration.item(other_row, col).flags() == qtc.Qt.NoItemFlags:
                                self.table_calibration.item(other_row, col).setCheckState(qtc.Qt.Checked)
                    elif self.check_option == 'par':
                        for other_row in self.par_rows_dict[par_label]:
                            self.table_calibration.item(other_row, col).setCheckState(qtc.Qt.Checked)
                    self.guard_status = False
                par.autocalibrate[pop_label] = True
            else:
                if not self.guard_status:   # Used here, this is a guard against recursion.
                    self.guard_status = True
                    if self.check_option == 'all':
                        for other_row in xrange(self.table_calibration.rowCount()):
                            if not self.table_calibration.item(other_row, col).flags() == qtc.Qt.NoItemFlags:
                                self.table_calibration.item(other_row, col).setCheckState(qtc.Qt.Unchecked)
                    elif self.check_option == 'par':
                        for other_row in self.par_rows_dict[par_label]:
                            self.table_calibration.item(other_row, col).setCheckState(qtc.Qt.Unchecked)
                    self.guard_status = False
                par.autocalibrate[pop_label] = False
                calib_prefix = 'un'

            try:
                new_val = float(new_val_str)
                if new_val < 0:
                    raise Exception('scaling factor is negative')
            except Exception as E:
                self.status = ('Status: Attempt to edit scaling factor failed because "%s"' % E.message)
                new_val = DEFAULT_YFACTOR
                self.refreshStatus()
                self.guard_status = True
                self.table_calibration.item(row, col).setText(str(new_val))
                self.guard_status = False
                return
            par.y_factor[pop_label] = new_val
            self.status = ('Status: Current edited parameter set has %smarked parameter "%s", population "%s", for autocalibration, with scaling factor "%f"' % (calib_prefix, par_label, pop_label, new_val))
        elif col == 1:
            calib_prefix = 'in'
            if self.table_calibration.item(row, col).checkState() == qtc.Qt.Checked:
                if not self.guard_status:   # Used here, this is a guard against recursion.
                    self.guard_status = True
                    if self.check_option == 'all':
                        for other_row in xrange(self.table_calibration.rowCount()):
                            if not self.table_calibration.item(other_row, col).flags() == qtc.Qt.NoItemFlags:
                                self.table_calibration.item(other_row, col).setCheckState(qtc.Qt.Checked)
                    elif self.check_option == 'par':
                        for other_row in self.par_rows_dict[par_label]:
                            self.table_calibration.item(other_row, col).setCheckState(qtc.Qt.Checked)
                    self.guard_status = False
                self.fitted_characs_dict[(par_label,pop_label)] = True
            else:
                if not self.guard_status:   # Used here, this is a guard against recursion.
                    self.guard_status = True
                    if self.check_option == 'all':
                        for other_row in xrange(self.table_calibration.rowCount()):
                            if not self.table_calibration.item(other_row, col).flags() == qtc.Qt.NoItemFlags:
                                self.table_calibration.item(other_row, col).setCheckState(qtc.Qt.Unchecked)
                    elif self.check_option == 'par':
                        for other_row in self.par_rows_dict[par_label]:
                            self.table_calibration.item(other_row, col).setCheckState(qtc.Qt.Unchecked)
                    self.guard_status = False
                try: del self.fitted_characs_dict[(par_label,pop_label)]
                except: pass
                calib_prefix = 'ex'

            self.status = ('Status: Current edited parameter set has %scluded parameter "%s", population "%s", in the autocalibration data-fitting metric' % (calib_prefix, par_label, pop_label))
        else:
            year = float(str(self.table_calibration.horizontalHeaderItem(col).text()))
            try:
                new_val = float(new_val_str)
            except:
                if new_val_str == '':
                    remove_success = par.removeValueAt(t=year, pop_label=pop_label)
                    if remove_success:
                        self.status = ('Status: Value successfully deleted from parameter set')
                        self.refreshStatus()
                        return
                    else: self.status = ('Status: Attempt to remove item in parameter set failed, at least one value per row required')
                else: self.status = ('Status: Attempt to edit item in parameter set failed, only numbers allowed')
                self.refreshStatus()
                self.guard_status = True
                self.table_calibration.item(row, col).setText(str(par.interpolate(tvec=[year], pop_label=pop_label)[0]))
                self.guard_status = False
                return
            par.insertValuePair(t=year, y=new_val, pop_label=pop_label)
            self.status = ('Status: Current edited parameter set uses value "%f" for parameter "%s", population "%s", year "%i"' % (new_val, par_label, pop_label, year))
        self.refreshStatus()


# %% GUI class for program reconciliation

class GUIReconciliation(GUIResultPlotterIntermediate):

    def __init__(self):
        super(GUIReconciliation, self).__init__()
        self.initUIReconciliation()


    def initUIReconciliation(self):
        self.setWindowTitle('Program reconciliation')
        self.resetAttributes()
        self.refreshVisibility()
        self.show()


    def resetAttributes(self):
        self.resetAttributesProjectManager()
        self.resetAttributesResultPlotter()
        
        self.progset = None             # This is the ProgramSet object that stores all edits made in the GUI.
                                        # It should be an (edited) copy, not a reference, to an existing Project ProgramSet.

        self.parset_name = None
        self.progset_name = None
        self.combo_parset_dict = {}     # Dictionary that maps ParameterSet names to indices used in combo boxes.
        self.combo_progset_dict = {}    # Dictionary that maps ProgramSet names to indices used in combo boxes.

        self.options = None                 # The options dictionary for running a budget scenario, i.e. a standard simulation with program-based parameter overwrites.
        self.reconciliation_id_dict = {}    # Dictionary that maps custom reconciliation-widget label to program label.
        self.sigma_dict = {}                # Dictionary that maps program labels with 'unit_cost', 'budget' and 'attribute' sigmas for program reconciliation.
        self.attribute_labels = []          # A list of attribute labels relevant to programs with a loaded ProgramSet.
        self.attribute_names = []           # A corresponding list of attribute names, assuming that each label/name pair is unique.
                                            # TODO: Reconsider this assumption.
                          
        self.check_option = 'one'   # String that denotes how progset checkboxes will be ticked.
                                    # Can be 'one' or 'all', depending whether selection is standard, grouped-across-populations or all-or-nothing.
                                   

    def refreshVisibility(self):
        self.refreshVisibilityProjectManager()
        self.refreshVisibilityResultPlotter()

        is_project_loaded = self.project is not None
        does_parset_exist = is_project_loaded and len(self.project.parsets) > 0
        does_progset_exist = is_project_loaded and len(self.project.progsets) > 0
#        do_options_exist = self.options is not None

        if is_project_loaded:
            self.refreshParsetComboBox()
            self.refreshProgsetComboBox()
        self.label_parset.setVisible(is_project_loaded)
        self.combo_parset.setVisible(is_project_loaded)
        self.label_progset.setVisible(is_project_loaded)
        self.combo_progset.setVisible(is_project_loaded)

        policy_min = qtw.QSizePolicy.Minimum
        policy_exp = qtw.QSizePolicy.Expanding
        if does_progset_exist:
            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_min)
            self.makeProgsetTable()
        else:
            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_exp)
        self.checkbox_align_sigmas.setVisible(does_progset_exist)
        self.table_reconciliation.setVisible(does_progset_exist)
        
        self.label_year_start.setVisible(does_progset_exist)
        self.edit_year_start.setVisible(does_progset_exist)
        self.label_reconcile.setVisible(does_parset_exist & does_progset_exist)
        self.edit_reconcile.setVisible(does_parset_exist & does_progset_exist)
        self.button_reconcile.setVisible(does_parset_exist & does_progset_exist)
        self.label_overwrite.setVisible(does_progset_exist)
        self.edit_overwrite.setVisible(does_progset_exist)
        self.button_overwrite.setVisible(does_progset_exist)
        self.label_model_run.setVisible(does_parset_exist & does_progset_exist)
        self.edit_model_run.setVisible(does_parset_exist & does_progset_exist)
        self.button_model_run.setVisible(does_parset_exist & does_progset_exist)


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


    # While UI initialisation can extend the interface, this method is where widgets for the core process should be set up.
    def developLayout(self, layout):
        self.developLayoutResultPlotter()

        # Widgets.
        self.label_parset = qtw.QLabel('Parameter set to use: ')
        self.combo_parset = qtw.QComboBox(self)
        self.combo_parset.activated[str].connect(self.loadCalibration)

        self.label_progset = qtw.QLabel('Program set to use: ')
        self.combo_progset = qtw.QComboBox(self)
        self.combo_progset.activated[str].connect(self.loadPrograms)

        self.label_year_start = qtw.QLabel('Year to reconcile calibration with fixed-budget programs... ')
        self.edit_year_start = qtw.QLineEdit()
        self.edit_year_start.returnPressed.connect(self.updateStartYear)
                
        self.checkbox_align_sigmas = qtw.QCheckBox('Duplicate sigma checkbox choices and values across programs')
        self.checkbox_align_sigmas.stateChanged.connect(self.checkOptionState)

        grid_check_option = qtw.QGridLayout()
        grid_check_option.setSpacing(10)
        grid_check_option.addWidget(self.checkbox_align_sigmas, 0, 0)
        
        self.label_reconcile = qtw.QLabel('Automatic reconciliation time in seconds... ')
        self.edit_reconcile = qtw.QLineEdit()
        self.edit_reconcile.setText(str(10))
        self.button_reconcile = qtw.QPushButton('Reconcile parameters and programs', self)
        self.button_reconcile.clicked.connect(self.reconcile)
        
        self.label_overwrite = qtw.QLabel('Save edits to... ')
        self.edit_overwrite = qtw.QLineEdit()
        self.button_overwrite = qtw.QPushButton('Save program set', self)
        self.button_overwrite.clicked.connect(self.savePrograms)

        self.label_model_run = qtw.QLabel('Run & save program-based simulation results as... ')
        self.edit_model_run = qtw.QLineEdit()
        self.button_model_run = qtw.QPushButton('Generate results', self)
        self.button_model_run.clicked.connect(self.runProgsetSim)

        # Layout.
        grid_progset_load = qtw.QGridLayout()
        grid_progset_load.setSpacing(10)

        grid_progset_load.addWidget(self.label_parset, 0, 0)
        grid_progset_load.addWidget(self.combo_parset, 0, 1)
        grid_progset_load.addWidget(self.label_progset, 1, 0)
        grid_progset_load.addWidget(self.combo_progset, 1, 1)

        grid_progset_load.addWidget(self.label_year_start, 2, 0)
        grid_progset_load.addWidget(self.edit_year_start, 2, 1)

        grid_progset_save = qtw.QGridLayout()
        grid_progset_save.setSpacing(10)

        grid_progset_save.addWidget(self.label_reconcile, 0, 0)
        grid_progset_save.addWidget(self.edit_reconcile, 0, 1)
        grid_progset_save.addWidget(self.button_reconcile, 0, 2)
        grid_progset_save.addWidget(self.label_overwrite, 1, 0)
        grid_progset_save.addWidget(self.edit_overwrite, 1, 1)
        grid_progset_save.addWidget(self.button_overwrite, 1, 2)
        grid_progset_save.addWidget(self.label_model_run, 2, 0)
        grid_progset_save.addWidget(self.edit_model_run, 2, 1)
        grid_progset_save.addWidget(self.button_model_run, 2, 2)
        
        self.table_reconciliation = qtw.QTableWidget()

        layout.addLayout(grid_progset_load)
        layout.addLayout(grid_check_option)
        layout.addWidget(self.table_reconciliation)
        layout.addLayout(grid_progset_save)
        
    def checkOptionState(self, state):
        if not state == qtc.Qt.Checked:
            self.status = ('Status: Unit-cost, budget and attribute sigmas for programs will be ticked, unticked and edited individually')
            self.check_option = 'one'
        elif state == qtc.Qt.Checked:
            self.status = ('Status: Unit-cost, budget and attribute sigmas, for the individual type being edited, will be duplicated across all programs in terms of selection and value')
            self.check_option = 'all'
        self.refreshStatus()

    def updateStartYear(self):
        try:    
            self.options['progs_start'] = float(str(self.edit_year_start.text()))
        except:
            self.status = ('Status: An invalid program start-year was chosen, so the value has reverted to default')
            self.options['progs_start'] = defaultOptimOptions(settings=self.project.settings, progset=self.progset)['progs_start']
            self.edit_year_start.setText(str(self.options['progs_start']))
        self.refreshVisibility()

    def refreshParsetComboBox(self):
        self.combo_parset_dict = {}
        self.combo_parset.clear()
        pid = 0
        for parset_name in self.project.parsets:
            self.combo_parset_dict[parset_name] = pid
            self.combo_parset.addItem(parset_name)
            pid += 1
        try: self.combo_parset.setCurrentIndex(self.combo_parset_dict[self.parset_name])
        except: pass

    def refreshProgsetComboBox(self):
        self.combo_progset_dict = {}
        self.combo_progset.clear()
        pid = 0
        for progset_name in self.project.progsets:
            self.combo_progset_dict[progset_name] = pid
            self.combo_progset.addItem(progset_name)
            pid += 1
        try: self.combo_progset.setCurrentIndex(self.combo_progset_dict[self.progset_name])
        except: pass

    def loadCalibration(self, parset_name, delay_refresh=False):
        self.parset_name = str(parset_name)
        self.status = ('Status: Parameter set "%s" selected for program reconciliation' % self.parset_name)
        if not delay_refresh:
            self.refreshVisibility()

    def loadPrograms(self, progset_name, delay_refresh=False):
        self.progset_name = str(progset_name)
        self.progset = dcp(self.project.progsets[self.progset_name])
        self.options = defaultOptimOptions(settings=self.project.settings, progset=self.progset)
        
        # Generate a zero-value sigma dictionary and progset attribute list for any progset reload.
        self.sigma_dict = {}
        attribute_label_dict = {}
        for prog in self.progset.progs:
            if not prog.func_specs['type'] == 'cost_only':
                self.sigma_dict[prog.label] = {'unit_cost':0.0, 'budget':0.0, 'attribute':0.0}
            for attribute_label in prog.attributes:
                if not attribute_label in attribute_label_dict:
                    attribute_label_dict[attribute_label] = self.project.settings.progtype_specs[prog.prog_type]['attribute_label_names'][attribute_label]
        self.attribute_labels = attribute_label_dict.keys()
        self.attribute_names = [attribute_label_dict[x] for x in self.attribute_labels]
                    
        self.edit_year_start.setText(str(self.options['progs_start']))
        self.status = ('Status: Program set "%s" selected for program reconciliation' % self.progset_name)
        if not delay_refresh:
            self.refreshVisibility()
            
    def savePrograms(self):
        progset_name = str(self.edit_overwrite.text())
        if progset_name == '':
            self.status = ('Status: Attempt to save program set failed, no name provided')
        else:
            if progset_name in self.project.progsets:
                self.status = ('Status: Program set "%s" successfully overwritten' % progset_name)
            else:
                self.status = ('Status: New program set "%s" added to project' % progset_name)
            self.progset.name = progset_name
            self.project.progsets[progset_name] = dcp(self.progset)
        self.refreshVisibility()
        
    def reconcile(self):
        try:
            reconciliation_time = float(str(self.edit_reconcile.text()))
            if reconciliation_time < 0:
                raise Exception('reconciliation time cannot be negative')
        except Exception as E:
            self.status = ('Status: Reconciliation process aborted because "%s"' % E.message)
            self.refreshStatus()
            return
        self.status = ('Status: Reconciling checked selection of program set "%s" with parameter set "%s" for %s seconds' % (self.progset.name, self.parset_name, str(reconciliation_time)))
        self.refreshStatus()
        try:
            self.progset, output = self.project.reconcile(parset_name=self.parset_name, progset=self.progset, reconcile_for_year = self.options['progs_start'], sigma_dict=self.sigma_dict, overwrite=True, max_time=reconciliation_time, save_progset=False)
            self.status = ('Status: Reconciliation process complete (but unsaved) for program set "%s"' % self.progset.name)
            
            # Print reconciliation output to a new window.
            # TODO: Consider redesign so that all extra-window widgets use common specifications.
            self.output_window = qtw.QMainWindow(self)
            self.output_window.setWindowTitle('Reconciliation Output')
            widget = qtw.QDesktopWidget()
            screen = widget.availableGeometry()
            self.output_window.resize(screen.width() * 4.0 / 10.0, screen.height() * 4.0 / 10.0)
            self.output_window.setGeometry((screen.width() - self.output_window.width()) / 2, (screen.height() - self.output_window.height()) / 2,
                                         self.output_window.width(), self.output_window.height())
            
            log_output = qtw.QTextEdit(self.output_window)
            log_output.setReadOnly(True)
            log_output.setLineWrapMode(qtw.QTextEdit.NoWrap)
            
            log_output.insertPlainText(output)
            
            self.output_window.setCentralWidget(log_output)
            self.output_window.show()
        except Exception as E:
            self.status = ('Status: Reconciliation process was unsuccessful because "%s"' % E.message)#, perhaps because no parameters were selected to calibrate or no parameter-associated data was chosen to fit against')
        self.refreshVisibility()
            
    # TODO: Remove magic numbers once layout design is more stable.
    def makeProgsetTable(self):
        self.table_reconciliation.setVisible(False)    # Resizing columns requires table to be hidden first.
        self.table_reconciliation.clear()

        # Disconnect the reconciliation table from cell change signals to avoid signal flooding during connection.
        try: self.table_reconciliation.cellChanged.disconnect()
        except: pass

        progset = self.progset
        row_count = len(progset.progs)
        self.table_reconciliation.setRowCount(row_count)
        self.table_reconciliation.setColumnCount(7+len(self.attribute_labels))        # Unit cost and sigma, total budget and sigma, effective coverage, historical coverage, impact sigma.

        row_id = 0
        custom_ids = []
        for prog in progset.progs:
            custom_ids.append(prog.name)
            if prog.name not in self.reconciliation_id_dict:
                self.reconciliation_id_dict[prog.name] = {}
            self.reconciliation_id_dict[prog.name]['prog_label'] = prog.label
#            self.reconciliation_id_dict[prog.name]['row'] = row_id
            for col_id in xrange(7):
                temp = qtw.QTableWidgetItem()
                if prog.func_specs['type'] == 'cost_only' and col_id != 2:
                    temp.setFlags(qtc.Qt.NoItemFlags)
                else:
                    temp.setTextAlignment(qtc.Qt.AlignCenter)
                    if col_id in [1,3,6]:
                        temp.setText(str(0))
                        temp.setFlags(qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable)# | qtc.Qt.ItemIsUserCheckable)
#                        temp.setCheckState(qtc.Qt.Unchecked)
                    elif col_id == 0:
                        temp.setText(str(prog.func_specs['pars']['unit_cost']))
                        temp.setFlags(qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable)
                    elif col_id == 2:
                        temp.setText(str(prog.getDefaultBudget(year=self.options['progs_start'])))
                        temp.setFlags(qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable)
                    elif col_id == 4:
                        temp.setText(str(prog.getCoverage(prog.getDefaultBudget(year=self.options['progs_start']))))
                        temp.setFlags(qtc.Qt.ItemIsEnabled)
                    elif col_id == 5:
                        temp.setText(str(prog.interpolate(tvec=[self.options['progs_start']], attributes=['cov'])['cov'][-1]))
                        temp.setFlags(qtc.Qt.ItemIsEnabled)
                self.table_reconciliation.setItem(row_id, col_id, temp)
            for attribute_id in xrange(len(self.attribute_labels)):
                attribute_label = self.attribute_labels[attribute_id]
                col_id = attribute_id + 7
                temp = qtw.QTableWidgetItem()
                temp.setTextAlignment(qtc.Qt.AlignCenter)
                temp.setFlags(qtc.Qt.ItemIsEnabled)
                if attribute_label in prog.attributes:
                    temp.setFlags(qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable)
                    temp.setText(str(prog.interpolate(tvec=[self.options['progs_start']], attributes=[attribute_label])[attribute_label][-1]))
                self.table_reconciliation.setItem(row_id, col_id, temp)
            row_id += 1
        self.table_reconciliation.setVerticalHeaderLabels(custom_ids)
        self.table_reconciliation.setHorizontalHeaderLabels(['Unit Cost','Unit Cost\nSigma','Program Budget\n(Year: %s)' % self.options['progs_start'],'Program Budget\nSigma','Effective Coverage\n(Year: %s)' % self.options['progs_start'],'Databook Coverage\n(Year: %s)' % self.options['progs_start'],'Attribute Sigma']+self.attribute_labels)
        self.table_reconciliation.resizeColumnsToContents()
        self.table_reconciliation.resizeRowsToContents()

        self.table_reconciliation.cellChanged.connect(self.updateProgset)
    
    # TODO: Uses the magic column numbers as above. If those values change, so should the ones below.
    def updateProgset(self, row, col):
        custom_id = str(self.table_reconciliation.verticalHeaderItem(row).text())
        prog_label = self.reconciliation_id_dict[custom_id]['prog_label']
        prog = self.progset.getProg(prog_label)
        
        col_to_sigma = {1:'unit_cost', 3:'budget', 6:'attribute'}

        new_val_str = str(self.table_reconciliation.item(row, col).text())
        if col in [1,3,6]:      # Hardcoded column values corresponding to checkable sigmas.
                  
            # Extract value from the table first, whether entered in or whether the accompanying checkbox was triggered.
            try:
                new_val = float(new_val_str)
                if new_val < 0:
                    raise Exception('sigma is negative')
            except Exception as E:
                self.status = ('Status: Attempt to edit sigma failed because "%s"' % E.message)
                new_val = 0.0
                self.refreshStatus()
                self.guard_status = True
                self.table_reconciliation.item(row, col).setText(str(new_val))
                self.guard_status = False
                return
            
            # Duplicate the sigma value if the relevant option is chosen.
            if not self.guard_status:   # Used here, this is a guard against recursion.
                self.guard_status = True
                if self.check_option == 'all':
                    for other_row in xrange(self.table_reconciliation.rowCount()):
                        if not self.table_reconciliation.item(other_row, col).flags() == qtc.Qt.NoItemFlags:
                            self.table_reconciliation.item(other_row, col).setText(str(new_val))
                self.guard_status = False
                
            self.sigma_dict[prog_label][col_to_sigma[col]] = new_val
            self.status = ('Status: Program reconciliation will allow program "%s" to vary in "%s" value by up to %f%%' % (prog_label, col_to_sigma[col], new_val*100))

        else:
            try:
                try: new_val = float(new_val_str)
                except: raise Exception('as only number inputs are allowed')
                if new_val < 0:
                    raise Exception('as values must be positive')
                if (col == 0 and new_val == 0):
                    raise Exception('as unit cost cannot be zero')
            except Exception as E:
                self.status = ('Status: Attempt to edit program set failed, %s' % E.message)
                self.refreshStatus()
                self.guard_status = True
                if col == 0:
                    self.table_reconciliation.item(row, col).setText(str(prog.func_specs['pars']['unit_cost']))
                elif col == 2:
                    self.table_reconciliation.item(row, col).setText(str(prog.getDefaultBudget(year=self.options['progs_start'])))
                else:
                    attribute_label = self.attribute_labels[col-7]
                    self.table_reconciliation.item(row, col).setText(str(prog.interpolate(tvec=[self.options['progs_start']], attributes=[attribute_label])[attribute_label][-1]))
                self.guard_status = False
                return
        
            if col == 0:
                prog.func_specs['pars']['unit_cost'] = new_val
                self.table_reconciliation.item(row, 4).setText(str(prog.getCoverage(prog.getDefaultBudget(year=self.options['progs_start']))))
                self.status = ('Status: Current edited program set uses unit cost "%f" for program "%s"' % (new_val, prog_label))
            elif col == 2:
                prog.insertValuePair(self.options['progs_start'], new_val, 'cost', rescale_after_year=True)
#                self.options['init_alloc'][prog_label] = new_val    # Update the options dictionary immediately rather than checking the table later.
                self.table_reconciliation.item(row, 4).setText(str(prog.getCoverage(prog.getDefaultBudget(year=self.options['progs_start']))))
                self.status = ('Status: Current edited program set associates program "%s" with a budget of "%f" in "%f"' % (prog_label, new_val, self.options['progs_start']))
            elif col > 6:
                attribute_label = self.attribute_labels[col-7]
                prog.insertValuePair(self.options['progs_start'], new_val, attribute_label, rescale_after_year=True)
                self.status = ('Status: Current edited program set associates program "%s" with an attribute "%s" of "%f" in "%f"' % (prog_label, attribute_label, new_val, self.options['progs_start']))
                
        self.refreshStatus()

    def runProgsetSim(self):
        if not self.updateOptions():
            self.status = ('Status: User-specified options could not be read into program-based simulation, so options are being reverted')
            self.edit_year_start.setText(str(self.options['progs_start']))
        else:
            self.status = ('Status: Running model for Parset "%s" and currently-displayed Progset "%s"' % (self.parset_name, self.progset.name))
            self.refreshStatus()
            result_name = str(self.edit_model_run.text())
            if result_name == '':
                result_name = None
            self.project.runSim(parset_name=self.parset_name, progset=self.progset, options=self.options, store_results=True, result_type='progsim', result_name=result_name)
            self.acknowledgeResults()
            self.status = ('Status: Model successfully processed for parameter set "%s" and program set "%s"' % (self.parset_name, self.progset_name))
        self.refreshVisibility()

    # Regenerates options dictionary if possible, which should update initial allocation based on ProgramSet edits.
    # Return True if successful or False if there was an error.
    def updateOptions(self):
        if self.options is None:
            return False
        else:
            try:
                progs_start = float(str(self.edit_year_start.text()))
                self.options = defaultOptimOptions(settings=self.project.settings, progset=self.progset, progs_start=progs_start)
            except: return False
            return True

# %% GUI class for running parameter scenarios

class GUIParameterScenario(GUIResultPlotterIntermediate):

    def __init__(self):
        super(GUIParameterScenario, self).__init__()
        self.initUIParameterScenario()


    def initUIParameterScenario(self):
        self.setWindowTitle('Parameter scenarios')
        self.resetAttributes()
        self.refreshVisibility()
        self.show()


    def resetAttributes(self):
        self.resetAttributesProjectManager()
        self.resetAttributesResultPlotter()

        self.parset_name = None
        self.scen_name = None
        self.combo_parset_dict = {}     # Dictionary that maps ParameterSet names to indices used in combo boxes.
        self.combo_scen_dict = {}    # Dictionary that maps Scen names to indices used in combo boxes.

        self.widget_pars_list = [ # List that maps cascade parameter names
                ('param_diag', 'Proportion of people diagnosed:'),
                ('param_link', 'Proportion of people linked to care:'),
                ('param_succ', 'Proportion of people successfully treated:'),
                ('param_fail', 'Proportion of people with treatment failure:'),
                ]

        self.cascade_pars = { # Dictionary for values
                'param_diag': 0.76, # TODO WARNING hard-coded
                'param_link': 0.76,
                'param_succ': 0.85,
                'param_fail': 0.05,
                }

        self.start_year = 2015.0 # WARNING, shouldn't hard-code


    def refreshVisibility(self):
        self.refreshVisibilityProjectManager()
        self.refreshVisibilityResultPlotter()

        is_project_loaded = self.project is not None
        does_parset_exist = is_project_loaded and len(self.project.parsets) > 0

        if is_project_loaded:
            self.refreshParsetComboBox()
            self.refreshScenComboBox()
        self.label_parset.setVisible(is_project_loaded)
        self.combo_parset.setVisible(is_project_loaded)
        self.label_scen.setVisible(is_project_loaded)
        self.combo_scen.setVisible(is_project_loaded)

#        policy_min = qtw.QSizePolicy.Minimum
#        policy_exp = qtw.QSizePolicy.Expanding
#        if is_project_loaded:
#            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_min)
#        else:
#            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_exp)
        self.scroll_area.setVisible(is_project_loaded)

        self.label_year_start.setVisible(is_project_loaded)
        self.edit_year_start.setVisible(is_project_loaded)
        self.label_model_run.setVisible(does_parset_exist)
        self.edit_model_run.setVisible(does_parset_exist)
        self.button_model_run.setVisible(does_parset_exist)

        # Update the visibility of widgets depending on if they have anything to show.


    def acknowledgeProject(self):
        self.acknowledgeProjectProjectManager()
        self.acknowledgeProjectResultPlotter()
        if not self.project is None:
            if len(self.project.parsets) < 1:
                self.project.makeParset(name='default')
            self.loadCalibration(self.project.parsets[0].name, delay_refresh=True)

        # If a project is loaded, do whatever initial pre-processing is needed.


    # While UI initialisation can extend the interface, this method is where widgets for the core process should be set up.
    def developLayout(self, layout):
        self.developLayoutResultPlotter()

        # Widgets.
        self.label_parset = qtw.QLabel('Parameter set to use: ')
        self.combo_parset = qtw.QComboBox(self)
        self.combo_parset.activated[str].connect(self.loadCalibration)

        self.label_scen = qtw.QLabel('Parameter scenario type: ')
        self.combo_scen = qtw.QComboBox(self)
#        self.combo_scen.activated[str].connect(self.loadPrograms)

        self.label_year_start = qtw.QLabel('Start year for parameter changes:')
        self.edit_year_start = qtw.QLineEdit()
        self.edit_year_start.setText(str(2015))

        self.label_model_run = qtw.QLabel('Run & save parameter scenario results as:')
        self.edit_model_run = qtw.QLineEdit()
        self.button_model_run = qtw.QPushButton('Generate results', self)
        self.button_model_run.clicked.connect(self.runParameterScenario)

        # Layout.
        grid_scen_load = qtw.QGridLayout()
        grid_scen_load.setSpacing(10)

        grid_scen_load.addWidget(self.label_parset, 0, 0)
        grid_scen_load.addWidget(self.combo_parset, 0, 1)
        grid_scen_load.addWidget(self.label_scen, 1, 0)
        grid_scen_load.addWidget(self.combo_scen, 1, 1)

        grid_scen_load.addWidget(self.label_year_start, 2, 0)
        grid_scen_load.addWidget(self.edit_year_start, 2, 1)

        grid_scen_save = qtw.QGridLayout()
        grid_scen_save.setSpacing(10)

        grid_scen_save.addWidget(self.label_model_run, 0, 0)
        grid_scen_save.addWidget(self.edit_model_run, 0, 1)
        grid_scen_save.addWidget(self.button_model_run, 0, 2)

        self.parscen_layout = qtw.QGridLayout()

        self.scroll_pars = qtw.QWidget()
        self.scroll_pars.setLayout(self.parscen_layout)

        self.scroll_area = qtw.QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setWidget(self.scroll_pars)

        layout.addLayout(grid_scen_load)
        layout.addWidget(self.scroll_area)
        layout.addLayout(grid_scen_save)


    # Updates all options-related widgets to display values from the options dictionary.
    # Generally should only be called when a default options dictionary is initialised.
    def refreshParsWidgets(self):

        # Clear out all widgets in the paraemters layout.
        for i in reversed(range(self.parscen_layout.count())):
            self.parscen_layout.itemAt(i).widget().setParent(None)

        for p, val in enumerate(self.widget_pars_list):
            par_label, par_name = val
            label_par = qtw.QLabel(par_name)
            edit_par = qtw.QLineEdit()
            self.parscen_layout.addWidget(label_par, p, 0)
            self.parscen_layout.addWidget(edit_par, p, 1)
            widget = self.parscen_layout.itemAtPosition(p, 1).widget()
            widget.setText(str(self.cascade_pars[par_label]))

        self.edit_year_start.setText(str(self.start_year))

    def refreshParsetComboBox(self):
        self.combo_parset_dict = {}
        self.combo_parset.clear()
        pid = 0
        for parset_name in self.project.parsets:
            self.combo_parset_dict[parset_name] = pid
            self.combo_parset.addItem(parset_name)
            pid += 1
        try: self.combo_parset.setCurrentIndex(self.combo_parset_dict[self.parset_name])
        except: pass

    def refreshScenComboBox(self):
        self.combo_scen_dict = ['Cascade']
        self.combo_scen.clear()
        for scen_name in self.combo_scen_dict:
            self.combo_scen.addItem(scen_name)
        self.combo_scen.setCurrentIndex(0)

    def loadCalibration(self, parset_name, delay_refresh=False):
        self.parset_name = str(parset_name)
        self.parset = dcp(self.project.parsets[self.parset_name])
        self.status = ('Status: Parameter set "%s" selected for parameter scenario' % self.parset_name)
        self.refreshParsWidgets()
        if not delay_refresh:
            self.refreshVisibility()

#    def loadScentype(self, scentype, delay_refresh=False): # CK: for setting different scenario types


    # Whatever auxiliary functions you require to link to widgets and so on.
    # E.g...
    def translateToParameterScenario(self):
        """
        
        Assumes self.project has been set. 

        # NOTE: Currently uses dummy values for diagnosis / linkage / treatment. Uncomment when values are linked in from GUI
        
        """
        # translate params in from scenario_dict
        # NOTE: Currently uses dummy values. Uncomment when values are linked in from GUI
        diag_param = self.cascade_pars['param_diag']
        link_param = self.cascade_pars['param_link']
        succ_param = self.cascade_pars['param_succ']
        fail_param = self.cascade_pars['param_fail']
        year_param = self.start_year

        # setup populations
        pops = self.project.data['pops']['label_names'].keys()

        # setup params
        active_succ_bases = ['s%s%ssuc_rate']
        active_link_bases = ['s%s%syes_rate']
        active_fail_bases = ['s%s%sno_rate']
        active_diag_bases = ['s%s%sdiag_rate']

        smear = ['p', 'n']
        strain = ['d', 'm', 'x']

        active_succ_params = [base % (sm, st) for base in active_succ_bases for sm in smear for st in strain]
        active_link_params = [base % (sm, st) for base in active_link_bases for sm in smear for st in strain]
        active_fail_params = [base % (sm, st) for base in active_fail_bases for sm in smear for st in strain]
        active_diag_params = [base % (sm, st) for base in active_diag_bases for sm in smear for st in strain]

        # setup year
        dt = self.project.settings.tvec_dt
        year_end = self.project.settings.tvec_end

        scenarios_years = np.arange(year_param, year_end, dt)

        # setup scenario values structure:
        values = odict()
        values['Care Cascade'] = {'active_diag': diag_param,
                                  'active_link': link_param,
                                  'active_succ': succ_param,
                                  'active_fail': fail_param}
        # if we wanted further scenarios, include here
        # ...

        # generate scenario values dictionary
        scen_values = odict()  # complete scenario values
        for scenario in values:
            scvalues = odict()  # the core of the parameter values for each scenario

            for paramclass in values[scenario]:
                paramlist = locals()['%s_params' % paramclass]

                for param in paramlist:
                    scvalues[param] = odict()

                    for pop_label in pops:

                        scvalues[param][pop_label] = odict()
                        scvalues[param][pop_label]['y_format'] = 'fraction'
                        scvalues[param][pop_label]['y_factor'] = DO_NOT_SCALE
                        scvalues[param][pop_label]['y'] = np.ones(scenarios_years.shape) * values[scenario][paramclass]
                        scvalues[param][pop_label]['t'] = scenarios_years

            scen_values[scenario] = { 'type': 'Parameter',
                                      'overwrite' : True,
                                      'run_scenario' : True,
                                      'scenario_values': scvalues}
        return scen_values


    # Scans through widgets and updates options dict appropriately.
    # Return True if successful or False if there was an error.
    def updatePars(self):
        if self.cascade_pars is None:
            return False
        else:
            try:
                self.start_year = float(str(self.edit_year_start.text()))

                for p, val in enumerate(self.widget_pars_list):
                    par_label, par_name = val
                    widget = self.parscen_layout.itemAtPosition(p, 1).widget()
                    self.cascade_pars[par_label] = float(str(widget.text()))

            except Exception as E:
                print('Could not update parameters: %s' % E.__repr__())
                return False
            return True


    def runParameterScenario(self):
        if not self.updatePars():
            self.status = ('Status: User-specified options could not be read into parameter scenario, so options are being reverted')
            self.refreshParsWidgets()
        else:
            self.status = ('Status: Running scenario for parameter set "%s"' % self.parset_name)
            self.refreshStatus()
            result_name = str(self.edit_model_run.text())
            if result_name == '':
                result_name = None

            scen_values = self.translateToParameterScenario()
            progset_name = None # TODO WARNING is this OK?

            self.project.createScenarios(scen_values)
            self.project.runScenarios(original_parset_name=self.parset_name, original_progset_name=progset_name,
                                          scenario_set_name=result_name, include_bau=False, save_results=False)

#            self.project.runSim(parset_name=self.parset_name, progset_name=self.progset_name, options=self.options, store_results=True, result_type='scen_budget', result_name=result_name)
            self.acknowledgeResults()
            self.status = ('Status: Scenarios successfully run for parameter set "%s"' % (self.parset_name))
        self.refreshVisibility()

# %% GUI class for running budget scenarios

class GUIBudgetScenario(GUIResultPlotterIntermediate):

    def __init__(self):
        super(GUIBudgetScenario, self).__init__()
        self.initUIBudgetScenario()


    def initUIBudgetScenario(self):
        self.setWindowTitle('Budget scenario')
        self.resetAttributes()
        self.refreshVisibility()
        self.show()


    def resetAttributes(self):
        self.resetAttributesProjectManager()
        self.resetAttributesResultPlotter()

        self.parset_name = None
        self.progset_name = None
        self.combo_parset_dict = {}     # Dictionary that maps ParameterSet names to indices used in combo boxes.
        self.combo_progset_dict = {}    # Dictionary that maps ProgramSet names to indices used in combo boxes.

        self.options = None     # The options dictionary for running a budget scenario.
        self.widget_budget_dict = {}    # Dictionary that maps Program names to indices marking their position in a scrolling list of budget widgets.

        # Initialise attributes specific to your budget scenario GUI.


    def refreshVisibility(self):
        self.refreshVisibilityProjectManager()
        self.refreshVisibilityResultPlotter()

        is_project_loaded = self.project is not None
        does_parset_exist = is_project_loaded and len(self.project.parsets) > 0
        does_progset_exist = is_project_loaded and len(self.project.progsets) > 0
        do_options_exist = self.options is not None

        if is_project_loaded:
            self.refreshParsetComboBox()
            self.refreshProgsetComboBox()
        self.label_parset.setVisible(is_project_loaded)
        self.combo_parset.setVisible(is_project_loaded)
        self.label_progset.setVisible(is_project_loaded)
        self.combo_progset.setVisible(is_project_loaded)

        policy_min = qtw.QSizePolicy.Minimum
        policy_exp = qtw.QSizePolicy.Expanding
        if do_options_exist:
            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_min)
        else:
            self.process_layout_stretch.changeSize(0, 0, policy_min, policy_exp)
        self.scroll_area.setVisible(do_options_exist)

        self.label_year_start.setVisible(do_options_exist)
        self.edit_year_start.setVisible(do_options_exist)
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
        self.label_parset = qtw.QLabel('Parameter set to use: ')
        self.combo_parset = qtw.QComboBox(self)
        self.combo_parset.activated[str].connect(self.loadCalibration)

        self.label_progset = qtw.QLabel('Program set to use: ')
        self.combo_progset = qtw.QComboBox(self)
        self.combo_progset.activated[str].connect(self.loadPrograms)

        self.label_year_start = qtw.QLabel('Start year for program budgets... ')
        self.edit_year_start = qtw.QLineEdit()

        self.label_model_run = qtw.QLabel('Run & save budget scenario results as... ')
        self.edit_model_run = qtw.QLineEdit()
        self.button_model_run = qtw.QPushButton('Generate results', self)
        self.button_model_run.clicked.connect(self.runBudgetScenario)

        # Layout.
        grid_progset_load = qtw.QGridLayout()
        grid_progset_load.setSpacing(10)

        grid_progset_load.addWidget(self.label_parset, 0, 0)
        grid_progset_load.addWidget(self.combo_parset, 0, 1)
        grid_progset_load.addWidget(self.label_progset, 1, 0)
        grid_progset_load.addWidget(self.combo_progset, 1, 1)

        grid_progset_load.addWidget(self.label_year_start, 2, 0)
        grid_progset_load.addWidget(self.edit_year_start, 2, 1)

        grid_progset_save = qtw.QGridLayout()
        grid_progset_save.setSpacing(10)

        grid_progset_save.addWidget(self.label_model_run, 0, 0)
        grid_progset_save.addWidget(self.edit_model_run, 0, 1)
        grid_progset_save.addWidget(self.button_model_run, 0, 2)

        self.budget_layout = qtw.QGridLayout()

        self.scroll_budgets = qtw.QWidget()
        self.scroll_budgets.setLayout(self.budget_layout)

        self.scroll_area = qtw.QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setWidget(self.scroll_budgets)


        layout.addLayout(grid_progset_load)
        layout.addWidget(self.scroll_area)
        layout.addLayout(grid_progset_save)


    # Updates all options-related widgets to display values from the options dictionary.
    # Generally should only be called when a default options dictionary is initialised.
    def refreshOptionWidgets(self):

        # Clear out all widgets in the budget layout.
        # TODO: Make more efficient by only clearing when absolutely necessary.
        for i in reversed(range(self.budget_layout.count())):
            self.budget_layout.itemAt(i).widget().setParent(None)

        self.widget_budget_dict = {}
        for prog_label in self.options['init_alloc']:
            prog_name = self.project.data['meta']['progs']['label_names'][prog_label]
            if prog_name not in self.widget_budget_dict:
                try: last_id = max(self.widget_budget_dict.values())    # TODO: Make more efficient by storing max id rather than calculating all the time.
                except: last_id = -1
                label_budget = qtw.QLabel(prog_name)
                edit_budget = qtw.QLineEdit()
                self.budget_layout.addWidget(label_budget, last_id + 1, 0)
                self.budget_layout.addWidget(edit_budget, last_id + 1, 1)
                self.widget_budget_dict[prog_name] = last_id + 1

            current_id = self.widget_budget_dict[prog_name]
            widget = self.budget_layout.itemAtPosition(current_id, 1).widget()
            widget.setText(str("%.0f" % self.options['init_alloc'][prog_label]))

        self.edit_year_start.setText(str(self.options['progs_start']))

    def refreshParsetComboBox(self):
        self.combo_parset_dict = {}
        self.combo_parset.clear()
        pid = 0
        for parset_name in self.project.parsets:
            self.combo_parset_dict[parset_name] = pid
            self.combo_parset.addItem(parset_name)
            pid += 1
        try: self.combo_parset.setCurrentIndex(self.combo_parset_dict[self.parset_name])
        except: pass

    def refreshProgsetComboBox(self):
        self.combo_progset_dict = {}
        self.combo_progset.clear()
        pid = 0
        for progset_name in self.project.progsets:
            self.combo_progset_dict[progset_name] = pid
            self.combo_progset.addItem(progset_name)
            pid += 1
        try: self.combo_progset.setCurrentIndex(self.combo_progset_dict[self.progset_name])
        except: pass

    def loadCalibration(self, parset_name, delay_refresh=False):
        self.parset_name = str(parset_name)
        self.status = ('Status: Parameter set "%s" selected for budget scenario' % self.parset_name)
        if not delay_refresh:
            self.refreshVisibility()

    def loadPrograms(self, progset_name, delay_refresh=False):
        self.progset_name = str(progset_name)
        self.options = defaultOptimOptions(settings=self.project.settings, progset=self.project.progsets[self.progset_name])
        self.refreshOptionWidgets()
        self.status = ('Status: Program set "%s" selected for budget scenario' % self.progset_name)
        if not delay_refresh:
            self.refreshVisibility()

    def runBudgetScenario(self):
        if not self.updateOptions():
            self.status = ('Status: User-specified options could not be read into budget scenario, so options are being reverted')
            self.refreshOptionWidgets()
        else:
            self.status = ('Status: Running model for Parset "%s" and Progset "%s"' % (self.parset_name, self.progset_name))
            self.refreshStatus()
            result_name = str(self.edit_model_run.text())
            if result_name == '':
                result_name = None
            self.project.runSim(parset_name=self.parset_name, progset_name=self.progset_name, options=self.options, store_results=True, result_type='scen_budget', result_name=result_name)
            self.acknowledgeResults()
            self.status = ('Status: Model successfully processed for parameter set "%s" and program set "%s"' % (self.parset_name, self.progset_name))
        self.refreshVisibility()

    # Scans through widgets and updates options dict appropriately.
    # Return True if successful or False if there was an error.
    def updateOptions(self):
        if self.options is None:
            return False
        else:
            try:
                self.options['progs_start'] = float(str(self.edit_year_start.text()))

                for prog_name in self.widget_budget_dict:
                    current_id = self.widget_budget_dict[prog_name]
                    widget = self.budget_layout.itemAtPosition(current_id, 1).widget()
                    prog_label = self.project.data['meta']['progs']['name_labels'][prog_name]
                    self.options['init_alloc'][prog_label] = float(str(widget.text()))

            except: return False
            return True



# %% This is the thing that's actually called

def runGUI():
    ''' Function that launches all available back-end GUIs as they are developed. '''
    app = qtw.QApplication(sys.argv)
    app.setApplicationName('Optima GUI')
    gui = GUI()
    sys.exit(app.exec_())
    return gui # Not super needed, but avoids pylint warning
