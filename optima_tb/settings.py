# %% Imports
from optima_tb.utils import odict
from optima_tb.cascade import loadCascadeSettingsFunc, plotCascadeFunc


import logging
logger = logging.getLogger(__name__)

import pylab as pl
import numpy as np
from matplotlib.ticker import FuncFormatter
import optima_tb.plotting2 as oplt


# %% Settings class (for data that is effectively static per epidemic context)

VALIDATION_IGNORE = 0
VALIDATION_WARN = 1
VALIDATION_ERROR = 2
VALIDATION_AVERT = 3

DEFAULT_YFACTOR = 1.
DO_NOT_SCALE = -1.

TOLERANCE = 1e-2

class Settings(object):
    '''
    An object for storing cascade metadata (loaded from a cascade workbook) and general defaults.
    A cascade structure comprises of compartments/nodes and transitions/links that detail the flow of people through stages of a disease.
    In general use, the attributes of this object are static regardless of project/simulation.

    Notes:
    Settings must be loaded into a project during its creation, as it details the framework for dealing with data and running the model.
    '''

    def __init__(self, cascade_path='../data/cascade.xlsx', **args):
        """
        
        
        """

        logging.info("Loading settings")

        try:
            self.validation = ValidationSettings(args['validation_level']).settings
        except:
            self.validation = ValidationSettings().settings

        try:
            self.plot_settings = PlottingSettings(args['plotting_level']).plotdict
        except:
            self.plot_settings = PlottingSettings().plotdict


        self.tvec_start = 2000.0  # Default start year for data input and simulations.
        self.tvec_end = 2035.0  # Default end year for data input and simulations.
        self.tvec_observed_end = 2015.0
        self.tvec_dt = 1.0 / 4  # Default timestep for simulations.

        self.recursion_limit = 100  # Limit for recursive references, primarily used in avoiding circular references for definitions using dependencies.
        self.fit_metric = 'wape'

        self.parser_debug = False  # Decomposes and evaluates functions written as strings, in accordance with a grammar defined within the parser object.

        # Separate autofit and optimization parameters so that can use the asd algorithm
        # with different settings.
        self.autofit_params = self.resetCalibrationParameters()
        self.optimization_params = self.resetOptimizationParameters()
        self.parallel_optimization_params = self.resetOptimizationParameters()
        # Settings for databooks / spreadsheets / workbooks / cascades:
        self.loadCascadeSettings(cascade_path)

        logging.info("Created settings based on cascade: %s" % cascade_path)

    @property
    def tvec(self):
        return np.arange(self.tvec_start, self.tvec_end + self.tvec_dt / 2, self.tvec_dt)

    def resetCascade(self):
        ''' Resets all cascade contents and settings that are fundamental to how a project is structured. '''

        self.node_specs = odict()  # Relates compartment code-labels with special tags, if they exist.
                                                # Key is a node label. Value is a dict including a 'dead' tag and networkx-related information.
        self.node_names = []  # A corresponding list of full names for compartments.
        self.junction_labels = []  # A list of labels for compartments for which inflows must immediately propagated as outflows.

        self.charac_specs = odict()  # Relates code-labels for defined characteristics (e.g. prevalence) with labels of compartments used in their definition.
                                                # Key is a characteristic label. Value is a dict containing characteristic name, a list of 'inclusions' and a normalising characteristic or compartment.
                                                # Be aware that inclusions/normalisations may refer to characteristics in the same odict.
        self.charac_name_labels = odict()  # Key is a characteristic name. Value is a characteristic label. (A partial reversed charac_specs.)

        self.charac_pop_count = 'auto_pop_count'  # The label for a characteristic that includes all 'transfer-enabled' compartments. Is overwritten if one is user-defined.
        self.charac_pop_count_name = 'Standard Compartment Sum'  # The name for this 'transfer-enabled' population-count characteristic.

        self.links = odict()  # Key is a tag. Value is a list of compartment-label tuple.
        self.linkpar_specs = odict()  # Key is a link-parameter label. Value is a dict including link tag, link-parameter name, default value.
        self.linkpar_name_labels = odict()  # Key is a link-parameter name. Value is a link-parameter label. (A partial reversed linkpar-specs.)
        self.linkpar_outputs = []  # A list of parameters that should be returned as result outputs along with characteristics.

        self.par_funcs = odict()  # A definition-ordered dictionary of parameters that are functionally defined in terms of other parameters/characteristics.
        self.par_deps = odict()  # A definition-ordered dictionary of 'untagged' parameters used as dependencies for other parameters.
        self.charac_deps = {}  # An unordered dictionary of characteristics that must be calculated at each model timestep due to being dependencies for other variables.
                                                # Should correspond to every item in charac_specs that has a 'par_dependency' tag.

        self.progtype_specs = odict()  # Relates program type code-labels with impact parameters, etc.
        self.progtype_name_labels = odict()  # Key is a program type name. Value is a program type label.

        # Project-data workbook metadata.
        self.databook = odict()
        self.databook['sheet_names'] = odict()
        self.databook['sheet_names']['pops'] = 'Population Definitions'
        self.databook['sheet_names']['contact'] = 'Population Contacts'
        self.databook['sheet_names']['transmat'] = 'Transfer Definitions'
        self.databook['sheet_names']['transval'] = 'Transfer Details'
        self.databook['sheet_names']['progmat'] = 'Program Definitions'
        self.databook['sheet_names']['progval'] = 'Program Details'
        self.databook['sheet_names']['charac'] = 'Epidemic Characteristics'
        self.databook['sheet_names']['linkpars'] = 'Cascade Parameters'
        self.databook['custom_sheet_names'] = odict()

        self.make_sheet_characs = True  # Tag for whether the default characteristics worksheet should be generated during databook creation.
        self.make_sheet_linkpars = True  # Tag for whether the default cascade parameters worksheet should be generated during databook creation.

        self.databook['suffix'] = odict()
        self.databook['suffix']['seed'] = ' [S]'  # Suffix for characteristics used as model seeds (i.e. for initialisation).
        self.databook['suffix']['output'] = ' [O]'  # Suffix for characteristics used solely as outputs for diagnostic and/or calibration purposes.
        self.databook['suffix']['par'] = ' [P]'  # Suffix for parameters that are used at every step of model calculations.

        self.databook['format'] = {'programs':{'max_lines_impact':0}}

    def resetCalibrationParameters(self):
        """
        Sets calibration parameters for use in ASD algorithm, 
        which is used for autofitting the calibration.
        For full list of calibration parameters, see asd.py > asd() function signature.
        """
        calibration = odict()
        calibration['stepsize'] = 0.1
        calibration['maxiters'] = 2000
        calibration['maxtime'] = 300.  # Time in seconds.

        calibration['sinc'] = 1.5
        calibration['sdec'] = 2.
        calibration['fulloutput'] = False

        calibration['useYFactor'] = True
#         calibration['useInitCompartments'] = True

        return calibration

    def resetOptimizationParameters(self):
        """
        Sets calibration parameters for use in ASD algorithm, 
        which is used for running optimizations
        For full list of calibration parameters, see asd.py > asd() function signature.
        """
        calibration = odict()
        calibration['stepsize'] = 0.1
        calibration['maxiters'] = 2000
        calibration['maxtime'] = 3600.  # Time in seconds.
        calibration['reltol'] = 1e-3
        calibration['sinc'] = 2.
        calibration['sdec'] = 2.
        calibration['fulloutput'] = True
        calibration['num_threads'] = 4
        calibration['block_iter'] = 10
        calibration['max_blocks'] = 10
#         calibration['useYFactor'] = True
#         calibration['useInitCompartments'] = True

        return calibration

    def loadCascadeSettings(self, cascade_path):
        ''' Generates node and link settings based on cascade spreadsheet. '''
        self.resetCascade()
        loadCascadeSettingsFunc(cascade_path, settings=self)

    def plotCascade(self,code_labels=True):
        ''' Plots cascade network. '''
        plotCascadeFunc(settings=self,code_labels=code_labels)

class ValidationSettings():
    """
    For each of the validation checks returned in getValidationTypes(),
    there is a validation setting of either
        IGNORE : do nothing
        AVERT  : try to modify values based on a reasonable assumption
        WARN   : continue performing the action, but notify occurrence of incorrect usage via a warn statement
        ERROR  : stop the action
        
    The exact process will vary on where and what the validation check is acting on, and should be noted
    in the method's usage. 
    
    Note that all validation checks are set to a level (default="warn"), but 
    that it is possible to set specific levels for different validation checks.
    
    Examples:
    
        from optima_tb.settings import ValidationSettings, VALIDATION_AVERT
        validation_settings = ValidationSettings(level="warn")
        print validation_settings.getValidationTypes()
        > ['negative_population', 'transition_fraction']
        validation_settings['transition_fraction'] = VALIDATION_AVERT
        
    """


    def __init__(self, level='warn'):

        self.defaultSettings()
        try:
            validationSettings = getattr(self, '%sSettings' % level)
            validationSettings()
            logging.info("Validation settings: %s" % level)

        except:
            logging.info("Could not load validation settings for level: %s" % level)


    def getValidationTypes(self):
        return ['negative_population',
                'transition_fraction',  # runs @ validation/checkTransitionFraction
                'databook_validation',  # runs @ databook/loadSpreadsheetFunc
                ]

    def defaultSettings(self):

        self.settings = odict()
        self.warnSettings()

    def setValidationLevel(self, validation_level):
        for v in self.getValidationTypes():
            self.settings[v] = validation_level

    def errorSettings(self):
        self.setValidationLevel(VALIDATION_ERROR)

    def ignoreSettings(self):
        self.setValidationLevel(VALIDATION_IGNORE)

    def warnSettings(self):
        self.setValidationLevel(VALIDATION_WARN)

    def avertSettings(self):
        self.setValidationLevel(VALIDATION_AVERT)

class PlottingSettings():



    def __init__(self, output='dev'):

        logging.info("Loading plotting settings: %s" % output)


        self.plotdict = {}  # holder
        self.defaultSettings()
        try:
            outputSettings = getattr(self, '%sSettings' % output)
            outputSettings()
        except:
            logging.info("Could not load rcParams for plotting for output: %s" % output)

    def KMSuffixFormatter(self, x, pos):
            'The two args are the value and tick position'
            if x >= 1e6:
                return '%1.1fM' % (x * 1e-6)
#             elif x >= 1e3:
#                 return '%1.1fK' % (x * 1e-3)
            else:
                return '%g' % x

    def PopSuffixFormatter(self, x, pos):
            'The two args are the value and tick position'
            if x >= 1e6:
                return '%1.1fM' % (x * 1e-6)
            elif x >= 1e3:
                return '%1.1fK' % (x * 1e-3)
            else:
                return '%i' % (x)

    def defaultSettings(self):

        pl.rcParams['font.size'] = 12
        pl.rcParams['font.family'] = 'sans-serif'

        pl.rcParams['savefig.dpi'] = 300
        pl.rcParams['savefig.format'] = 'png'
        pl.rcParams['savefig.transparent'] = 'True'

        pl.rcParams['figure.max_open_warning'] = 40

        pl.rcParams['xtick.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['xtick.major.size'] = 3
        pl.rcParams['xtick.minor.size'] = 3
        pl.rcParams['xtick.major.width'] = 1
        pl.rcParams['xtick.minor.width'] = 1
        pl.rcParams['ytick.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['ytick.major.size'] = 3
        pl.rcParams['ytick.minor.size'] = 3
        pl.rcParams['ytick.major.width'] = 1
        pl.rcParams['ytick.minor.width'] = 1

        pl.rcParams['legend.frameon'] = False
        pl.rcParams['legend.loc'] = 'center left'
        pl.rcParams['legend.fontsize'] = pl.rcParams['font.size']

        pl.rcParams['axes.linewidth'] = 2
        pl.rcParams['axes.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['axes.titlesize'] = 1.5 * pl.rcParams['font.size']

        pl.rcParams['lines.linewidth'] = 3
        pl.rcParams['lines.marker'] = 'None'

        # The following requires matplotlib 2.X
#         pl.rcParams['hatch.linewidth'] = 1.

        # Non-standard list of parameters used in plotting
        self.plotdict = {  # scatter plotting values
                         'marker' : 'o',
                         'facecolors' : 'none',
                         'marker_color': 'k',
                         'hatch_bg' : 'white',
                         's' : 40,
                         # axes format
                         'year_inc':5,
                         # colormapping for category lists
                         'colormapping_order':'alternate3',  # as we have triplets in undiagnosed --> diagnosed --> on treatment
                         'formatter': FuncFormatter(self.KMSuffixFormatter) ,
                         'barwidth': 0.8,
                         'bar_offset': 0.2,
                         # alpha for fill-between
                         'alpha': 0.3,
                         # linestyle to be used as default
                         'default_linestyle' : '-',
                         # box width of plot and offset
                         'box_width' : 0.8,
                         'box_offset' : 0.,
                         # legend
                         'legend_off' : False, # I am legend
                         'legendsettings': {'loc':'center left',
                                            'bbox_to_anchor':(1.05, 0.5),
                                            'ncol':1},
                         # labels
                         'use_full_labels' : True,
                         'effective_rate': "Effective number",
                         'default_ylabel' : "Number of cases",
                         'default_pops' : "All populations",
                         'default_figname' : "Figname",
                         'default_year' : 2016., # TODO change so that it's last year,
                         # 'title' : None,
                         # relative plot settings
                         'relative_yticks' : np.arange(0., 1.1, 0.2),
                         'y_intercept_line': {'colors': '#AAAAAA',
                                         'linewidth': 3.,
                                         'alpha': 0.5,
                                         'linestyle':':'}
                         }

    def devSettings(self):
        pl.rcParams['figure.figsize'] = (10, 8)
        pl.rcParams['savefig.dpi'] = 100
        pl.rcParams['savefig.transparent'] = 'False'  # relax

    def printSettings(self):

        pl.rcParams['figure.figsize'] = (9, 7)
        pl.rcParams['savefig.dpi'] = 200
        pl.rcParams['savefig.transparent'] = 'True'  # enforce
        pl.rcParams['font.size'] = 14
        pl.rcParams['lines.linewidth'] = 4

        oplt.settings['separate_legend'] = False

        self.plotdict['legend_off'] = False
        self.plotdict['use_full_labels'] = True
        self.plotdict['title'] = ''  # No title when we have presentation quality


        pl.rcParams['axes.linewidth'] = 3
        pl.rcParams['axes.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['axes.titlesize'] = pl.rcParams['font.size']

        pl.rcParams['xtick.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['xtick.major.size'] = 3
        pl.rcParams['xtick.minor.size'] = 0
        pl.rcParams['xtick.major.width'] = 2
        pl.rcParams['xtick.minor.width'] = 0
        pl.rcParams['ytick.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['ytick.major.size'] = 3
        pl.rcParams['ytick.minor.size'] = 0
        pl.rcParams['ytick.major.width'] = 2
        pl.rcParams['ytick.minor.width'] = 0

    def presentationSettings(self):
        pl.rcParams['font.size'] = 14
        pl.rcParams['figure.figsize'] = (9, 7)
        pl.rcParams['savefig.dpi'] = 200

        pl.rcParams['lines.linewidth'] = 5
        pl.rcParams['lines.marker'] = 'None'

        pl.rcParams['axes.linewidth'] = 3
        pl.rcParams['axes.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['axes.titlesize'] = pl.rcParams['font.size']

        pl.rcParams['xtick.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['xtick.major.size'] = 3
        pl.rcParams['xtick.minor.size'] = 0
        pl.rcParams['xtick.major.width'] = 2
        pl.rcParams['xtick.minor.width'] = 0
        pl.rcParams['ytick.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['ytick.major.size'] = 3
        pl.rcParams['ytick.minor.size'] = 0
        pl.rcParams['ytick.major.width'] = 2
        pl.rcParams['ytick.minor.width'] = 0

        pl.rcParams['savefig.transparent'] = 'True'  # enforce

        oplt.settings['separate_legend'] = True

        self.plotdict['legend_off'] = True # plots separate legend
        self.plotdict['title'] = ''  # No title when we have presentation quality
        self.plotdict['num_cols'] = 1
        self.plotdict['use_full_labels'] = True

    def guiSettings(self):

        pl.rcParams['lines.linewidth'] = 4
        pl.rcParams['axes.linewidth'] = 1.5
        pl.rcParams['axes.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['axes.titlesize'] = pl.rcParams['font.size']


#         self.plotdict['legendsettings'] = {'loc': 'center left' , 'ncol':1} # best ##{'loc': 4 } # lower right
#         self.plotdict['box_width'] = 0.9
#         self.plotdict['box_offset'] = 0.1
#         self.plotdict['use_full_labels'] = True




