# %% Imports
from optima_tb.utils import odict
from optima_tb.cascade import loadCascadeSettingsFunc, plotCascadeFunc


import logging
logger = logging.getLogger(__name__)

import pylab as pl
import numpy as np
from matplotlib.ticker import FuncFormatter



# %% Settings class (for data that is effectively static per epidemic context)

VALIDATION_IGNORE = 0
VALIDATION_WARN = 1
VALIDATION_ERROR = 2
VALIDATION_AVERT = 3

DEFAULT_YFACTOR = 1.
DO_NOT_SCALE = -1.

TOLERANCE = 1e-6

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
        self.tvec_end = 2030.0  # Default end year for data input and simulations.
        self.tvec_observed_end = 2015.0
        self.tvec_dt = 1.0 / 4  # Default timestep for simulations.

        self.recursion_limit = 100  # Limit for recursive references, primarily used in avoiding circular references for definitions using dependencies.
        self.fit_metric = 'wape'

        self.parser_debug = False  # Decomposes and evaluates functions written as strings, in accordance with a grammar defined within the parser object.

        # Separate autofit and optimization parameters so that can use the asd algorithm
        # with different settings.
        self.autofit_params = self.resetCalibrationParameters()
        self.optimization_params = self.resetCalibrationParameters()

        # Settings for databooks / spreadsheets / workbooks / cascades:
        self.loadCascadeSettings(cascade_path)
        self.initCustomDatabookFramework()  # Creates self.countrybook.
                                                        # NOTE: Databook will hopefully one day be capable of replicating countrybook, making the two different
                                                        # unnecessary.

        logging.info("Created settings based on cascade: %s" % cascade_path)



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
        which is used for autofitting the calibration and running optimizations
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
        calibration['useInitCompartments'] = True

        return calibration

    def loadCascadeSettings(self, cascade_path):
        ''' Generates node and link settings based on cascade spreadsheet. '''
        self.resetCascade()
        loadCascadeSettingsFunc(cascade_path, settings=self)

    def plotCascade(self):
        ''' Plots cascade network. '''
        plotCascadeFunc(settings=self)



    def initCustomDatabookFramework(self):
        """
        Settings for country data book
        """

        self.countrybook = odict()
        self.countrybook['sheet_names'] = odict()
        self.countrybook['sheet_names']['populations'] = 'Populations'
        self.countrybook['sheet_names']['demographics'] = 'Demographics'
        self.countrybook['sheet_names']['prevalence'] = 'TB incidence & prevalence'
        self.countrybook['sheet_names']['notifications'] = 'Notifications'
        self.countrybook['sheet_names']['smear'] = 'Smear status'
        self.countrybook['sheet_names']['comorbidity'] = 'Comorbidities'
        self.countrybook['sheet_names']['morbidity'] = 'TB morbidity'
        self.countrybook['sheet_names']['testing_treatment_latent'] = 'Outcomes - Latent TB'
        self.countrybook['sheet_names']['testing_treatment_active'] = 'Outcomes - Active TB'
        # self.countrybook['sheet_names']['programs'] = 'Programs'
        # self.countrybook['sheet_names']['cost_coverage'] = 'Cost and coverage'
        # self.countrybook['sheet_names']['unitcost'] = 'Unit costs'
        # self.countrybook['sheet_names']['poptransitions'] = 'Population transitions'

        # headers for special sheets (i.e. headers aren't years)
        self.countrybook['headers'] = odict()
        self.countrybook['headers']['populations'] = ['Name', 'Minimum Age', 'Maximum Age']
        # self.countrybook['headers']['programs'] = ['Name', 'Short name', 'Intervention class', 'Coverage indicator (annual)', 'Duration of treatment (days per person on average)', 'Frequency of intervention (in years)']

        # labels for each sheet
        self.countrybook['labels'] = {'populations': '',
                                      'demographics':'',
                                      'prevalence'     : 'Estimated incidence and prevalence of TB',
                                      'notifications'  : 'TB notifications: please enter the number of notifications of active TB per population and drug-resistant strain per year',
                                      'smear'          : '',
                                      'comorbidity':'Comorbidities',
                                      'morbidity' : 'Number of TB-related deaths per year',
                                      'testing_treatment_latent': 'Testing and treatment outcomes for latent TB',
                                      'testing_treatment_active'        : 'Testing and treatment outcomes for latent TB'
                                      }

        # info
        self.countrybook['disaggregations'] = odict()
        # other values
        self.countrybook['disaggregations']['strains'] = ['DS-TB', 'MDR-TB', 'XDR-TB']
        self.countrybook['disaggregations']['smears'] = ['Smear-', 'Smear+']  # also potentially 'Extrapulmonary'
        # disagg for non-total smear and strain
        self.countrybook['disaggregations']['nt_strains'] = ['DS-TB', 'MDR-TB', 'XDR-TB']
        self.countrybook['disaggregations']['nt_smears'] = ['Smear-', 'Smear+']  # also potentially 'Extrapulmonary'
        self.countrybook['disaggregations']['populations'] = []  # determined dynamically at runtime
        self.countrybook['disaggregations']['regimens'] = ['DS-TB regimen', 'MDR-TB regimen', 'XDR-TB regimen']
        self.countrybook['disaggregations']['programs'] = []  # determined dynamically at runtime
        self.countrybook['disaggregations']['total_pop'] = ['Total population']
        self.countrybook['disaggregations']['total_birth'] = ['Total number of births']

        # for univalue sheets, includes information on how data should be disaggregated
        self.countrybook['sheet_classes'] = odict()
        self.countrybook['sheet_classes']['univalue'] = odict()
#         self.countrybook['sheet_classes']['univalue']['population_sizes'] = ['populations']
#         self.countrybook['sheet_classes']['univalue']['total_cases'] = ['populations', 'smears', 'strains']
        self.countrybook['sheet_classes']['univalue']['notifications'] = ['populations', 'smears', 'strains']
        self.countrybook['sheet_classes']['univalue']['morbidity'] = ['populations', 'strains']

        # sheet specific values
        self.countrybook['sheet_values'] = odict()

        self.countrybook['sheet_values']['demographics'] = odict()
        self.countrybook['sheet_values']['demographics']['Population size per year'] = ['populations']
        self.countrybook['sheet_values']['demographics']['Number of births per year'] = ['total_birth']
        self.countrybook['sheet_values']['demographics']['Percentage of people vaccinated per year'] = ['populations']
        self.countrybook['sheet_values']['demographics']['Percentage of people who die from non-TB related causes per year'] = ['populations']

        self.countrybook['sheet_values']['prevalence'] = odict()
        self.countrybook['sheet_values']['prevalence']['Estimated TB incidence (per 100,000)'] = ['populations']
        self.countrybook['sheet_values']['prevalence']['Estimated active TB prevalence'] = ['populations']
        self.countrybook['sheet_values']['prevalence']['Estimated MDR TB prevalence'] = ['populations']
        self.countrybook['sheet_values']['prevalence']['Estimated XDR TB prevalence'] = ['populations']
        self.countrybook['sheet_values']['prevalence']['Estimated latent TB prevalence'] = ['populations']

        self.countrybook['sheet_values']['smear'] = odict()
        self.countrybook['sheet_values']['smear']['Smear status by drug-resistant strain'] = ['nt_smears', 'nt_strains']
        self.countrybook['sheet_values']['smear']['Smear status by population'] = ['nt_smears', 'populations']


        self.countrybook['sheet_values']['comorbidity'] = odict()
#         self.countrybook['sheet_values']['comorbidity']['HIV prevalence'] = ['populations', 'smears', 'strains']
#         self.countrybook['sheet_values']['comorbidity']['Diabetes prevalence'] = ['populations', 'smears', 'strains']

        self.countrybook['sheet_values']['testing_treatment_latent'] = odict()
        self.countrybook['sheet_values']['testing_treatment_latent']['Percentage of population tested for latent TB per year'] = ['populations']
        self.countrybook['sheet_values']['testing_treatment_latent']['Number of people initiating treatment for latent TB per year'] = ['populations']
        self.countrybook['sheet_values']['testing_treatment_latent']['Number of people lost to follow up for latent TB per year'] = ['populations']
        self.countrybook['sheet_values']['testing_treatment_latent']['Number of people who succesfully completed treatment for latent TB'] = ['populations']

        self.countrybook['sheet_values']['testing_treatment_active'] = odict()
        self.countrybook['sheet_values']['testing_treatment_active']['Percentage of population tested for active TB per year'] = ['populations']
        self.countrybook['sheet_values']['testing_treatment_active']['Number of people initiating treatment for active TB per year'] = ['regimens', 'populations']
        self.countrybook['sheet_values']['testing_treatment_active']['Number of people lost to follow up for active TB per year'] = ['regimens', 'populations']
        self.countrybook['sheet_values']['testing_treatment_active']['Number of people who failed treatment for active TB'] = ['regimens', 'populations']
        self.countrybook['sheet_values']['testing_treatment_active']['Number of people who successfully completed treatment for active TB'] = ['regimens', 'populations']


        self.countrybook['constants'] = {'spacing_interpopulation':2,
                                         'spacing_intrapopulation':1,
                                         'spacing_interproperty'  :4,
                                         'spacing_multivalue_label':2,
                                         'total_strains': 'Total',  # All strains',
                                         'total_smears' : 'Total',
                                         'num_default_programs':28,
                                         'row_index_start':2,  # for when there are no disaggregations, etc.
                                         'col_index_start':1}  #

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
        return ['negative_population',  # runs @ validation/checkNegativePopulation
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
            elif x >= 1e3:
                return '%1.1fK' % (x * 1e-3)
            else:
                return x

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
                         'formatter': None, # FuncFormatter(self.PopSuffixFormatter), # KMSuffixFormatter) ,
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
                         'use_full_labels' : False,
                         'effective_rate': "Effective number",
                         'default_ylabel' : "Number of cases",
                         'default_pops' : "All populations",
                         'default_figname' : "Figname",
                         'default_year' : 2016., # TODO change so that it's last year,
                         'title' : None,
                         # relative plot settings
                         'relative_yticks' : np.arange(0., 1.1, 0.2)
                         }

    def devSettings(self):
        pl.rcParams['figure.figsize'] = (10, 8)
        pl.rcParams['savefig.dpi'] = 100
        pl.rcParams['savefig.transparent'] = 'False'  # relax

    def printSettings(self):

        pl.rcParams['figure.figsize'] = (15, 10)
        pl.rcParams['savefig.dpi'] = 200
        pl.rcParams['savefig.transparent'] = 'True'  # enforce
        self.plotdict['legend_off'] = False
        self.plotdict['use_full_labels'] = True

    def presentationSettings(self):
        pl.rcParams['font.size'] = 14
        pl.rcParams['figure.figsize'] = (9, 7)
        pl.rcParams['savefig.dpi'] = 300

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

        self.plotdict['legend_off'] = True
        self.plotdict['title'] = ''  # No title when we have presentation quality
        self.plotdict['num_cols'] = 1
        self.plotdict['use_full_labels'] = True

    def guiSettings(self):

        pl.rcParams['lines.linewidth'] = 2
        pl.rcParams['axes.linewidth'] = 1.5
        pl.rcParams['axes.labelsize'] = pl.rcParams['font.size']
        pl.rcParams['axes.titlesize'] = pl.rcParams['font.size']

        pl.rcParams['legend.frameon'] = False

        self.plotdict['legendsettings'] = {'loc': 0 } # best ##{'loc': 4 } # lower right
        self.plotdict['box_width'] = 0.9
        self.plotdict['box_offset'] = 0.1
        self.plotdict['use_full_labels'] = True




