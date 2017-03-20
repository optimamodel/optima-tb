#%% Imports
from optima_tb.utils import odict
from optima_tb.cascade import loadCascadeSettingsFunc, plotCascadeFunc


import logging
logger = logging.getLogger(__name__)

import pylab as pl
import numpy as np
from matplotlib.ticker import FuncFormatter



#%% Settings class (for data that is effectively static per epidemic context)

VALIDATION_IGNORE = 0
VALIDATION_WARN = 1
VALIDATION_ERROR = 2
VALIDATION_AVERT = 3

DEFAULT_YFACTOR = 1.
DO_NOT_SCALE = -1.

TOLERANCE = 1e-9

class Settings(object):
    '''
    An object for storing cascade metadata (loaded from a cascade workbook) and general defaults.
    A cascade structure comprises of compartments/nodes and transitions/links that detail the flow of people through stages of a disease.
    In general use, the attributes of this object are static regardless of project/simulation.

    Notes:
    Settings must be loaded into a project during its creation, as it details the framework for dealing with data and running the model.
    '''    
    
    def __init__(self, cascade_path = '../data/cascade.xlsx', **args):
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
            
            

        self.tvec_start = 2000.0     # Default start year for data input and simulations.
        self.tvec_end = 2030.0       # Default end year for data input and simulations.
        self.tvec_observed_end = 2015.0
        self.tvec_dt = 1.0/4         # Default timestep for simulations.
        
        self.recursion_limit = 100      # Limit for recursive references, primarily used in avoiding circular references for definitions using dependencies.
        self.fit_metric = 'wape'
        
        self.parser_debug = False      # Decomposes and evaluates functions written as strings, in accordance with a grammar defined within the parser object.

        # Separate autofit and optimization parameters so that can use the asd algorithm
        # with different settings.
        self.autofit_params = self.resetCalibrationParameters()
        self.optimization_params = self.resetCalibrationParameters()
        
        # Settings for databooks / spreadsheets / workbooks / cascades:
        self.loadCascadeSettings(cascade_path)
        self.initCustomDatabookFramework()              # Creates self.countrybook.
                                                        # NOTE: Databook will hopefully one day be capable of replicating countrybook, making the two different
                                                        # unnecessary.
                                                         
        logging.info("Created settings based on cascade: %s"%cascade_path)
        
        
        
    def resetCascade(self):
        ''' Resets all cascade contents and settings that are fundamental to how a project is structured. '''
            
        self.node_specs = odict()               # Relates compartment code-labels with special tags, if they exist. 
                                                # Key is a node label. Value is a dict including a 'dead' tag and networkx-related information.
        self.node_names = []                    # A corresponding list of full names for compartments.
        self.junction_labels = []               # A list of labels for compartments for which inflows must immediately propagated as outflows.
        
        self.charac_specs = odict()             # Relates code-labels for defined characteristics (e.g. prevalence) with labels of compartments used in their definition.
                                                # Key is a characteristic label. Value is a dict containing characteristic name, a list of 'inclusions' and a normalising characteristic or compartment.
                                                # Be aware that inclusions/normalisations may refer to characteristics in the same odict.
        self.charac_name_labels = odict()       # Key is a characteristic name. Value is a characteristic label. (A partial reversed charac_specs.)
        
        self.charac_pop_count = 'auto_pop_count'                 # The label for a characteristic that includes all 'transfer-enabled' compartments. Is overwritten if one is user-defined.
        self.charac_pop_count_name = 'Standard Compartment Sum'  # The name for this 'transfer-enabled' population-count characteristic.
        
        self.links = odict()                    # Key is a tag. Value is a list of compartment-label tuple.
        self.linkpar_specs = odict()            # Key is a link-parameter label. Value is a dict including link tag, link-parameter name, default value.
        self.linkpar_name_labels = odict()      # Key is a link-parameter name. Value is a link-parameter label. (A partial reversed linkpar-specs.)
        
        self.par_funcs = odict()                # A definition-ordered dictionary of parameters that are functionally defined in terms of other parameters/characteristics.
        self.par_deps = odict()                 # A definition-ordered dictionary of 'untagged' parameters used as dependencies for other parameters.
        self.charac_deps = {}                   # An unordered dictionary of characteristics that must be calculated at each model timestep due to being dependencies for other variables.
                                                # Should correspond to every item in charac_specs that has a 'par_dependency' tag.
        
        self.progtype_specs = odict()           # Relates program type code-labels with impact parameters, etc.
        self.progtype_name_labels = odict()     # Key is a program type name. Value is a program type label.
        
        # Project-data workbook metadata.
        self.databook = odict()
        self.databook['sheet_names'] = odict()
        self.databook['sheet_names']['pops'] =      'Population Definitions'
        self.databook['sheet_names']['contact'] =   'Population Contacts'
        self.databook['sheet_names']['transmat'] =  'Transfer Definitions'
        self.databook['sheet_names']['transval'] =  'Transfer Details'
        self.databook['sheet_names']['progmat'] =  'Program Definitions'
        self.databook['sheet_names']['progval'] =  'Program Details'
        self.databook['sheet_names']['charac'] =    'Epidemic Characteristics'
        self.databook['sheet_names']['linkpars'] =  'Cascade Parameters'
        self.databook['custom_sheet_names'] = odict()
        
        self.make_sheet_characs = True      # Tag for whether the default characteristics worksheet should be generated during databook creation.
        self.make_sheet_linkpars = True     # Tag for whether the default cascade parameters worksheet should be generated during databook creation.
        
        self.databook['suffix'] = odict()        
        self.databook['suffix']['seed'] =   ' [S]'  # Suffix for characteristics used as model seeds (i.e. for initialisation).
        self.databook['suffix']['output'] = ' [O]'  # Suffix for characteristics used solely as outputs for diagnostic and/or calibration purposes.
        self.databook['suffix']['par'] =    ' [P]'  # Suffix for parameters that are used at every step of model calculations.
        
        self.databook['format'] = {'programs':{'max_lines_impact':0}}
    
    def resetCalibrationParameters(self):
        """
        Sets calibration parameters for use in ASD algorithm, 
        which is used for autofitting the calibration and running optimizations
        For full list of calibration parameters, see asd.py > asd() function signature.
        """
        calibration = odict()
        calibration['stepsize'] = 0.1
        calibration['MaxIter'] = 2000
        calibration['timelimit'] = 300.     # Time in seconds.
        
        calibration['sinc'] = 1.5
        calibration['sdec'] = 2.
        calibration['fulloutput'] = False
        
        calibration['useYFactor'] = True
        calibration['useInitCompartments'] = True
        
        return calibration

    def loadCascadeSettings(self, cascade_path):
        ''' Generates node and link settings based on cascade spreadsheet. '''
        self.resetCascade()
        loadCascadeSettingsFunc(cascade_path, settings = self)
        
    def plotCascade(self):
        ''' Plots cascade network. '''
        plotCascadeFunc(settings = self)
        


    def initCustomDatabookFramework(self):
        """
        Settings for country data book
        """
        
        self.countrybook = odict()
        self.countrybook['sheet_names'] = odict()
        self.countrybook['sheet_names']['populations'] = 'Populations'
        self.countrybook['sheet_names']['population_sizes'] = 'Population size'
        self.countrybook['sheet_names']['total_cases'] = 'Total cases'
        self.countrybook['sheet_names']['incident_cases'] = 'Incident cases'
        self.countrybook['sheet_names']['other_epidemiology'] = 'Other epidemiology'
        self.countrybook['sheet_names']['comorbidity'] = 'Comorbidity'
        self.countrybook['sheet_names']['testing_treatment'] = 'Testing and Treatment'
        self.countrybook['sheet_names']['programs'] = 'Programs'
        self.countrybook['sheet_names']['cost_coverage'] = 'Cost and coverage'
        #self.countrybook['sheet_names']['unitcost'] = 'Unit costs'
        #self.countrybook['sheet_names']['poptransitions'] = 'Population transitions'
        
        # headers for special sheets (i.e. headers aren't years)
        self.countrybook['headers'] = odict()
        self.countrybook['headers']['populations'] = ['Name','Minimum Age','Maximum Age']
        self.countrybook['headers']['programs'] = ['Name','Short name','Intervention class','Coverage indicator (annual)','Duration of treatment (days per person on average)','Frequency of intervention (in years)']
        
        # labels for each sheet
        self.countrybook['labels'] = {'populations': '',
                                      'population_sizes':'Population size: please enter population values for each year',
                                      'total_cases'     : 'TB total cases: please enter number of TB cases per year for each population and strain',
                                      'incident_cases'  : 'TB Incident cases: please enter number of new cases per year for each population group and TB strain',
                                      'other_epidemiology':'Other epidemiological data:',
                                      'comorbidity' : 'Comorbidity data:',
                                      'testing_treatment': 'Testing and treatment data: ',
                                      'programs'        : 'Enter program data:',
                                      'cost_coverage'   : 'Cost and coverage data:',
                                      }
        
        # info
        self.countrybook['disaggregations'] = odict() 
        # other values
        self.countrybook['disaggregations']['strains'] = ['DS-TB','MDR-TB','XDR-TB']
        self.countrybook['disaggregations']['smears'] = ['Smear-','Smear+'] # also potentially 'Extrapulmonary'
        self.countrybook['disaggregations']['populations'] = [] # determined dynamically at runtime
        self.countrybook['disaggregations']['regimens'] = ['DS-TB Regimen','MDR-TB Regimen','XDR-TB Regimen']
        self.countrybook['disaggregations']['programs'] = [] # determined dynamically at runtime
        self.countrybook['disaggregations']['total_pop'] = ['Total population']
        
        # for univalue sheets, includes information on how data should be disaggregated 
        self.countrybook['sheet_classes'] = odict()
        self.countrybook['sheet_classes']['univalue'] = odict()
        self.countrybook['sheet_classes']['univalue']['population_sizes'] = ['populations']
        self.countrybook['sheet_classes']['univalue']['total_cases'] = ['populations','smears','strains']
        self.countrybook['sheet_classes']['univalue']['incident_cases'] = ['populations','smears','strains']
        
        # sheet specific values
        self.countrybook['sheet_values'] = odict()
        self.countrybook['sheet_values']['other_epidemiology'] = odict()
        self.countrybook['sheet_values']['other_epidemiology']['Percentage of people vaccinated per year'] = ['populations']
        self.countrybook['sheet_values']['other_epidemiology']['Percentage of people who die from non-TB-related causes per year'] = ['populations']
        self.countrybook['sheet_values']['other_epidemiology']['Percentage of people who die from TB-related deaths per year'] = ['populations','smears','strains']
        self.countrybook['sheet_values']['other_epidemiology']['Birth rate (births per woman per year)'] = ['total_pop']
        
        self.countrybook['sheet_values']['comorbidity'] = odict()
        self.countrybook['sheet_values']['comorbidity']['HIV prevalence'] = ['populations','smears','strains']
        self.countrybook['sheet_values']['comorbidity']['Diabetes prevalence'] = ['populations','smears','strains']
        
        self.countrybook['sheet_values']['testing_treatment'] = odict()
        self.countrybook['sheet_values']['testing_treatment']['Background testing rates'] = ['populations']
        self.countrybook['sheet_values']['testing_treatment']['Testing for latent TB'] = odict()
        self.countrybook['sheet_values']['testing_treatment']['Testing for latent TB']['Percentage of population tested for latent TB per year'] = ['populations']
        self.countrybook['sheet_values']['testing_treatment']['Testing for latent TB']['Number of people initiating treatment for latent TB each year'] = ['populations']
        self.countrybook['sheet_values']['testing_treatment']['Testing for latent TB']['Number of people lost to follow up for latent TB each year'] = ['populations']
        self.countrybook['sheet_values']['testing_treatment']['Testing for latent TB']['Number of people successfully completing treatment for latent TB each year'] = ['populations']

        self.countrybook['sheet_values']['testing_treatment']['Testing for active TB'] = odict()
        self.countrybook['sheet_values']['testing_treatment']['Testing for active TB']['Percentage of population tested for active TB per year'] = ['populations']
        self.countrybook['sheet_values']['testing_treatment']['Testing for active TB']['Number of people initiating treatment for active TB each year'] = ['regimens','populations']
        self.countrybook['sheet_values']['testing_treatment']['Testing for active TB']['Number of people lost to follow up for active TB each year'] = ['regimens','populations']
        self.countrybook['sheet_values']['testing_treatment']['Testing for active TB']['Number of people successfully completing treatment for active TB each year'] = ['regimens','populations']
             
        
        self.countrybook['constants'] = {'spacing_interpopulation':2,
                                         'spacing_intrapopulation':1,
                                         'spacing_interproperty'  :4,
                                         'spacing_multivalue_label':2,
                                         'total_strains': 'Total',#All strains',
                                         'total_smears' : 'Total',
                                         'num_default_programs':28,
                                         'row_index_start':2, # for when there are no disaggregations, etc.
                                         'col_index_start':1} # 

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
    
    
    def __init__(self,level='warn'):
        
        self.defaultSettings()
        try:
            validationSettings = getattr(self, '%sSettings'%level)
            validationSettings()
            logging.info("Validation settings: %s"%level)
        
        except:
            logging.info("Could not load validation settings for level: %s"%level)
        
        
    def getValidationTypes(self):
        return ['negative_population', # runs @ validation/checkNegativePopulation
                'transition_fraction', # runs @ validation/checkTransitionFraction
                'databook_validation', # runs @ databook/loadSpreadsheetFunc
                ]        
            
    def defaultSettings(self):
        
        self.settings = odict()
        self.warnSettings()
     
    def setValidationLevel(self,validation_level):   
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
    
    
    
    def __init__(self,output='dev'):
        
        logging.info("Loading plotting settings")
        
        
        self.plotdict = {} # holder
        self.defaultSettings()
        try:
            outputSettings = getattr(self, '%sSettings'%output)
            outputSettings()
        except:
            logging.info("Could not load rcParams for plotting for output: %s"%output)
        
    def KMSuffixFormatter(self,x,pos):
            'The two args are the value and tick position'
            if x >= 1e7:
                return '%1.1fM' % (x*1e-6)
            elif x >= 1e6:
                return '%1.fM' % (x*1e-6)
            elif x >= 1e3:
                return '%1.fK' % (x*1e-3) 
            else:
                return x
            
    def defaultSettings(self):
        
        pl.rcParams['font.size'] = 12
        pl.rcParams['font.family'] = 'sans-serif'
        
        pl.rcParams['savefig.dpi'] = 300
        pl.rcParams['savefig.format'] = 'png'
        pl.rcParams['savefig.transparent'] =  'True'
        
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
        pl.rcParams['axes.titlesize'] = 1.5*pl.rcParams['font.size']
    
        pl.rcParams['lines.linewidth'] = 3
        pl.rcParams['lines.marker'] = 'None'
        
        # Non-standard list of parameters used in plotting
        self.plotdict = {# scatter plotting values
                         'marker' : 'o',
                         'facecolors' : 'none',
                         's' : 40,
                         # axes format
                         'year_inc':5,
                         # colormapping for category lists
                         'colormapping_order':'alternate3',# as we have triplets in undiagnosed --> diagnosed --> on treatment
                         'formatter': FuncFormatter(self.KMSuffixFormatter) } 


    def devSettings(self):
        pl.rcParams['figure.figsize'] = (10, 8)
        pl.rcParams['savefig.dpi'] = 100
        pl.rcParams['savefig.transparent'] =  'False' # relax
        
    def printSettings(self):
        
        pl.rcParams['figure.figsize'] = (15, 10)
        pl.rcParams['savefig.dpi'] = 300
        pl.rcParams['savefig.transparent'] =  'True' # enforce
        self.plotdict['legend_off'] = True
        
    def presentationSettings(self):
        pl.rcParams['font.size'] = 16
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
        
        pl.rcParams['savefig.transparent'] =  'True' # enforce
        
        self.plotdict['legend_off'] = True
        
        
        
