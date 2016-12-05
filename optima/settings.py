#%% Imports
from utils import odict
from parsing import FunctionParser
from cascade import loadCascadeSettingsFunc, plotCascadeFunc


import logging
logger = logging.getLogger(__name__)

import pylab as pl
import numpy as np



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

        self.tvec_start = 2000.0     # Default start year for data input and simulations.
        self.tvec_end = 2030.0       # Default end year for data input and simulations.
        self.tvec_observed_end = 2015.0
        self.tvec_dt = 1.0/4         # Default timestep for simulations.
        
        self.recursion_limit = 100      # Limit for recursive references, primarily used in avoiding circular references for definitions using dependencies.
        self.fit_metric = 'wape'
        
        self.parser = FunctionParser(debug = False)      # Decomposes and evaluates functions written as strings, in accordance with a grammar defined within the parser object.
                
#        # Initialize all cascade and databook parameters for fresh import, via reset
#        self.resetCascade()
        
        # Settings for databooks / spreadsheets / workbooks / cascades:
        self.loadCascadeSettings(cascade_path)
        self.initCustomDatabookFramework()              # creates self.countrybook
                                                        # NOTE: databook will hopefully one day be capable of replicating countrybook, making the two different
                                                        # unnecessary.
                                                         
        logging.info("Created settings based on cascade: %s"%cascade_path)
        
        self.plot_settings = PlottingSettings()        
        
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
        
        self.links = odict()                    # Key is a tag. Value is a list of compartment-label tuple.
        self.linkpar_specs = odict()            # Key is a link-parameter label. Value is a dict including link tag, link-parameter name, default value.
        self.linkpar_name_labels = odict()      # Key is a link-parameter name. Value is a link-parameter label. (A partial reversed linkpar-specs.)
        
        self.par_funcs = odict()                # A definition-ordered dictionary of parameters that are functionally defined in terms of other parameters/characteristics.
        self.par_deps = odict()                 # A definition-ordered dictionary of 'untagged' parameters used as dependencies for other parameters.
        self.charac_deps = {}                   # An unordered dictionary of characteristics that must be calculated at each model timestep due to being dependencies for other variables.
                                                # Should correspond to every item in charac_specs that has a 'par_dependency' tag.
        
        # Project-data workbook metadata.
        self.databook = odict()
        self.databook['sheet_names'] = odict()
        self.databook['sheet_names']['pops'] =      'Population Definitions'
        self.databook['sheet_names']['transmat'] =  'Transfer Definitions'
        self.databook['sheet_names']['transval'] =  'Transfer Details'
        self.databook['sheet_names']['charac'] =    'Epidemic Characteristics'
        self.databook['sheet_names']['linkpars'] =  'Cascade Parameters'
        self.databook['custom_sheet_names'] = odict()
        self.make_sheet_characs = True      # Tag for whether the default characteristics worksheet should be generated during databook creation.
        self.make_sheet_linkpars = True     # Tag for whether the default cascade parameters worksheet should be generated during databook creation.

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
    
    
    
    def __init__(self,level='warn'):
        
        self.defaultSettings()
        try:
            validationSettings = getattr(self, '%sSettings'%level)
            validationSettings()
            logging.info("Validation settings: %s"%level)
        
        except:
            logging.info("Could not load validation settings for level: %s"%level)
        
        
    def getValidationTypes(self):
        return ['negative_population']        
            
    def defaultSettings(self):
        
        self.settings = odict()
        self.errorSettings()
     
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
        self.defaultSettings()
        try:
            outputSettings = getattr(self, '%sSettings'%output)
            outputSettings()
        except:
            logging.info("Could not load rcParams for plotting for output: %s"%output)
        
        
    def defaultSettings(self):
        
        pl.rcParams['font.size'] = 12
        pl.rcParams['font.family'] = 'sans-serif'
        
        pl.rcParams['savefig.dpi'] = 300
        
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
        pl.rcParams['lines.markersize'] = 40
        pl.rcParams['lines.markeredgewidth'] = 3
    


    def devSettings(self):
        pl.rcParams['figure.figsize'] = (10, 8)
        pl.rcParams['savefig.dpi'] = 300
        
    def printSettings(self):
        
        pl.rcParams['figure.figsize'] = (15, 10)