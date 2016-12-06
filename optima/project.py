#%% Imports
import logging
import logging.config

logging.config.fileConfig('./logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

from utils import tic, toc, odict, OptimaException
from model import runModel
from settings import Settings
from parameters import ParameterSet
from plotting import Plotter
from databook import makeSpreadsheetFunc, loadSpreadsheetFunc
from calibration import makeManualCalibration, calculateFitFunc, performAutofit

from uuid import uuid4 as uuid
from numpy import max


#%% Project class (i.e. one self-contained geographical unit)

class Project(object):
    ''' The main Optima project class. Almost all Optima functionality is provided by this class. '''

    def __init__(self, name = 'default', cascade_path = './data/cascade.xlsx', **args):
        ''' Initialize project. '''

        self.name = name
        self.uid = uuid()
        
        self.settings = Settings(cascade_path = cascade_path, **args)        
        self.data = odict()

        self.parsets = odict()
        self.results = odict()
        
        self.plotter = Plotter({})
        
        logger.info("Created project: %s"%self.name)
        
    def setYear(self, yearRange,observed_data=True):
        '''
        
        @param yearRange: tuple or list 
        @param observed_data: bool indicating whether to change observed date end or simulation date end 
        '''
        self.settings.tvec_start = yearRange[0]
        if observed_data:
            self.settings.tvec_observed_end = yearRange[1]
        else:
            self.settings.tvec_end = yearRange[1]
    
    
    def runSim(self, parset_name = 'default', parameterset = None, plot = False):
        ''' Run model using a selected parset and store/return results. '''
        
        if parameterset is not None:
            parset = parameterset
        elif len(self.parsets) < 1: 
            raise OptimaException('ERROR: Project "%s" appears to have no parameter sets. Cannot run model.' % self.name)
        else:
            try: parset = self.parsets[parset_name]
            except: raise OptimaException('ERROR: Project "%s" is lacking a parset named "%s". Cannot run model.' % (self.name, parset_name))

        tm = tic()
        #results, sim_settings, outputs = runModel(settings = self.settings, parset = parset)
        robj, results, sim_settings, outputs = runModel(settings = self.settings, parset = parset)
        toc(tm, label = 'running %s model' % self.name)
        
        if plot:
            tp = tic()
            self.plotter.updateData(self.data)
            self.plotter.plotProjectResults(results, outputs, sim_settings, self.settings.charac_specs, title = self.name.title())
            toc(tp, label = 'plotting %s' % self.name)
        
        return results, outputs, sim_settings, robj
        
    
    def makeSpreadsheet(self, databook_path = None, num_pops = 5, num_migrations = 2):
        ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''
        
        if databook_path is None: databook_path = '../data/' + self.name + '-data.xlsx'
        logging.info("Attempting to create databook %s"%databook_path)
        
        makeSpreadsheetFunc(settings = self.settings, databook_path = databook_path, num_pops = num_pops, num_migrations = num_migrations)        
        
    
    def loadSpreadsheet(self, databook_path = None):
        ''' Load data spreadsheet into Project data dictionary. '''
        if databook_path is None: databook_path = '../data/' + self.name + '-data.xlsx'
        logging.info("Attempting to load databook %s"%databook_path)
        
        self.data = loadSpreadsheetFunc(settings = self.settings, databook_path = databook_path) 
        


    def makeParset(self, name = 'default'):
        ''' Transform project data into a set of parameters that can be used in model simulations. '''

        if not self.data: raise OptimaException('ERROR: No data exists for project "%s".' % self.name)
        self.parsets[name] = ParameterSet(name = name)
        self.parsets[name].makePars(self.data)
        
        
    def makeManualCalibration(self, parset_name, rate_dict):
        ''' Take in dictionary of updated values for rates that can be used for manual calibration
            Update values dict: {pop_name : {parameter : value} '''
        if not parset_name in self.parsets.keys():
            self.makeParset(name=parset_name)
        paramset = self.parsets[parset_name]
        logging.info("Updating parameter values in parset=%s"%(parset_name))
        
        makeManualCalibration(paramset,rate_dict)
    
    
    def calculateFit(self,results,metric=None):
        '''
        Calculates the score for the fit during manual calibration and prints to output. 
        
        Params:
            results    resultSet object
            metric     type of metric used (defaults to default specified in settings).
        
        Future: can consider saving results into resultset and accessing results by name / parset name
        '''
        if metric is None:
            metric = self.settings.fit_metric
        
        if results is None:
            raise OptimaException('ERROR: no result is specified. Cannot calculate fit.')
        
        if self.data is None or len(self.data)==0:
            raise OptimaException('ERROR: no data is specified. Cannot calculate fit.')
        
        datapoints = results.getCharacteristicDatapoints()
        score = calculateFitFunc(datapoints,results.t_observed_data,self.data['characs'],metric)
        logger.info("Calculated scores for fit using %s: largest value=%.2f"%(metric,max(score)))
      
      
    def runAutofitCalibration(self,new_parset_name = None, old_parset_name="default"):
        """
        
        TODO: implement changes with old and new parameter set name
        """
        
        if not old_parset_name in self.parsets.keys():
            self.makeParset(name=old_parset_name)
        paramset = self.parsets[old_parset_name]
        
        if new_parset_name is None:
            # TODO: check that autofit doesn't already exist; if so, add suffix
            new_parset_name = "autofit" 
        
        logger.info("About to run autofit on parameters using parameter set = %s"%old_parset_name)
        new_parset = performAutofit(self,paramset,new_parset_name=new_parset_name,**self.settings.autofit_params)
        logger.info("Created new parameter set '%s' using autofit"%new_parset_name)
        self.parsets[new_parset_name] = new_parset
        
         