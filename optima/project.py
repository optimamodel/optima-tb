#%% Imports

from utils import tic, toc, odict, OptimaException
from model import runModel
from settings import Settings
from parameters import ParameterSet
from plotting import Plotter
from databook import makeSpreadsheetFunc, loadSpreadsheetFunc

from uuid import uuid4 as uuid



#%% Project class (i.e. one self-contained geographical unit)

class Project(object):
    ''' The main Optima project class. Almost all Optima functionality is provided by this class. '''

    def __init__(self, name = 'default', cascade_path = '../data/cascade.xlsx'):
        ''' Initialize project. '''

        self.name = name
        self.uid = uuid()
        
        self.settings = Settings(cascade_path = cascade_path)        
        self.data = odict()

        self.parsets = odict()
        self.results = odict()
        
        self.plotter = Plotter({})
        
        
    def runSim(self, parset_name = 'default', plot = False):
        ''' Run model using a selected parset and store/return results. '''
        
        if len(self.parsets) < 1: raise OptimaException('ERROR: Project "%s" appears to have no parameter sets. Cannot run model.' % self.name)
        try: parset = self.parsets[parset_name]
        except: raise OptimaException('ERROR: Project "%s" is lacking a parset named "%s". Cannot run model.' % (self.name, parset_name))

        tm = tic()
        results, sim_settings, outputs = runModel(settings = self.settings, parset = parset)
        toc(tm, label = 'running %s model' % self.name)
        
        if plot:
            tp = tic()
            self.plotter.updateData(self.data)
            self.plotter.plotProjectResults(results,outputs,sim_settings,self.settings.charac_specs,title=self.name.title())
            toc(tp, label = 'plotting %s' % self.name)
        
        return results, outputs
        
    
    def makeSpreadsheet(self, num_pops = 5, num_migrations = 2):
        ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''
        
        databook_path = '../data/' + self.name + '-data.xlsx'
        makeSpreadsheetFunc(settings = self.settings, databook_path = databook_path, num_pops = num_pops, num_migrations = num_migrations)        
        
    
    def loadSpreadsheet(self, databook_path = None):
        ''' Load data spreadsheet into Project data dictionary. '''
        
        if databook_path is None: databook_path = '../data/' + self.name + '-data.xlsx'
        self.data = loadSpreadsheetFunc(settings = self.settings, databook_path = databook_path) 
        
    


    def makeParset(self, name = 'default'):
        ''' Transform project data into a set of parameters that can be used in model simulations. '''

        if not self.data: raise OptimaException('ERROR: No data exists for project "%s".' % self.name)
        self.parsets[name] = ParameterSet(name = name)
        self.parsets[name].makePars(self.data)