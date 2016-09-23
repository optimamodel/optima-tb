#%% Imports

from utils import tic, toc, odict, OptimaException
from model import model
from settings import Settings
from parameters import ParameterSet
from plotting import gridColorMap
from databook import makeSpreadsheetFunc, loadSpreadsheetFunc

import pylab as pl
from uuid import uuid4 as uuid
from copy import deepcopy as dcp



#%% Project class (i.e. one self-contained geographical unit)

class Project(object):
    ''' The main Optima project class. Almost all Optima functionality is provided by this class. '''

    def __init__(self, name = 'default', cascade_path = './cascade.xlsx'):
        ''' Initialize project. '''

        self.name = name
        self.uid = uuid()
        
        self.settings = Settings(cascade_path = cascade_path)        
        self.data = odict()

        self.parsets = odict()
        self.results = odict()
        
        
    def runSim(self, parset_name = 'default'):
        ''' Run model using a selected parset and store/return results. '''
        
        if len(self.parsets) < 1: raise OptimaException('ERROR: Project %s appears to have no parameter sets. Cannot run model.' % self.name)
        try: parset = self.parsets[parset_name]
        except: raise OptimaException('ERROR: Project %s is lacking a parset named %s. Cannot run model.' % (self.name, parset_name))

        tm = tic()
        results, sim_settings = model(settings = self.settings, parset = parset)
        toc(tm, label = 'running %s model' % self.name)
        
        tp = tic()
        for pop_oid in results:
            pop = results[pop_oid]
            
            fig, ax = pl.subplots(figsize=(15,10))
            colors = gridColorMap(len(pop.nodes))
            bottom = 0*sim_settings['tvec']
            
            k = 0
            for node in pop.nodes:
                top = bottom + node.popsize
                
                ax.fill_between(sim_settings['tvec'], bottom, top, facecolor=colors[k], alpha=1, lw=0)
                ax.plot((0, 0), (0, 0), color=colors[k], linewidth=10)
                bottom = dcp(top)
                k += 1
            
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])   
            
            legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5), 'ncol':2}
            ax.set_title('%s Cascade - %s' % (self.name.title(), pop.name.title()))
            ax.set_xlabel('Year')
            ax.set_ylabel('People')
            ax.set_xlim((sim_settings['tvec'][0], sim_settings['tvec'][-1]))
            ax.set_ylim((0, max(top)))
            cascadenames = [node.name for node in pop.nodes]
            ax.legend(cascadenames, **legendsettings)
        toc(tp, label = 'plotting %s' % self.name)
        
        return results
        
    
    def makeSpreadsheet(self, num_pops = 5):
        ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''
        
        databook_path = './' + self.name + '-data.xlsx'
        makeSpreadsheetFunc(settings = self.settings, databook_path = databook_path, num_pops = num_pops)        
        
    
    def loadSpreadsheet(self, databook_path = None):
        ''' Load data spreadsheet into Project data dictionary. '''
        
        if databook_path is None: databook_path = './' + self.name + '-data.xlsx'
        self.data = loadSpreadsheetFunc(settings = self.settings, databook_path = databook_path) 


    def makeParset(self, name = 'default'):
        ''' Transform project data into a set of parameters that can be used in model simulations. '''

        if not self.data: raise OptimaException('ERROR: No data exists for project %s.' % self.name)
        self.parsets[name] = ParameterSet(name = name)
        self.parsets[name].makePars(self.data)