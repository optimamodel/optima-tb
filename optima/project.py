#%% Imports

from utils import tic, toc, odict
from model import model
from settings import Settings
from plotting import gridColorMap

import pylab as pl
from uuid import uuid4 as uuid
from copy import deepcopy as dcp

import xlsxwriter as xw
import numpy as np



#%% Project class (i.e. one self-contained geographical unit)

class Project(object):
    ''' The main Optima project class. Almost all Optima functionality is provided by this class. '''

    def __init__(self, name = 'default', cascade_path = './cascade.xlsx'):
        ''' Initialize project. '''

        self.name = name
        self.uid = uuid()
        
        self.settings = Settings(cascade_path = cascade_path)        
        self.data = {}

        self.parsets = odict()
        self.results = odict()

        return None
        
        
    def runSim(self):
        ''' Run model using a selected parset and store/return results. '''

        tm = tic()
        results, sim_settings = model(settings = self.settings)
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
            
            legendsettings = {'loc':'center left', 'bbox_to_anchor':(1.05, 0.5)}
            ax.set_title('%s Cascade - %s' % (self.name.title(), pop.name.title()))
            ax.set_xlabel('Year')
            ax.set_ylabel('People')
            ax.set_xlim((sim_settings['tvec'][0], sim_settings['tvec'][-1]))
            ax.set_ylim((0, max(top)))
            cascadenames = [node.name for node in pop.nodes]
            ax.legend(cascadenames, **legendsettings)
        toc(tp, label = 'plotting %s' % self.name)
        
        return results
        
        
    def makeSpreadsheet(self):
        ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''
        
        databook_path = './' + self.name + '-data.xlsx'
        workbook = xw.Workbook(databook_path)
        ws_linkpars = workbook.add_worksheet('Transition Parameters')
        
        data_tvec = np.arange(self.settings.tvec_start, self.settings.tvec_end + 1.0/2)
        
        row_id = 0
        for linkpar_label in self.settings.linkpar_specs.keys():
            ws_linkpars.write(row_id, 0, self.settings.linkpar_specs[linkpar_label]['name'])
            for k in xrange(len(data_tvec)):
                ws_linkpars.write(row_id, k+1, data_tvec[k])
            row_id += 1
            
            row_id += 2
        
        workbook.close()
        
        
        
#%% Test project functionality
tt = tic()
ts1 = tic()
p1 = Project(name = 'simple-test', cascade_path = './cascade-simple.xlsx')
toc(ts1, label = 'creating %s project' % p1.name)
r1 = p1.runSim()
p1.settings.plotCascade()
p1.makeSpreadsheet()
ts2 = tic()
p2 = Project(name = 'standard-test')
toc(ts2, label = 'creating %s project' % p2.name)
p2.runSim()
p2.settings.plotCascade()
p2.makeSpreadsheet()
toc(tt, label = 'entire process')