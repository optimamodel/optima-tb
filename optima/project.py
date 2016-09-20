#%% Imports

from utils import tic, toc, odict, OptimaException
from model import model
from settings import Settings
from parameters import ParameterSet
from plotting import gridColorMap

import pylab as pl
from uuid import uuid4 as uuid
from copy import deepcopy as dcp

import xlrd
import xlsxwriter as xw
import numpy as np
from xlsxwriter.utility import xl_rowcol_to_cell as rc



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
        
        
    def makeSpreadsheet(self, num_pops = 5):
        ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''
        
        databook_path = './' + self.name + '-data.xlsx'
        workbook = xw.Workbook(databook_path)
        ws_pops = workbook.add_worksheet(self.settings.databook['sheet_names']['pops'])
        ws_linkpars = workbook.add_worksheet(self.settings.databook['sheet_names']['linkpars'])
        
        data_tvec = np.arange(self.settings.tvec_start, self.settings.tvec_end + 1.0/2)
        offset_tvec = 3     # Offset to denote at which column the time vector begins in spreadsheet.
        
        # Population names sheet.
        ws_pops_width = 15
        ws_pops.write(0, 0, 'Name')
        ws_pops.write(0, 1, 'Abbreviation')
        temp_pop_names = []
        for k in xrange(num_pops):
            temp_pop_name = 'Population '+str(k+1)
            temp_pop_names.append(temp_pop_name)
            ws_pops.write(k+1, 0, temp_pop_name)
            ws_pops.write(k+1, 1, '=LEFT(%s,3)&"%i"' % (rc(k+1,0), k+1))
        ws_pops.set_column(0, 1, ws_pops_width)
        
        # Transition parameters sheet.
        ws_linkpars_width = 40
        row_id = 0
        for link_name in self.settings.linkpar_name_labels.keys():
            ws_linkpars.write(row_id, 0, link_name)
            ws_linkpars.write(row_id, 1, 'Constant')
            for k in xrange(len(data_tvec)):
                ws_linkpars.write(row_id, k+offset_tvec, data_tvec[k])
            for pid in xrange(num_pops):
                row_id += 1
                ws_linkpars.write(row_id, 0, "='Population Definitions'!%s" % rc(pid+1,0), None, temp_pop_names[pid])
                ws_linkpars.write(row_id, 1, '=IF(SUMPRODUCT(--(%s:%s<>""))=0,0.0,"N.A.")' % (rc(row_id,offset_tvec), rc(row_id,offset_tvec+len(data_tvec)-1)))
                ws_linkpars.write(row_id, 2, 'OR')
                
                # Values for all extra populations default to first population values.
                if pid > 0:
                    for k in xrange(len(data_tvec)):
                        ws_linkpars.write(row_id, k+offset_tvec, '=IF(%s="","",%s)' % (rc(row_id-pid,k+offset_tvec),rc(row_id-pid,k+offset_tvec)))
            
            row_id += 2
        ws_linkpars.set_column(0, 0, ws_linkpars_width)
        
        workbook.close()
        
    
    # NOTE: This needs so much quality-assurance testing. Need to ensure that input data sheet aligns with cascade settings.
    #       It would also be good to align metadata between making and loading spreadsheets.
    def loadSpreadsheet(self, databook_path = None):
        ''' Load data spreadsheet into Project data dictionary. '''
        
        if databook_path is None: databook_path = './' + self.name + '-data.xlsx'
        try: workbook = xlrd.open_workbook(databook_path)
        except: raise OptimaException('ERROR: Project data workbook was unable to be loaded from... %s' % databook_path)
        ws_pops = workbook.sheet_by_name(self.settings.databook['sheet_names']['pops'])
        ws_linkpars = workbook.sheet_by_name(self.settings.databook['sheet_names']['linkpars'])

        self.data = odict()
        self.data['pops'] = odict()
        self.data['pops']['name_labels'] = odict()
        self.data['pops']['label_names'] = odict()      # A reverse of the previous dictionary.
        for row_id in xrange(ws_pops.nrows):
            if row_id > 0 and ws_pops.cell_value(row_id, 0) not in ['']:
                pop_label = 'pop' + str(row_id)
                self.data['pops']['name_labels'][str(ws_pops.cell_value(row_id, 0))] = pop_label
                self.data['pops']['label_names'][pop_label] = str(ws_pops.cell_value(row_id, 0))
                
        self.data['linkpars'] = odict()
        current_linkpar_name = None
        current_linkpar_label = None
        for row_id in xrange(ws_linkpars.nrows):
            val = str(ws_linkpars.cell_value(row_id, 0))
            if val in ['']:
                current_linkpar_name = None
            elif current_linkpar_name is None:
                current_linkpar_name = val
                current_linkpar_label = self.settings.linkpar_name_labels[val]
                self.data['linkpars'][current_linkpar_label] = odict()
            else:
                current_pop_label = self.data['pops']['name_labels'][val]
                self.data['linkpars'][current_linkpar_label][current_pop_label] = odict()
                self.data['linkpars'][current_linkpar_label][current_pop_label]['t'] = None
                self.data['linkpars'][current_linkpar_label][current_pop_label]['y'] = ws_linkpars.cell_value(row_id, 1)


    def makeParset(self, name = 'default'):
        ''' Transform project data into a set of parameters that can be used in model simulations. '''

        if not self.data: raise OptimaException('ERROR: No data exists for project %s.' % self.name)
        self.parsets[name] = ParameterSet(name = name)
        self.parsets[name].makePars(self.data)