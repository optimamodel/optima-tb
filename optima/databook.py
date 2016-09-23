#%% Imports

from utils import odict, OptimaException

import xlrd
import xlsxwriter as xw
import numpy as np
from xlsxwriter.utility import xl_rowcol_to_cell as rc
from numbers import Number



#%% Constants used across spreadsheet functions

default_path = './project-data.xlsx'    # Default path if not called via a project method.
offset_tvec = 3     # Offset to denote at which column the time vector begins in spreadsheet.

# NOTE: Should cascade metadata like spreadsheet titles be moved here too...?


#%% Function for generating project databook

# NOTE: Comment this better later, especially with the fact that all Excel formulae need a fall-back value.
def makeSpreadsheetFunc(settings, databook_path = default_path, num_pops = 5):
    ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''
    
    workbook = xw.Workbook(databook_path)
    ws_pops = workbook.add_worksheet(settings.databook['sheet_names']['pops'])
    ws_linkpars = workbook.add_worksheet(settings.databook['sheet_names']['linkpars'])
    
    data_tvec = np.arange(settings.tvec_start, settings.tvec_end + 1.0/2)
    ws_pops_width = 15
    ws_linkpars_width = 40
    assumption_width = 10
    
    # Population names sheet.
    ws_pops.write(0, 0, 'Name')
    ws_pops.write(0, 1, 'Abbreviation')
    temp_pop_names = []
    for k in xrange(num_pops):
        temp_pop_name = 'Population '+str(k+1)
        temp_pop_names.append(temp_pop_name)
        ws_pops.write(k+1, 0, temp_pop_name)
        ws_pops.write(k+1, 1, '=LEFT(%s,3)&"%i"' % (rc(k+1,0), k+1), None, temp_pop_name[0:3]+str(k+1))
    ws_pops.set_column(0, 1, ws_pops_width)
    
    # Transition parameters sheet.
    row_id = 0
    for link_name in settings.linkpar_name_labels.keys():
        link_label = settings.linkpar_name_labels[link_name]
        def_val = 0
        if 'default' in settings.linkpar_specs[link_label]:
            def_val = settings.linkpar_specs[link_label]['default']
            
        ws_linkpars.write(row_id, 0, link_name)
        ws_linkpars.write(row_id, 1, 'Assumption')
        for k in xrange(len(data_tvec)):
            ws_linkpars.write(row_id, k+offset_tvec, data_tvec[k])
        for pid in xrange(num_pops):
            row_id += 1
            ws_linkpars.write(row_id, 0, "='Population Definitions'!%s" % rc(pid+1,0), None, temp_pop_names[pid])
            ws_linkpars.write(row_id, 1, '=IF(SUMPRODUCT(--(%s:%s<>""))=0,%f,"N.A.")' % (rc(row_id,offset_tvec), rc(row_id,offset_tvec+len(data_tvec)-1), def_val), None, def_val)
            ws_linkpars.write(row_id, 2, 'OR')
            
            # Values for all extra populations default to first population values.
            if pid > 0:
                for k in xrange(len(data_tvec)):
                    ws_linkpars.write(row_id, k+offset_tvec, '=IF(%s="","",%s)' % (rc(row_id-pid,k+offset_tvec),rc(row_id-pid,k+offset_tvec)), None, '')
        
        row_id += 2
    ws_linkpars.set_column(0, 0, ws_linkpars_width)
    ws_linkpars.set_column(1, 1, assumption_width)
    
    workbook.close()
    

#%% Function for loading project databook

# NOTE: This needs so much quality-assurance testing. Need to ensure that input data sheet aligns with cascade settings.
def loadSpreadsheetFunc(settings, databook_path = None):
    ''' Load data spreadsheet into Project data dictionary. '''
    
    if databook_path is None: databook_path = default_path
    try: workbook = xlrd.open_workbook(databook_path)
    except: raise OptimaException('ERROR: Project data workbook was unable to be loaded from... %s' % databook_path)
    ws_pops = workbook.sheet_by_name(settings.databook['sheet_names']['pops'])
    ws_linkpars = workbook.sheet_by_name(settings.databook['sheet_names']['linkpars'])

    # Population names sheet.
    data = odict()
    data['pops'] = odict()
    data['pops']['name_labels'] = odict()
    data['pops']['label_names'] = odict()      # A reverse of the previous dictionary.
    for row_id in xrange(ws_pops.nrows):
        if row_id > 0 and ws_pops.cell_value(row_id, 0) not in ['']:
            pop_label = 'pop' + str(row_id)
            data['pops']['name_labels'][str(ws_pops.cell_value(row_id, 0))] = pop_label
            data['pops']['label_names'][pop_label] = str(ws_pops.cell_value(row_id, 0))
    
    # Transition parameters sheet.
    data['linkpars'] = odict()
    current_linkpar_name = None
    current_linkpar_label = None
    pop_id = 0
    for row_id in xrange(ws_linkpars.nrows):
        val = str(ws_linkpars.cell_value(row_id, 0))
        if val in ['']:
            current_linkpar_name = None
            pop_id = 0
        elif current_linkpar_name is None:
            current_linkpar_name = val
            current_linkpar_label = settings.linkpar_name_labels[val]
            data['linkpars'][current_linkpar_label] = odict()
        else:
            current_pop_label = data['pops']['name_labels'][val]
            if current_pop_label != data['pops']['label_names'].keys()[pop_id]:
                raise OptimaException('ERROR: Somewhere in the transition parameters sheet, populations are not ordered as in the population definitions sheet.')
            data['linkpars'][current_linkpar_label][current_pop_label] = odict()
            
            # Run through the rows beneath the year range, but only if there is not a number in the cell corresponding to assumption.
            list_t = []
            list_y = []
            for col_id in xrange(ws_linkpars.ncols):
                if col_id > 0 and isinstance(ws_linkpars.cell_value(row_id, col_id), Number):
                    list_y.append(float(ws_linkpars.cell_value(row_id, col_id)))
                    if not isinstance(ws_linkpars.cell_value(row_id-1-pop_id, col_id), Number):
                        list_t.append(float(ws_linkpars.cell_value(row_id-1-pop_id, offset_tvec)))
                        break
                    else:
                        list_t.append(float(ws_linkpars.cell_value(row_id-1-pop_id, col_id)))
            data['linkpars'][current_linkpar_label][current_pop_label]['t'] = np.array(list_t)
            data['linkpars'][current_linkpar_label][current_pop_label]['y'] = np.array(list_y)                
            
            pop_id += 1
            
    return data