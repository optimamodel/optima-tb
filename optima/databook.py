#%% Imports

from utils import odict, OptimaException

import xlrd
import xlsxwriter as xw
import numpy as np
from xlsxwriter.utility import xl_rowcol_to_cell as rc
from numbers import Number
from copy import deepcopy as dcp



#%% Constants used across spreadsheet functions

default_path = './project-data.xlsx'    # Default path if not called via a project method.

# NOTE: Should cascade metadata like spreadsheet titles be moved here too...?


#%% Utility functions to generate sub-blocks of the project databook

def makeValueEntryArrayBlock(worksheet, at_row, at_col, num_arrays, tvec, assumption = 0.0, data_formats = None):
    '''
    Create a block where users choose data-entry format and enter values either as an assumption or time-dependent array.
    
    Args:
        worksheet       -   The worksheet in which to produce the block.
        at_row          -   The row at which the top-left corner of the block will be fixed.
        at_col          -   The column at which the top-left corner of the block will be fixed.
        num_arrays      -   The number of input arrays to generate.
        tvec            -   The year range for which values are (optionally) to be entered.
        assumption      -   The default value to enter in the assumption column.
        data_formats    -   A list of formats that data can be entered as.
                            This is not data validation, but will determine the form of calculations during model processing.
    '''
    
    if data_formats is None: data_formats = ['Probability', 'Fraction', 'Number']
    offset = at_col + 3     # This is the column at which the input year range and corresponding values should be written.
    
    worksheet.write(at_row, at_col, 'Format')
    worksheet.write(at_row, at_col + 1, 'Assumption')
    for k in xrange(len(tvec)):
        worksheet.write(at_row, offset + k, tvec[k])  
    
    for aid in xrange(num_arrays):
        row_id = at_row + aid + 1
        offset = at_col + 3
        worksheet.write(row_id, at_col, data_formats[0])      # Default choice for data format.  
        worksheet.data_validation('%s' % rc(row_id,at_col), {'validate': 'list', 'source': data_formats, 'ignore_blank': False})       
        worksheet.write(row_id, at_col + 1, '=IF(SUMPRODUCT(--(%s:%s<>""))=0,%f,"N.A.")' % (rc(row_id,offset), rc(row_id,offset+len(tvec)-1), assumption), None, assumption)
        worksheet.write(row_id, at_col + 2, 'OR')
        
        # Make changes to the first row be mirrored in the other rows.
        if aid > 0:
            for k in xrange(len(tvec)):
                worksheet.write(row_id, offset + k, '=IF(%s="","",%s)' % (rc(row_id-aid,offset+k),rc(row_id-aid,offset+k)), None, '')
               
               
def makeConnectionMatrix(worksheet, at_row, at_col, labels, formula_labels = None, allowed_vals = None):
    '''
    Create a matrix where users tag connections from a (left-column positioned) list to the same (top-row positioned) list.
    Diagonals (self-connections) cannot be filled with values.
    
    Args:
        worksheet       -   The worksheet in which to produce the matrix.
        at_row          -   The row at which the top-left corner of the matrix will be fixed.
        at_col          -   The column at which the top-left corner of the matrix will be fixed.
        labels          -   A list of labels for objects to be connected (e.g. population groups).
        formula_labels  -   An optional list corresponding to labels that contains Excel formulas.
                            The labels list is still required for default values.
                            In the case that the databook is not opened in Excel prior to loading, the defaults are read in.
        allowed_vals    -   A list of allowed values that matrix elements can be tagged by.
    '''
    
    if allowed_vals is None: allowed_vals = ['n', 'y']
    if formula_labels is None: formula_labels = dcp(labels)
    if len(formula_labels) != len(labels): raise OptimaException('ERROR: The numbers of formula-based and default labels passed to a matrix-writing process do not match.')
    
    for k in xrange(len(labels)):
        worksheet.write(at_row + k + 1, at_col, formula_labels[k], None, labels[k])
        worksheet.write(at_row, at_col + k + 1, formula_labels[k], None, labels[k])
        
    for row_pre in xrange(len(labels)):
        row_id = at_row + row_pre + 1
        for col_pre in xrange(len(labels)):
            col_id = at_col + col_pre + 1
            
            if row_pre != col_pre:
                worksheet.write(row_id, col_id, allowed_vals[0])      # Default choice for data format.  
                worksheet.data_validation('%s' % rc(row_id,col_id), {'validate': 'list', 'source': allowed_vals, 'ignore_blank': False})
            else:
                worksheet.write(row_id, col_id, '')
                worksheet.data_validation('%s' % rc(row_id,col_id), {'validate': 'list', 'source': [''], 'ignore_blank': False})
    
    


#%% Function for generating project databook

# NOTE: Comment this better later, especially with the fact that all Excel formulae need a fall-back value.
def makeSpreadsheetFunc(settings, databook_path = default_path, num_pops = 5, num_migrations = 2):
    ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''
    
    workbook = xw.Workbook(databook_path)
    ws_pops = workbook.add_worksheet(settings.databook['sheet_names']['pops'])
    ws_poptrans = workbook.add_worksheet(settings.databook['sheet_names']['poptrans'])
    ws_linkpars = workbook.add_worksheet(settings.databook['sheet_names']['linkpars'])
    
    data_tvec = np.arange(settings.tvec_start, settings.tvec_end + 1.0/2)
    ws_pops_width = 15
    ws_poptrans_width = 15
    ws_linkpars_width = 40
    assumption_width = 10
    
    # Population names sheet.
    ws_pops.write(0, 0, 'Name')
    ws_pops.write(0, 1, 'Abbreviation')
    ws_pops.write(0, 2, 'Minimum Age')
    ws_pops.write(0, 3, 'Maximum Age')
    
    # While writing default population names and labels, they are stored for future reference as well.
    pop_names_default = []
    pop_labels_default = []
    for pid in xrange(num_pops):
        pop_name = 'Population '+str(pid+1)
        pop_label = pop_name[0:3]+str(pid+1)
        pop_names_default.append(pop_name)
        pop_labels_default.append(pop_label)
        ws_pops.write(pid+1, 0, pop_name)
        ws_pops.write(pid+1, 1, '=LEFT(%s,3)&"%i"' % (rc(pid+1,0), pid+1), None, pop_label)
    ws_pops.set_column(0, 4, ws_pops_width)
    
    # Excel formulae strings that point to population names and labels are likewise stored.
    pop_names_formula = []
    pop_labels_formula = []
    for pid in xrange(num_pops):
        pop_names_formula.append("='%s'!%s" % (settings.databook['sheet_names']['pops'], rc(pid+1,0)))
        pop_labels_formula.append("='%s'!%s" % (settings.databook['sheet_names']['pops'], rc(pid+1,1)))
    
    # Inter-population transitions sheet.
    # Produce aging matrix.
    ws_poptrans.write(0, 0, 'Aging')
    makeConnectionMatrix(worksheet = ws_poptrans, at_row = 0, at_col = 0, labels = pop_labels_default, formula_labels = pop_labels_formula)
    
    # Produce an extra matrix for each 'migration type' (e.g. prisoner-transfers, migration mixing, HIV infection, etc.).
    row_offset = num_pops + 2
    for mid in xrange(num_migrations):
        ws_poptrans.write(row_offset, 0, 'Migration Type '+str(mid+1))
        makeConnectionMatrix(worksheet = ws_poptrans, at_row = row_offset, at_col = 0, labels = pop_labels_default, formula_labels = pop_labels_formula)
        row_offset += num_pops + 2
    
    ws_poptrans.set_column(0, 0, ws_poptrans_width)
        
#            ws_poptrans.write(row_id, 0, link_name)
#        ws_linkpars.write(row_id, 1, 'Assumption')
#        for k in xrange(len(data_tvec)):
#            ws_linkpars.write(row_id, k+offset_tvec, data_tvec[k])
#        for pid in xrange(num_pops):
#            row_id += 1
#            ws_linkpars.write(row_id, 0, "='%s'!%s" % (settings.databook['sheet_names']['pops'], rc(pid+1,0)), None, temp_pop_names[pid])
#            ws_linkpars.write(row_id, 1, '=IF(SUMPRODUCT(--(%s:%s<>""))=0,%f,"N.A.")' % (rc(row_id,offset_tvec), rc(row_id,offset_tvec+len(data_tvec)-1), def_val), None, def_val)
#            ws_linkpars.write(row_id, 2, 'OR')
#            
#            # Values for all extra populations default to first population values.
#            if pid > 0:
#                for k in xrange(len(data_tvec)):
#                    ws_linkpars.write(row_id, k+offset_tvec, '=IF(%s="","",%s)' % (rc(row_id-pid,k+offset_tvec),rc(row_id-pid,k+offset_tvec)), None, '')
    
    # Cascade parameters sheet.
    row_id = 0
    for link_name in settings.linkpar_name_labels.keys():
        link_label = settings.linkpar_name_labels[link_name]
        def_val = 0
        if 'default' in settings.linkpar_specs[link_label]:
            def_val = settings.linkpar_specs[link_label]['default']
        
        # Make the parameter-specific data-entry block.
        makeValueEntryArrayBlock(worksheet = ws_linkpars, at_row = row_id, at_col = 1, num_arrays = num_pops, tvec = data_tvec, assumption = def_val)        
        
        # Make the parameter-specific population references.
        ws_linkpars.write(row_id, 0, link_name)
        for pid in xrange(num_pops):
            row_id += 1
            ws_linkpars.write(row_id, 0, "='%s'!%s" % (settings.databook['sheet_names']['pops'], rc(pid+1,0)), None, pop_names_default[pid])
        
        row_id += 2

    ws_linkpars.set_column(0, 0, ws_linkpars_width)
    ws_linkpars.set_column(1, 2, assumption_width)
    
    workbook.close()
    

#%% Function for loading project databook

# NOTE: This needs so much quality-assurance testing. Need to ensure that input data sheet aligns with cascade settings.
def loadSpreadsheetFunc(settings, databook_path = None):
    ''' Load data spreadsheet into Project data dictionary. '''
    
    if databook_path is None: databook_path = default_path
    try: workbook = xlrd.open_workbook(databook_path)
    except: raise OptimaException('ERROR: Project data workbook was unable to be loaded from... %s' % databook_path)
    ws_pops = workbook.sheet_by_name(settings.databook['sheet_names']['pops'])
    ws_poptrans = workbook.sheet_by_name(settings.databook['sheet_names']['poptrans'])
    ws_linkpars = workbook.sheet_by_name(settings.databook['sheet_names']['linkpars'])

    # Population names sheet.
    data = odict()
    data['pops'] = odict()
    data['pops']['name_labels'] = odict()
    data['pops']['label_names'] = odict()      # A reverse of the previous dictionary.
    data['pops']['ages'] = odict()
    for row_id in xrange(ws_pops.nrows):
        if row_id > 0:
            if ws_pops.cell_value(row_id, 0) not in ['']:
                if ws_pops.cell_value(row_id, 1) not in ['']:
                    pop_label = str(ws_pops.cell_value(row_id, 1))
                else:
                    pop_label = 'pop' + str(row_id)
                data['pops']['name_labels'][str(ws_pops.cell_value(row_id, 0))] = pop_label
                data['pops']['label_names'][pop_label] = str(ws_pops.cell_value(row_id, 0))
                age_min = ws_pops.cell_value(row_id, 2)
                age_max = ws_pops.cell_value(row_id, 3)
                if isinstance(age_min, Number) and isinstance(age_max, Number):
                    data['pops']['ages'][pop_label] = {'min':float(age_min), 'max':float(age_max), 'range':1+float(age_max)-float(age_min)}

    # Inter-population transitions sheet.
    # Migration matrices must be divided from each other by an empty row.
    # NOTE: Needs some way for users to know that only the first filled-in cell per row will be noted for aging.
    data['pops']['age_trans'] = odict()
    mig_specified = False
    mig_type = None
    for row_id in xrange(ws_poptrans.nrows):
        zero_col = ws_poptrans.cell_value(row_id, 0)
        
        # If parser finds an empty row (technically empty first cell) migration-parsing resets, ready for the next custom type.
        if mig_specified and zero_col in ['']:
            mig_specified = False
            
        # Parse through migration matrix.
        if mig_specified:
            pop_source = str(zero_col)
            for col_id in xrange(ws_poptrans.ncols):
                if col_id > 0:
                    val = ws_poptrans.cell_value(row_id, col_id)
                    if val in ['y']:
                        pop_sink = str(ws_poptrans.cell_value(0, col_id))
                        if 'range' not in data['pops']['ages'][pop_source].keys():
                            raise OptimaException('ERROR: An age transition has been flagged for a source population group with no age range.')
                        data['pops']['age_trans'][pop_source] = pop_sink
                        break   # Only the first tag in a row gets counted currently!
        
        # First row after a blank one must contain the new migration type as its first element. Parser re-activated.
        if not mig_specified and zero_col not in ['']:
            mig_specified = True
            mig_type = str(zero_col)
    
    # Cascade parameters sheet.
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
            # NOTE: Somewhat hard-coded. Improve.
            list_t = []
            list_y = []
            for col_id in xrange(ws_linkpars.ncols):
                if col_id == 1:
                    data['linkpars'][current_linkpar_label][current_pop_label]['y_format'] = str(ws_linkpars.cell_value(row_id, col_id))
                if col_id > 1 and isinstance(ws_linkpars.cell_value(row_id, col_id), Number):
                    list_y.append(float(ws_linkpars.cell_value(row_id, col_id)))
                    if not isinstance(ws_linkpars.cell_value(row_id-1-pop_id, col_id), Number):
                        list_t.append(float(ws_linkpars.cell_value(row_id-1-pop_id, col_id+2)))
                        break
                    else:
                        list_t.append(float(ws_linkpars.cell_value(row_id-1-pop_id, col_id)))
            data['linkpars'][current_linkpar_label][current_pop_label]['t'] = np.array(list_t)
            data['linkpars'][current_linkpar_label][current_pop_label]['y'] = np.array(list_y)                
            
            pop_id += 1
            
    return data