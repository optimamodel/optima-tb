#%% Imports

from utils import odict, OptimaException

import xlrd
import xlsxwriter as xw
import numpy as np
from xlsxwriter.utility import xl_rowcol_to_cell as rc
from numbers import Number
from copy import deepcopy as dcp



#%% Utility functions to generate sub-blocks of the project databook

def makeValueEntryArrayBlock(worksheet, at_row, at_col, num_arrays, tvec, assumption = 0.0, data_formats = None, print_conditions = None):
    '''
    Create a block where users choose data-entry format and enter values either as an assumption or time-dependent array.
    
    Args:
        worksheet           -   The worksheet in which to produce the block.
        at_row              -   The row at which the top-left corner of the block will be fixed.
        at_col              -   The column at which the top-left corner of the block will be fixed.
        num_arrays          -   The number of input arrays to generate.
        tvec                -   The year range for which values are (optionally) to be entered.
        assumption          -   The default value to enter in the assumption column.
        data_formats        -   A list of formats that data can be entered as.
                                This is not data validation, but will determine the form of calculations during model processing.
        print_conditions    -   A list of Excel-string conditions that are used to test whether data-entry arrays should be shown.
                                List must be of num_arrays length, but can include values of None to allow default printing behaviour for certain rows.
    '''
    
    if data_formats is None: data_formats = ['Fraction']#, 'Number']#, 'Probability', 'Number']
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
        if not print_conditions is None and not print_conditions[aid] is None:
            worksheet.write(row_id, at_col, '=IF(%s,"%s","")' % (print_conditions[aid], data_formats[0]), None, '')      # Default choice for data format.
            worksheet.write(row_id, at_col + 1, '=IF(%s,IF(SUMPRODUCT(--(%s:%s<>""))=0,%f,"N.A."),"")' % (print_conditions[aid], rc(row_id,offset), rc(row_id,offset+len(tvec)-1), assumption), None, '')
            worksheet.write(row_id, at_col + 2, '=IF(%s,"OR","")' % print_conditions[aid], None, '')
        else:
            worksheet.write(row_id, at_col, data_formats[0])      # Default choice for data format.
            worksheet.write(row_id, at_col + 1, '=IF(SUMPRODUCT(--(%s:%s<>""))=0,%f,"N.A.")' % (rc(row_id,offset), rc(row_id,offset+len(tvec)-1), assumption), None, assumption)
            worksheet.write(row_id, at_col + 2, 'OR')
        
#        # Make changes to the first row be mirrored in the other rows.
#        if aid > 0:
#            for k in xrange(len(tvec)):
#                worksheet.write(row_id, offset + k, '=IF(%s="","",%s)' % (rc(row_id-aid,offset+k),rc(row_id-aid,offset+k)), None, '')
               
               
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
def makeSpreadsheetFunc(settings, databook_path, num_pops = 5, num_migrations = 2):
    ''' Generate a data-input spreadsheet (e.g. for a country) corresponding to the loaded cascade settings. '''
    
    workbook = xw.Workbook(databook_path)
    ws_pops = workbook.add_worksheet(settings.databook['sheet_names']['pops'])
    ws_transmat = workbook.add_worksheet(settings.databook['sheet_names']['transmat'])
    ws_transval = workbook.add_worksheet(settings.databook['sheet_names']['transval'])
    
    # Regarding cascade parameters and characteristics, store sheets and corresponding row ids for further writing.
    ws_params = odict()
    for custom_label in settings.databook['custom_sheet_names']:
        ws_params[custom_label] = {'row_id':0, 'ws':workbook.add_worksheet(settings.databook['custom_sheet_names'][custom_label])}
    
    # Only produce standard sheets if there are any parameters left over that are not allocated to custom sheets.
    if settings.make_sheet_charac:
        ws_params['charac'] = {'row_id':0, 'ws':workbook.add_worksheet(settings.databook['sheet_names']['charac'])}
    if settings.make_sheet_linkpars:
        ws_params['linkpars'] = {'row_id':0, 'ws':workbook.add_worksheet(settings.databook['sheet_names']['linkpars'])}
    
    data_tvec = np.arange(settings.tvec_start, settings.tvec_end + 1.0/2)
    ws_pops_width = 15
    ws_transmat_width = 15
    ws_transval_width = 15
    name_width = 40
    assumption_width = 10
    
    #%% Population names sheet.
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
    
    #%% Inter-population transfers matrix sheet (from node to corresponding node).
    # Produce aging matrix.
    ws_transmat.write(0, 0, 'Aging')
    makeConnectionMatrix(worksheet = ws_transmat, at_row = 0, at_col = 0, labels = pop_labels_default, formula_labels = pop_labels_formula)
    
    # Produce an extra matrix for each 'migration type' (e.g. prisoner-transfers, migration mixing, HIV infection, etc.).
    # Store type names and formulae strings that reference them. Also store the rows at which migration matrices start.
    row_offset = num_pops + 2
    mig_types_default = []
    mig_types_formula = []
    mig_matrix_rows = []
    for mid in xrange(num_migrations):
        mig_type = 'Migration Type '+str(mid+1)
        mig_types_default.append(mig_type)
        mig_matrix_rows.append(row_offset)
        mig_types_formula.append("='%s'!%s" % (settings.databook['sheet_names']['transmat'], rc(row_offset,0)))
        ws_transmat.write(row_offset, 0, mig_type)
        makeConnectionMatrix(worksheet = ws_transmat, at_row = row_offset, at_col = 0, labels = pop_labels_default, formula_labels = pop_labels_formula)
        row_offset += num_pops + 2
    
    ws_transmat.set_column(0, 0, ws_transmat_width)
    
    #%% Inter-population transfer details sheet.
    row_id = 0
    for mid in xrange(num_migrations):
        ws_transval.write(row_id, 0, mig_types_formula[mid], None, mig_types_default[mid])
        print_conditions = []
        for k in xrange(num_pops*(num_pops-1)):
            print_conditions.append('%s<>"..."' % rc(row_id+k+1,0))
        makeValueEntryArrayBlock(worksheet = ws_transval, at_row = row_id, at_col = 3, num_arrays = num_pops*(num_pops-1), tvec = data_tvec, print_conditions = print_conditions)
        for source_id in xrange(num_pops):
            for target_id in xrange(num_pops):
                if source_id != target_id:
                    row_id += 1
                    r = mig_matrix_rows[mid] + source_id + 1
                    c = target_id + 1
                    ws_transval.write(row_id, 0, "=IF('%s'!%s=%s,%s,%s)" % (settings.databook['sheet_names']['transmat'],rc(r,c),'"y"',pop_names_formula[source_id][1:],'"..."'), None, '...')
                    ws_transval.write(row_id, 1, "=IF('%s'!%s=%s,%s,%s)" % (settings.databook['sheet_names']['transmat'],rc(r,c),'"y"','"--->"','""'), None, '')
                    ws_transval.write(row_id, 2, "=IF('%s'!%s=%s,%s,%s)" % (settings.databook['sheet_names']['transmat'],rc(r,c),'"y"',pop_names_formula[target_id][1:],'""'), None, '')
        
        row_id += 2
        
    ws_transval.set_column(0, 0, ws_transval_width)
    ws_transval.set_column(2, 2, ws_transval_width)
    ws_transval.set_column(3, 4, assumption_width)
    
    #%% Combine characteristics and parameters into one dictionary, then print out all elements to appropriate datasheets.
    all_specs = dcp(settings.charac_specs)
    all_specs.update(settings.linkpar_specs)
    specs_ordered = sorted(dcp(all_specs.keys()), key=lambda x: all_specs[x]['databook_order'])
    for def_label in specs_ordered:
        if not all_specs[def_label]['databook_order'] < 0:      # Do not print out characteristic/parameter if databook order is negative.        
            ws = ws_params[all_specs[def_label]['sheet_label']]['ws']
            row_id = ws_params[all_specs[def_label]['sheet_label']]['row_id']            
            
            def_name = all_specs[def_label]['name']
            default_val = 0.0
            if 'default' in all_specs[def_label]:
                default_val = all_specs[def_label]['default']
            
            # Make the data-entry blocks.
            # First handle 'count' and 'percentage' characteristics, as well as transitions for junctions.
            if def_label in settings.charac_specs.keys():
                if 'plot_percentage' in all_specs[def_label]:
                    data_formats = ['Fraction']
                else:
                    data_formats = ['Number']
            elif def_label in settings.linkpar_specs.keys():
                data_formats = None
                if 'tag' in settings.linkpar_specs[def_label]:
                    for pair in settings.links[settings.linkpar_specs[def_label]['tag']]:
                        if 'junction' in settings.node_specs[pair[0]].keys():
                            data_formats = ['Proportion']
            makeValueEntryArrayBlock(worksheet = ws, at_row = row_id, at_col = 1, num_arrays = num_pops, tvec = data_tvec, assumption = default_val, data_formats = data_formats)
            
            # Make the population references.
            ws.write(row_id, 0, def_name)
            for pid in xrange(num_pops):
                row_id += 1
                ws.write(row_id, 0, pop_names_formula[pid], None, pop_names_default[pid])
        
            row_id += 2
            ws_params[all_specs[def_label]['sheet_label']]['row_id'] = row_id

    # Adjust widths for all custom sheets.
    for ws_label in ws_params:
        ws = ws_params[ws_label]['ws']
        ws.set_column(0, 0, name_width)
        ws.set_column(1, 2, assumption_width)    
    
    
    workbook.close()
    

#%% Function for loading project databook

# NOTE: This needs so much quality-assurance testing. Need to ensure that input data sheet aligns with cascade settings.
def loadSpreadsheetFunc(settings, databook_path):
    ''' Load data spreadsheet into Project data dictionary. '''
    
    try: workbook = xlrd.open_workbook(databook_path)
    except: raise OptimaException('ERROR: Project data workbook was unable to be loaded from... %s' % databook_path)
    ws_pops = workbook.sheet_by_name(settings.databook['sheet_names']['pops'])
    ws_transmat = workbook.sheet_by_name(settings.databook['sheet_names']['transmat'])
    ws_transval = workbook.sheet_by_name(settings.databook['sheet_names']['transval'])
    
    # Regarding cascade parameters and characteristics, store sheets and corresponding row ids for further writing.
    ws_params = odict()
    for custom_label in settings.databook['custom_sheet_names'].keys():
        ws_params[custom_label] = workbook.sheet_by_name(settings.databook['custom_sheet_names'][custom_label])
    
    # Only produce standard sheets if there are any parameters left over that are not allocated to custom sheets.
    if settings.make_sheet_charac:
        ws_params['charac'] = workbook.sheet_by_name(settings.databook['sheet_names']['charac'])
    if settings.make_sheet_linkpars:
        ws_params['linkpars'] = workbook.sheet_by_name(settings.databook['sheet_names']['linkpars'])
#    ws_charac = workbook.sheet_by_name(settings.databook['sheet_names']['charac'])
#    ws_linkpars = workbook.sheet_by_name(settings.databook['sheet_names']['linkpars'])

    #%% Population names sheet.
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

    #%% Inter-population transitions sheet.
    # Migration matrices must be divided from each other by an empty row.
    data['transfers'] = odict()
    mig_specified = False
    mig_type = None
    for row_id in xrange(ws_transmat.nrows):
        zero_col = ws_transmat.cell_value(row_id, 0)
        
        # If parser finds an empty row (technically empty first cell) migration-parsing resets, ready for the next custom type.
        if mig_specified and zero_col in ['']:
            mig_specified = False
            
        # Parse through migration matrix.
        if mig_specified:
            pop_source = str(zero_col)
            for col_id in xrange(ws_transmat.ncols):
                if col_id > 0:
                    val = ws_transmat.cell_value(row_id, col_id)
                    if val in ['y']:
                        pop_target = str(ws_transmat.cell_value(0, col_id))
                        if pop_source not in data['transfers'][mig_type].keys():
                            data['transfers'][mig_type][pop_source] = odict()
                        data['transfers'][mig_type][pop_source][pop_target] = odict()
                        if mig_type == 'aging':
                            if 'range' not in data['pops']['ages'][pop_source].keys():
                                raise OptimaException('ERROR: An age transition has been flagged for a source population group with no age range.')
                            else:
                                data['transfers'][mig_type][pop_source][pop_target]['t'] = np.array([settings.tvec_start])
                                data['transfers'][mig_type][pop_source][pop_target]['y'] = np.array([float(1/data['pops']['ages'][pop_source]['range'])])
                                data['transfers'][mig_type][pop_source][pop_target]['y_format'] = 'Fraction'.lower()
                            if len(data['transfers'][mig_type][pop_source].keys()) > 1:
                                raise OptimaException('ERROR: There are too many outgoing "%s" transitions listed for population "%s".' % (mig_type,pop_source))
        
        # First row after a blank one must contain the new migration type as its first element. Parser re-activated.
        if not mig_specified and zero_col not in ['']:
            mig_specified = True
            mig_type = str(zero_col).lower().replace(' ','_')
            data['transfers'][mig_type] = odict()
    
    #%% Inter-population transfer details sheet.
    mig_specified = False
    mig_type = None
    array_id = 0
    for row_id in xrange(ws_transval.nrows):
        zero_col = ws_transval.cell_value(row_id, 0)
        
        # If parser finds an empty row (technically empty first cell) migration-parsing resets, ready for the next custom type.
        if mig_specified and zero_col in ['']:
            mig_specified = False
            
        # Parse through migration-specific data rows.
        if mig_specified:
            pop_source = str(zero_col)
            if pop_source not in ['','...']:
                pop_source_label = data['pops']['name_labels'][pop_source]
                pop_target = ws_transval.cell_value(row_id, 2)
                pop_target_label = data['pops']['name_labels'][pop_target]
                list_t = []
                list_y = []
                for col_id in xrange(ws_transval.ncols):
                    if col_id == 3:
                        data['transfers'][mig_type][pop_source_label][pop_target_label]['y_format'] = str(ws_transval.cell_value(row_id, col_id)).lower()
                    if col_id > 3 and isinstance(ws_transval.cell_value(row_id, col_id), Number):
                        val = ws_transval.cell_value(row_id, col_id)
                        
                        list_y.append(float(val))
                        if not isinstance(ws_transval.cell_value(row_id-1-array_id, col_id), Number):
                            list_t.append(float(ws_transval.cell_value(row_id-1-array_id, col_id+2)))
                            break
                        else:
                            list_t.append(float(ws_transval.cell_value(row_id-1-array_id, col_id)))
                                
                data['transfers'][mig_type][pop_source_label][pop_target_label]['t'] = np.array(list_t)
                data['transfers'][mig_type][pop_source_label][pop_target_label]['y'] = np.array(list_y)
            array_id += 1


        
        # First row after a blank one must contain the new migration type as its first element. Parser re-activated.
        if not mig_specified and zero_col not in ['']:
            mig_specified = True
            mig_type = str(zero_col).lower().replace(' ','_')
            array_id = 0


    
    #%% Gather data value inputs for epidemic characteristics and cascade parameters sheets, be they custom or standard.
    
    name_to_label = dcp(settings.charac_name_labels)
    name_to_label.update(settings.linkpar_name_labels)
    data['characs'] = odict()
    data['linkpars'] = odict()

    for ws_label in ws_params.keys():
        ws = ws_params[ws_label]
        current_def_name = None
        current_def_label = None
        data_label = None
        pop_id = 0
        for row_id in xrange(ws.nrows):
            val = str(ws.cell_value(row_id, 0))
            if val in ['']:
                current_def_name = None
                pop_id = 0
            elif current_def_name is None:
                current_def_name = val
                current_def_label = name_to_label[val]
                if current_def_label in settings.charac_specs.keys():
                    data_label = 'characs'
                elif current_def_label in settings.linkpar_specs.keys():
                    data_label = 'linkpars'
                else:
                    raise OptimaException('ERROR: Data entry in the "%s" sheet includes "%s", which has no defined characteristic/parameter specifications.' % (ws_label, current_def_name))
                data[data_label][current_def_label] = odict()
            else:
                current_pop_label = data['pops']['name_labels'][val]
                if current_pop_label != data['pops']['label_names'].keys()[pop_id]:
                    raise OptimaException('ERROR: Somewhere in the "%s" parameters sheet, populations are not ordered as in the population definitions sheet.' % data_label)
                data[data_label][current_def_label][current_pop_label] = odict()
                
                # Run through the rows beneath the year range, but only if there is not a number in the cell corresponding to assumption.
                # NOTE: Somewhat hard-coded. Improve.
                list_t = []
                list_y = []
                for col_id in xrange(ws.ncols):
                    if col_id == 1:
                        data[data_label][current_def_label][current_pop_label]['y_format'] = str(ws.cell_value(row_id, col_id)).lower()
                    if col_id > 1 and isinstance(ws.cell_value(row_id, col_id), Number):
                        list_y.append(float(ws.cell_value(row_id, col_id)))
                        if not isinstance(ws.cell_value(row_id-1-pop_id, col_id), Number):
                            list_t.append(float(ws.cell_value(row_id-1-pop_id, col_id+2)))
                            break
                        else:
                            list_t.append(float(ws.cell_value(row_id-1-pop_id, col_id)))
                data[data_label][current_def_label][current_pop_label]['t'] = np.array(list_t)
                data[data_label][current_def_label][current_pop_label]['y'] = np.array(list_y)                
                
                pop_id += 1
    
#    # Epidemic characteristics and cascade parameters sheet.
#    data_labels = ['characs', 'linkpars']
#    ws_list = [ws_charac, ws_linkpars]
#    name_to_label_list = [settings.charac_name_labels, settings.linkpar_name_labels]
#    for k in xrange(2):
#        data_label = data_labels[k]
#        ws = ws_list[k]
#        name_to_label = name_to_label_list[k]
#    
#        data[data_label] = odict()
#        current_def_name = None
#        current_def_label = None
#        pop_id = 0
#        for row_id in xrange(ws.nrows):
#            val = str(ws.cell_value(row_id, 0))
#            if val in ['']:
#                current_def_name = None
#                pop_id = 0
#            elif current_def_name is None:
#                current_def_name = val
#                current_def_label = name_to_label[val]
#                data[data_label][current_def_label] = odict()
#            else:
#                current_pop_label = data['pops']['name_labels'][val]
#                if current_pop_label != data['pops']['label_names'].keys()[pop_id]:
#                    raise OptimaException('ERROR: Somewhere in the "%s" parameters sheet, populations are not ordered as in the population definitions sheet.' % data_label)
#                data[data_label][current_def_label][current_pop_label] = odict()
#                
#                # Run through the rows beneath the year range, but only if there is not a number in the cell corresponding to assumption.
#                # NOTE: Somewhat hard-coded. Improve.
#                list_t = []
#                list_y = []
#                for col_id in xrange(ws.ncols):
#                    if col_id == 1:
#                        data[data_label][current_def_label][current_pop_label]['y_format'] = str(ws.cell_value(row_id, col_id)).lower()
#                    if col_id > 1 and isinstance(ws.cell_value(row_id, col_id), Number):
#                        list_y.append(float(ws.cell_value(row_id, col_id)))
#                        if not isinstance(ws.cell_value(row_id-1-pop_id, col_id), Number):
#                            list_t.append(float(ws.cell_value(row_id-1-pop_id, col_id+2)))
#                            break
#                        else:
#                            list_t.append(float(ws.cell_value(row_id-1-pop_id, col_id)))
#                data[data_label][current_def_label][current_pop_label]['t'] = np.array(list_t)
#                data[data_label][current_def_label][current_pop_label]['y'] = np.array(list_y)                
#                
#                pop_id += 1
                
    # All parameters must be defined whether they are in the project databook or not.
    for label in settings.linkpar_specs.keys():
        if label not in data['linkpars'].keys():
            def_format = 'Fraction'.lower()     # NOTE: Hard-coded format assumption. Improve at some stage when allowing for default formats.
            def_val = np.nan
            if 'default' in settings.linkpar_specs[label]:
                def_val = settings.linkpar_specs[label]['default']
            if 'f_stack' in settings.linkpar_specs[label]:
                def_val = np.nan
            else:
                print('WARNING: Project data sheet does not contain required cascade parameter "%s".\n         Using default format "%s" and default value %f.' % (label, def_format, def_val))
            data['linkpars'][label] = odict()
            for pop in data['pops']['label_names'].keys():
                data['linkpars'][label][pop] = odict()
                data['linkpars'][label][pop]['y_format'] = def_format
                data['linkpars'][label][pop]['t'] = np.array([settings.tvec_start])
                data['linkpars'][label][pop]['y'] = np.array([float(def_val)]) 
                
            
    return data