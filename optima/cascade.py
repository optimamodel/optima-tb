#%% Imports

from utils import odict, OptimaException, flattenDict
import settings as project_settings

import logging
logger = logging.getLogger(__name__)

import xlrd
import numpy as np
import pylab as pl
from copy import deepcopy as dcp



#%% Function to convert a cascade workbook into a framework to store in project settings

def loadCascadeSettingsFunc(cascade_path, settings):
    ''' Generates node and link settings based on cascade spreadsheet. 
    
    Note: at the moment, this function is intended to be called from settings.py, and takes in settings as a parameter, 
    updates the relevants fields and then returns settings. This isn't the most elegant way, but was better than returning 11 parameters.
    It can - and should - be improved.
    Also, do not call this function directly as the cascade reset is not part of this function. Unexpected behaviour is guaranteed.
    '''    
    
    try: workbook = xlrd.open_workbook(cascade_path)
    except: raise OptimaException('ERROR: Cannot find cascade workbook from which to load model structure.')
    ws_nodes =      workbook.sheet_by_name('Compartments')
    ws_links =      workbook.sheet_by_name('Transitions')
    ws_characs =    workbook.sheet_by_name('Characteristics')
    ws_pars =       workbook.sheet_by_name('Parameters')
    try: ws_sheetnames = workbook.sheet_by_name('Databook Sheet Names')
    except: ws_sheetnames = None
    
    
    #%% Optional sheet: Databook Sheet Names
    # This worksheet, if it exists, defines custom sheets to produce when generating a project databook.
    
    if ws_sheetnames:
        cid_label = None
        cid_name = None
        for col_id in xrange(ws_sheetnames.ncols):
            if ws_sheetnames.cell_value(0, col_id) == 'Sheet Label': cid_label = col_id
            if ws_sheetnames.cell_value(0, col_id) == 'Sheet Name': cid_name = col_id
                
        for row_id in xrange(ws_sheetnames.nrows):
            if row_id > 0 and ws_sheetnames.cell_value(row_id, cid_label) not in ['']:
                sheet_label = str(ws_sheetnames.cell_value(row_id, cid_label))
                sheet_name = str(ws_sheetnames.cell_value(row_id, cid_name))
                if sheet_label in settings.databook['sheet_names'].keys():
                    raise OptimaException('ERROR: Custom sheet label "%s" is already used as a standard label when generating project databooks.' % sheet_label)
                settings.databook['custom_sheet_names'][sheet_label] = sheet_name        
    
    #%% First core sheet: Compartments
    # This worksheet defines nodes/compartments of the cascade network, with label and name required at a bare minimum.
    
    cid_label = None
    cid_name = None
    cid_coords = None
    cid_no_plot = None
    cid_birth = None
    cid_dead = None
    cid_junction = None
    for col_id in xrange(ws_nodes.ncols):
        if ws_nodes.cell_value(0, col_id) == 'Code Label': cid_label = col_id
        if ws_nodes.cell_value(0, col_id) == 'Full Name': cid_name = col_id
        if ws_nodes.cell_value(0, col_id) == 'Plot Coordinates': cid_coords = col_id
        if ws_nodes.cell_value(0, col_id) == 'No Plot': cid_no_plot = col_id
        if ws_nodes.cell_value(0, col_id) == 'Birth Tag': cid_birth = col_id
        if ws_nodes.cell_value(0, col_id) == 'Dead Tag': cid_dead = col_id
        if ws_nodes.cell_value(0, col_id) == 'Junction': cid_junction = col_id
    if None in [cid_label, cid_name]:
        raise OptimaException('ERROR: Cascade compartment worksheet does not have correct column headers.')
    
    std_norm_nodes = []    # A list for compartments that can be involved in inter-population transfers.
    
    # Append node labels and full names to relevant odicts and lists. Labels are crucial.
    for row_id in xrange(ws_nodes.nrows):
        if row_id > 0 and ws_nodes.cell_value(row_id, cid_label) not in ['']:
            node_label = str(ws_nodes.cell_value(row_id, cid_label))
            settings.node_specs[node_label] = odict()
            settings.node_names.append(str(ws_nodes.cell_value(row_id, cid_name)))
            good_for_transfer = True    # A flag for whether a compartment is a valid node to 'migrate' from/into.
            
            # Store optional graph coordinates for this node in the cascade network. Used when plotting schematic.
            if not cid_coords is None:
                val = str(ws_nodes.cell_value(row_id, cid_coords))
                if val not in ['']:
                    try: coords = val.strip('()').split(',')
                    except: raise OptimaException('ERROR: Plot Coordinates for "%s" are not written in "(x,y)" format.' % node_label)
                    settings.node_specs[node_label]['coords'] = (float(coords[0]),float(coords[1]))
            
            # Store optional information about whether to avoid plotting node in the schematic.
            if not cid_no_plot is None:
                val = str(ws_nodes.cell_value(row_id, cid_no_plot))
                if val not in ['']:
                    settings.node_specs[node_label]['tag_no_plot'] = val
            
            # Store optional information about whether this node is a compartment representing births.
            if not cid_birth is None:
                val = str(ws_nodes.cell_value(row_id, cid_birth))
                if val not in ['']:
                    settings.node_specs[node_label]['tag_birth'] = val
                    good_for_transfer = False

            # Store optional information about whether this node is a compartment for the dead.
            if not cid_dead is None:
                val = str(ws_nodes.cell_value(row_id, cid_dead))
                if val not in ['']:
                    settings.node_specs[node_label]['tag_dead'] = val
                    good_for_transfer = False

            # Optionally note down whether the compartment is just a junction.
            # Junctions empty themselves via outflows at the end of each model timestep.
            if not cid_junction is None:
                val = str(ws_nodes.cell_value(row_id, cid_junction))
                if val not in ['']:
                    settings.node_specs[node_label]['junction'] = val
                    settings.junction_labels.append(node_label)
                    
            if good_for_transfer:
                std_norm_nodes.append(node_label)
            
#    print std_norm_nodes
                    
    


    #%% Second core sheet: Transitions
    # This worksheet describes links/transitions of the cascade network with respect to the nodes/compartments.

    # Validate bijectivity of compartments defined in previous sheet and those that mark the edges of the adjacency matrix in this sheet.
    test = {}
    for row_id in xrange(ws_links.nrows):
        val = str(ws_links.cell_value(row_id, 0))
        if row_id > 0 and val not in ['']:
            if val not in settings.node_specs.keys():
                raise OptimaException('ERROR: Cascade transitions worksheet has a row header ("%s") that is not a known compartment code label.' % val)
            test[val] = True
    for node_label in settings.node_specs.keys():
        if node_label not in test.keys():
            raise OptimaException('ERROR: Compartment code label "%s" is not represented in row headers of transitions worksheet.' % node_label)
    test = {}
    for col_id in xrange(ws_links.ncols):
        val = str(ws_links.cell_value(0, col_id))
        if col_id > 0 and val not in ['']:
            if val not in settings.node_specs.keys():
                raise OptimaException('ERROR: Cascade transitions worksheet has a column header ("%s") that is not a known compartment code label.' % val)
            test[val] = True
    for node_label in settings.node_specs.keys():
        if node_label not in test.keys():
            raise OptimaException('ERROR: Compartment code label "%s" is not represented in column headers of transitions worksheet.' % node_label)
    
    # Store linked node tuples by tag.
    for row_id in xrange(ws_links.nrows):
        for col_id in xrange(ws_links.ncols):
            if row_id > 0 and col_id > 0:
                n1 = str(ws_links.cell_value(row_id, 0))
                n2 = str(ws_links.cell_value(0, col_id))
                val = str(ws_links.cell_value(row_id, col_id))
                if not '' in [n1,n2,val]:
                    if not val in settings.links.keys():
                        settings.links[val] = []
                    settings.links[val].append((n1, n2))
                    
                
    #%% Third core sheet: Characteristics
    # This worksheet describes characteristics of the model (i.e. they can be considered as outputs).
    # A label, name and set of included compartments is required.
                
    # Sweep through column headers to make sure the right tags exist. Basically checking spreadsheet format.
    cid_label = None
    cid_name = None
    cid_sheet = None
    cid_plot = None
    cid_percentage = None
    cid_denom = None
    cid_order = None
    cid_default = None
    cid_entry = None
    cid_include_start = None
    cid_include_end = None
    cid_yfactor = None
    for col_id in xrange(ws_characs.ncols):
        if ws_characs.cell_value(0, col_id) == 'Code Label': cid_label = col_id
        if ws_characs.cell_value(0, col_id) == 'Full Name': cid_name = col_id
        if ws_characs.cell_value(0, col_id) == 'Sheet Label': cid_sheet = col_id
        if ws_characs.cell_value(0, col_id) == 'Plot Value': cid_plot = col_id
        if ws_characs.cell_value(0, col_id) == 'Plot Percentage': cid_percentage = col_id
        if ws_characs.cell_value(0, col_id) == 'Denominator': cid_denom = col_id
        if ws_characs.cell_value(0, col_id) == 'Databook Order': cid_order = col_id
        if ws_characs.cell_value(0, col_id) == 'Default Value': cid_default = col_id
        if ws_characs.cell_value(0, col_id) == 'Entry Point': cid_entry = col_id
        if ws_characs.cell_value(0, col_id) == 'Includes': cid_include_start = col_id
        if ws_characs.cell_value(0, col_id) == 'Calibrate?': cid_yfactor = col_id
    
    
    
    # Work out where the 'include' columns end when defining cascade characteristics.
    cid_list = np.array(sorted([cid_label, cid_name, cid_sheet, cid_plot, cid_percentage, cid_denom, cid_order, cid_default, cid_entry, ws_characs.ncols]))
    cid_include_end = cid_list[sum(cid_include_start > cid_list)] - 1
    
    if None in [cid_label, cid_name, cid_include_start, cid_include_end]:
        raise OptimaException('ERROR: Cascade characteristics worksheet does not have correct column headers.')
    
    custom_std_norm_charac = False    
    
    # Actually append characteristic labels and full names to relevant odicts and lists. Labels are crucial.
    entry_dict = {}
    standard_sheet_count = ws_characs.nrows - 1      # All characteristics are printed to the standard databook sheet to begin with.
    for row_id in xrange(ws_characs.nrows):
        if row_id > 0 and ws_characs.cell_value(row_id, cid_label) not in ['']:
            charac_label = str(ws_characs.cell_value(row_id, cid_label))
            if charac_label == settings.charac_std_norm: raise OptimaException('ERROR: Characteristic cannot be named %s as that is a reserved special label.' % settings.charac_std_norm)
            charac_name = str(ws_characs.cell_value(row_id, cid_name))
            
            settings.charac_name_labels[charac_name] = charac_label
            settings.charac_specs[charac_label] = odict()
            settings.charac_specs[charac_label]['name'] = charac_name
            settings.charac_specs[charac_label]['includes'] = []
            settings.charac_specs[charac_label]['sheet_label'] = 'charac'    # Default label of databook sheet to print to.
            
            # Attribute a custom sheet label to this characteristic if available.
            if not cid_sheet is None:
                val = str(ws_characs.cell_value(row_id, cid_sheet))
                if val not in ['']:
                    if not val in settings.databook['custom_sheet_names'].keys():
                        raise OptimaException('ERROR: Project databook sheet-label "%s" for characteristic "%s" has not been previously defined.' % (val, charac_label))
                    settings.charac_specs[charac_label]['sheet_label'] = val
                    standard_sheet_count -= 1
            
            # Work out which compartments/characteristics this characteristic is a sum of.
            for col_id_inc in xrange(cid_include_end - cid_include_start + 1):
                col_id = cid_include_start + col_id_inc
                val = str(ws_characs.cell_value(row_id, col_id))
                if val not in ['']:
                    if val not in settings.node_specs.keys() + settings.charac_specs.keys()[:-1]:
                        raise OptimaException('ERROR: Cascade characteristic "%s" is being defined with an inclusion reference to "%s", which has not been defined yet.' % (charac_label, val))
                    settings.charac_specs[charac_label]['includes'].append(val)
            
            flat_list, dep_list = flattenDict(input_dict = settings.charac_specs, base_key = charac_label, sub_keys = ['includes'], limit = settings.recursion_limit)
            if len(flat_list) != len(set(flat_list)):
                raise OptimaException('ERROR: Cascade characteristic "%s" contains duplicate references to a compartment (in its numerator) when all the recursion is flattened out.' % (charac_label))
            numerator_ref_list = dcp(flat_list)   # Useful for checking what compartments the numerator characteristic references.
            
            # If a user-defined characteristic includes all transfer-enabled compartments, overwrite the settings reference to it.
            if set(flat_list) == set(std_norm_nodes):
                settings.charac_std_norm = charac_label
                custom_std_norm_charac = True
            
            # Work out which compartment/characteristic this characteristic is normalised by.
            if not cid_denom is None:
                val = str(ws_characs.cell_value(row_id, cid_denom))
                if val not in ['']:
                    if val not in settings.node_specs.keys() + settings.charac_specs.keys()[:-1]:
                        raise OptimaException('ERROR: Cascade characteristic "%s" is being defined with reference to denominator "%s", which has not been defined yet.' % (charac_label, val))
                    settings.charac_specs[charac_label]['denom'] = val
            
            # Store whether we should plot this value as a default or not
            if not cid_plot is None:
                val = str(ws_characs.cell_value(row_id, cid_plot))
                if val not in ['']:
                    settings.charac_specs[charac_label]['plot_characteristic'] = val
            
            # Store whether characteristic should be converted to percentages when plotting.
            if not cid_percentage is None:
                val = str(ws_characs.cell_value(row_id, cid_percentage))
                if val not in ['']:
                    settings.charac_specs[charac_label]['plot_percentage'] = val
            
            # Store order that characteristics should be printed in project databook.
            # Defaults to high value (i.e. low priority), even if column does not exist.
            settings.charac_specs[charac_label]['databook_order'] = ws_characs.nrows + 1   # Any value missing from a databook order column means their printouts are lowest priority.
            if not cid_order is None:
                val = ws_characs.cell_value(row_id, cid_order)
                if val not in ['']:
                    settings.charac_specs[charac_label]['databook_order'] = int(val)
                    if int(val) < 0:
                        standard_sheet_count -= 1   # Unprinted characteristics do not end up on the standard sheet.
                    
            # Store characteristic default value if available.
            if not cid_default is None:
                val = str(ws_characs.cell_value(row_id, cid_default))
                if val not in ['']:
                    settings.charac_specs[charac_label]['default'] = float(val)
                    
            # Store characteristic entry point if available.
            # Without this, the characteristic will not be converted into an initial value for the model.
            if not cid_entry is None:
                val = str(ws_characs.cell_value(row_id, cid_entry))
                if val not in ['']:
                    settings.charac_specs[charac_label]['entry_point'] = val
                    
                    if settings.charac_specs[charac_label]['databook_order'] < 0:
                        raise OptimaException('ERROR: Characteristic "%s" cannot be used to seed a model (i.e. have an entry point) if it does not feature in the databook (i.e. if it has a negative databook order).' % (charac_label))
                    if 'denom' in settings.charac_specs[charac_label].keys():
                        denom_label = settings.charac_specs[charac_label]['denom']
                        if denom_label not in settings.charac_specs.keys()[:-1] or 'entry_point' not in settings.charac_specs[denom_label].keys():
                            raise OptimaException('ERROR: At this time, characteristic "%s" cannot be used to seed a model (i.e. have an entry point) if its denominator "%s" is not a model-seeding characteristic.' % (charac_label, denom_label))
                    if val in entry_dict.keys():
                        raise OptimaException('ERROR: Cascade characteristic "%s" is not the first to use "%s" as an entry point.' % (charac_label, val))
                    if val not in settings.node_specs.keys():
                        raise OptimaException('ERROR: There is no compartment named "%s" that can be used as an entry point by cascade characteristic "%s".' % (val, charac_label))                        
                    numerator_ref_list = set(numerator_ref_list)
                    try: numerator_ref_list.remove(val)
                    except: raise OptimaException('ERROR: Entry point "%s" for characteristic "%s" must be a compartment it includes (in its numerator), either directly or via characteristic reference.' % (val, charac_label))
                    entry_dict[val] = dcp(list(numerator_ref_list))
                    
            # Store whether we should use calibrate the initial value or not (only applicable when there's also an entry point)
            if not cid_yfactor is None:
                val = str(ws_characs.cell_value(row_id, cid_yfactor))
                if val.lower() == 'n' or val == '-1':
                    settings.charac_specs[charac_label]['y_factor'] = project_settings.DO_NOT_SCALE
                elif val not in ['']:
                    settings.charac_specs[charac_label]['y_factor'] = float(val)
                else:
                    settings.charac_specs[charac_label]['y_factor'] = project_settings.DEFAULT_YFACTOR
        
        # Make sure empty space rows do not get counted when deciding if there are characteristics left to populate a default databook sheet.
        elif row_id > 0:
            standard_sheet_count -= 1
                    
    # If a characteristic that includes all transfer-enabled compartments is not defined, create one.
    if not custom_std_norm_charac:
        settings.charac_name_labels[settings.charac_std_norm_name] = settings.charac_std_norm
        settings.charac_specs[settings.charac_std_norm] = odict()
        settings.charac_specs[settings.charac_std_norm]['name'] = settings.charac_std_norm_name
        settings.charac_specs[settings.charac_std_norm]['includes'] = std_norm_nodes
        settings.charac_specs[settings.charac_std_norm]['databook_order'] = -1      # This special 'standard' characteristic is not initialisable or used for calibration unless user-defined.
                    

    
    
    # If all characteristics in this sheet are to be printed to custom databook sheets, no need for a default.
    if standard_sheet_count <= 0:
        settings.make_sheet_characs = False
                        
    # Ensuring that no two characteristics with entry-points include each other's entry-points (or more complex referencing cycles).
    # This allows for disaggregated characteristic values to be partitioned from aggregate characteristics during model-seeding.
    for entry_label in entry_dict.keys():
        try:
            flattenDict(input_dict = entry_dict, base_key = entry_label, limit = settings.recursion_limit)
        except OptimaException:
            raise OptimaException('Characteristic "%s" references an entry point for another characteristic that, via some chain, eventually references the entry point of characteristic "%s". Alternatively, maximum recursion depth is set too small.' % (entry_label, entry_label))
    
        
    #%% Fourth core sheet: Parameters
    # Sweep through column headers to make sure the right tags exist. Basically checking spreadsheet format.
    cid_tag = None
    cid_label = None
    cid_name = None
    cid_sheet = None
    cid_default = None
    cid_order = None
    cid_function = None
    cid_yfactor = None
    for col_id in xrange(ws_pars.ncols):
        if ws_pars.cell_value(0, col_id) == 'Transition Tag': cid_tag = col_id
        if ws_pars.cell_value(0, col_id) == 'Code Label': cid_label = col_id
        if ws_pars.cell_value(0, col_id) == 'Full Name': cid_name = col_id
        if ws_pars.cell_value(0, col_id) == 'Sheet Label': cid_sheet = col_id
        if ws_pars.cell_value(0, col_id) == 'Default Value': cid_default = col_id
        if ws_pars.cell_value(0, col_id) == 'Databook Order': cid_order = col_id
        if ws_pars.cell_value(0, col_id) == 'Function': cid_function = col_id
        if ws_pars.cell_value(0, col_id) == 'Calibrate': cid_yfactor = col_id
    if None in [cid_tag, cid_label, cid_name]:
        raise OptimaException('ERROR: Cascade transition-parameters worksheet does not have correct column headers.')
    
    # Store transition details in settings and make sure there is tag bijectivity between this sheet and the transition matrix.
    standard_sheet_count = ws_pars.nrows - 1      # All parameters are printed to the standard databook sheet to begin with.
    for row_id in xrange(ws_pars.nrows):
        tag = str(ws_pars.cell_value(row_id, cid_tag))
        label = str(ws_pars.cell_value(row_id, cid_label))
        name = str(ws_pars.cell_value(row_id, cid_name))
        if row_id > 0 and label not in ['']:
            settings.linkpar_specs[label] = {'name':name, 'sheet_label':'linkpars'}
            settings.linkpar_name_labels[name] = label
            if tag not in ['']:
                if tag not in settings.links:
                    raise OptimaException('ERROR: Cascade transition-parameter worksheet has a tag (%s) that is not in the transition matrix.' % tag)
                settings.linkpar_specs[label]['tag'] = tag
            else:
                settings.par_deps[label] = True     # Untagged parameters are dependencies for actual tagged transitions.
            
            # Attribute a custom sheet label to this parameter if available.
            if not cid_sheet is None:
                val = str(ws_pars.cell_value(row_id, cid_sheet))
                if val not in ['']:
                    if not val in settings.databook['custom_sheet_names'].keys():
                        raise OptimaException('ERROR: Project databook sheet-label "%s" for parameter "%s" has not been previously defined.' % (val, label))
                    settings.linkpar_specs[label]['sheet_label'] = val
                    standard_sheet_count -= 1
            
            # Store order that cascade parameters should be printed in project databook.
            settings.linkpar_specs[label]['databook_order'] = ws_pars.nrows+1   # Any value missing from a databook order column means their printouts are lowest priority.
            if not cid_order is None:
                val = ws_pars.cell_value(row_id, cid_order)
                if val not in ['']:
                    settings.linkpar_specs[label]['databook_order'] = int(val)
                    if int(val) < 0:    
                        standard_sheet_count -= 1   # Unprinted parameters do not end up on the standard sheet.
            
            # Store parameter default value if available.
            if not cid_default is None:
                val = str(ws_pars.cell_value(row_id, cid_default))
                if val not in ['']:
                    settings.linkpar_specs[label]['default'] = float(val)
                    
            # Store the token stack corresponding to a custom function that defines this parameter, if available.
            if not cid_function is None:
                val = str(ws_pars.cell_value(row_id, cid_function))
                if val not in ['']:
                    if 'default' in settings.linkpar_specs[label]:
                        raise OptimaException('ERROR: Parameter "%s" is a custom function of other parameters and characteristics. Specifying a default is thus restricted, so as to avoid user confusion.' % label)
                    if settings.linkpar_specs[label]['databook_order'] >= 0:
                        raise OptimaException('ERROR: Parameter "%s" is a custom function of other parameters and characteristics. Tag "Databook Order" column with a negative number so that conflicts with user-provided values do not arise.' % label)
                    settings.par_funcs[label] = True
                    expr_stack, var_dict = settings.parser.produceStack(val)
                    settings.linkpar_specs[label]['f_stack'] = expr_stack
                    settings.linkpar_specs[label]['deps'] = var_dict
                    for var in var_dict.keys():
                        if not var in settings.charac_specs.keys():
                            if not var in settings.linkpar_specs.keys():
                                raise OptimaException('ERROR: Dependency "%s" has not been defined by the time "%s" is loaded into settings.' % (var, label))
                        else:
                            if not var in settings.charac_deps.keys():
                                flat_list, dep_list = flattenDict(input_dict = settings.charac_specs, base_key = var, sub_keys = ['includes','denom'], limit = settings.recursion_limit)
                                for dep in dep_list:    # NOTE: Dependencies are presumably provided in proper order by flattenDict due to recursion. Could a counterexample be engineered that breaks the assumption...? 
                                    settings.charac_specs[dep]['par_dependency'] = True
                                    settings.charac_deps[dep] = True
                                    
            # Store whether we should use calibrate the initial value or not         
            if not cid_yfactor is None:
                val = str(ws_pars.cell_value(row_id, cid_yfactor))
                if val.lower() == 'n' or val == '-1':
                    settings.linkpar_specs[label]['y_factor'] = project_settings.DO_NOT_SCALE
                elif val not in ['']:
                    settings.linkpar_specs[label]['y_factor'] = float(val)
                else:
                    settings.linkpar_specs[label]['y_factor'] = project_settings.DEFAULT_YFACTOR
                    
        # Make sure empty space rows do not get counted when deciding if there are parameters left to populate a default databook sheet.
        elif row_id > 0:
            standard_sheet_count -= 1
    
    # If the default/overwritten population-count characteristic is not a dependency by now, along with its own dependencies, make it one at the end of the charac_deps odict.
    if settings.charac_std_norm not in settings.charac_deps.keys():
        flat_list, dep_list = flattenDict(input_dict = settings.charac_specs, base_key = settings.charac_std_norm, sub_keys = ['includes','denom'], limit = settings.recursion_limit)
        for dep in dep_list:    # NOTE: Dependencies are presumably provided in proper order by flattenDict due to recursion. Could a counterexample be engineered that breaks the assumption...? 
            settings.charac_specs[dep]['par_dependency'] = True
            settings.charac_deps[dep] = True
        settings.charac_specs[settings.charac_std_norm]['par_dependency'] = True
        settings.charac_deps[settings.charac_std_norm] = True
        
#    print settings.charac_deps
    
    # If all parameters in this sheet are to be printed to custom databook sheets, no need for a default.
    if standard_sheet_count <= 0:
        settings.make_sheet_linkpars = False
    
    # Final validations.            
    for tag in settings.links.keys():
        if tag not in [x['tag'] for x in settings.linkpar_specs[:] if 'tag' in x]:
            raise OptimaException('ERROR: Transition matrix tag "%s" is not represented in transition-parameter worksheet.' % tag)
    
    test_labels = settings.linkpar_specs.keys() + settings.charac_specs.keys()       
    if len(test_labels) != len(set(test_labels)):
        raise OptimaException('ERROR: Cascade workbook appears to have duplicate characteristic/parameter code labels.')
    
    test_names = settings.linkpar_name_labels.keys() + settings.charac_name_labels.keys()
    if len(test_names) != len(set(test_names)):
        raise OptimaException('ERROR: Cascade workbook appears to have duplicate characteristic/parameter full names.')
    


#%% Function to plot a cascade framework loaded into settings
    
def plotCascadeFunc(settings):
    
    import networkx as nx
    
    fig, ax = pl.subplots(figsize=(10,10))
    G = nx.DiGraph()
    plottable_nodes = [nid for nid in settings.node_specs.keys() if 'tag_no_plot' not in settings.node_specs[nid]]
    plottable_links = [link for lid in settings.links for link in settings.links[lid] if (link[0] in plottable_nodes and link[1] in plottable_nodes)]
    G.add_nodes_from(plottable_nodes)
    G.add_edges_from(plottable_links)

    # Use plot coordinates if stored and arrange the rest of the cascade out in a unit circle.
    pos = {}
    num_nodes = len(plottable_nodes)
    k = 0
    for node in plottable_nodes:
        try: pos[node] = (settings.node_specs[node]['coords'][0], settings.node_specs[node]['coords'][1])
        except: pos[node] = (np.sin(2.0*np.pi*k/num_nodes), np.cos(2.0*np.pi*k/num_nodes))
        k += 1
    
    # Generate edge label dictionary with tags from spreadsheet.
    el = {}
    for par_name in settings.linkpar_specs.keys():
        if 'tag' in settings.linkpar_specs[par_name]:
            for link in settings.links[settings.linkpar_specs[par_name]['tag']]:
                el[link] = settings.linkpar_specs[par_name]['tag']

    nx.draw_networkx_nodes(G, pos, node_shape = 'o', nodelist = [x for x in G.nodes() if not 'junction' in settings.node_specs[x].keys()], node_size = 1250, node_color = 'w')
    nx.draw_networkx_nodes(G, pos, node_shape = 's', nodelist = [x for x in G.nodes() if 'junction' in settings.node_specs[x].keys()], node_size = 750, node_color = 'w')
    ax.axis('tight')
    nx.draw_networkx_labels(G, pos)
    nx.draw_networkx_edges(G, pos)
#        nx.draw_networkx_edge_labels(G, pos, edge_labels = el, label_pos = 0.25, font_size = 14)
    
    [sp.set_visible(False) for sp in ax.spines.values()]
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Cascade Schematic')
    pl.show()