#%% Imports

from utils import flattenDict, odict, OptimaException
from parsing import FunctionParser

import xlrd
import pylab as pl
import numpy as np
from copy import deepcopy as dcp


#%% Settings class (for data that is effectively static per epidemic context)

class Settings(object):
    '''
    An object for storing cascade metadata (loaded from a cascade workbook) and general defaults.
    A cascade structure comprises of compartments/nodes and transitions/links that detail the flow of people through stages of a disease.
    In general use, the attributes of this object are static regardless of project/simulation.

    Notes:
    Settings must be loaded into a project during its creation, as it details the framework for dealing with data and running the model.
    '''    
    
    def __init__(self, cascade_path = '../data/cascade.xlsx'):
        self.tvec_start = 2000.0     # Default start year for data input and simulations.
        self.tvec_end = 2030.0       # Default end year for data input and simulations.
        self.tvec_dt = 1.0/4         # Default timestep for simulations.
        
        self.recursion_limit = 100      # Limit for recursive references, primarily used in avoiding circular references for definitions using dependencies.
        
        self.parser = FunctionParser(debug = False)      # Decomposes and evaluates functions written as strings, in accordance with a grammar defined within the parser object.
        
        # Project-data workbook metadata. Constant regardless of cascade structure.
        self.databook = odict()
        self.databook['sheet_names'] = odict()
        self.databook['sheet_names']['pops'] = 'Population Definitions'
        self.databook['sheet_names']['transmat'] = 'Transfer Definitions'
        self.databook['sheet_names']['transval'] = 'Transfer Details'
        self.databook['sheet_names']['charac'] = 'Epidemic Characteristics'
        self.databook['sheet_names']['linkpars'] = 'Cascade Parameters'
        
        self.loadCustomDatabookFramework()      # NOTE: The cascade workbook will hopefully one day be capable of replicating this, making it unnecessary.
        self.plot_settings = PlottingSettings()        
        
        self.loadCascadeSettings(cascade_path)

    
    def startFresh(self):
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
        ''' Resets, then generates node and link settings based on cascade spreadsheet. '''
        
        self.startFresh()
        
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
                    if sheet_label in self.databook['sheet_names'].keys():
                        raise OptimaException('ERROR: Custom sheet label "%s" is already used as a standard label when generating project databooks.' % sheet_label)
                    self.databook['custom_sheet_names'][sheet_label] = sheet_name        
        
        
        #%% First core sheet: Compartments
        # This worksheet defines nodes/compartments of the cascade network, with label and name required at a bare minimum.
        
        cid_label = None
        cid_name = None
        cid_coords = None
        cid_no_plot = None
        cid_dead = None
        cid_junction = None
        for col_id in xrange(ws_nodes.ncols):
            if ws_nodes.cell_value(0, col_id) == 'Code Label': cid_label = col_id
            if ws_nodes.cell_value(0, col_id) == 'Full Name': cid_name = col_id
            if ws_nodes.cell_value(0, col_id) == 'Plot Coordinates': cid_coords = col_id
            if ws_nodes.cell_value(0, col_id) == 'No Plot': cid_no_plot = col_id
            if ws_nodes.cell_value(0, col_id) == 'Dead Tag': cid_dead = col_id
            if ws_nodes.cell_value(0, col_id) == 'Junction': cid_junction = col_id
        if None in [cid_label, cid_name]:
            raise OptimaException('ERROR: Cascade compartment worksheet does not have correct column headers.')
        
        # Append node labels and full names to relevant odicts and lists. Labels are crucial.
        for row_id in xrange(ws_nodes.nrows):
            if row_id > 0 and ws_nodes.cell_value(row_id, cid_label) not in ['']:
                node_label = str(ws_nodes.cell_value(row_id, cid_label))
                self.node_specs[node_label] = odict()
                self.node_names.append(str(ws_nodes.cell_value(row_id, cid_name)))
                
                # Store optional graph coordinates for this node in the cascade network. Used when plotting schematic.
                if not cid_coords is None:
                    val = str(ws_nodes.cell_value(row_id, cid_coords))
                    if val not in ['']:
                        try: coords = val.strip('()').split(',')
                        except: raise OptimaException('ERROR: Plot Coordinates for "%s" are not written in "(x,y)" format.' % node_label)
                        self.node_specs[node_label]['coords'] = (float(coords[0]),float(coords[1]))
                
                # Store optional information about whether to avoid plotting node in the schematic.
                if not cid_no_plot is None:
                    val = str(ws_nodes.cell_value(row_id, cid_no_plot))
                    if val not in ['']:
                        self.node_specs[node_label]['tag_no_plot'] = val
                
                # Store optional information about whether this node is a compartment for the dead.
                if not cid_dead is None:
                    val = str(ws_nodes.cell_value(row_id, cid_dead))
                    if val not in ['']:
                        self.node_specs[node_label]['tag_dead'] = val

                # Optionally note down whether the compartment is just a junction.
                # Junctions empty themselves via outflows at the end of each model timestep.
                if not cid_junction is None:
                    val = str(ws_nodes.cell_value(row_id, cid_junction))
                    if val not in ['']:
                        self.node_specs[node_label]['junction'] = val
                        self.junction_labels.append(node_label)


        #%% Second core sheet: Transitions
        # This worksheet describes links/transitions of the cascade network with respect to the nodes/compartments.

        # Validate bijectivity of compartments defined in previous sheet and those that mark the edges of the adjacency matrix in this sheet.
        test = {}
        for row_id in xrange(ws_links.nrows):
            val = str(ws_links.cell_value(row_id, 0))
            if row_id > 0 and val not in ['']:
                if val not in self.node_specs.keys():
                    raise OptimaException('ERROR: Cascade transitions worksheet has a row header ("%s") that is not a known compartment code label.' % val)
                test[val] = True
        for node_label in self.node_specs.keys():
            if node_label not in test.keys():
                raise OptimaException('ERROR: Compartment code label "%s" is not represented in row headers of transitions worksheet.' % node_label)
        test = {}
        for col_id in xrange(ws_links.ncols):
            val = str(ws_links.cell_value(0, col_id))
            if col_id > 0 and val not in ['']:
                if val not in self.node_specs.keys():
                    raise OptimaException('ERROR: Cascade transitions worksheet has a column header ("%s") that is not a known compartment code label.' % val)
                test[val] = True
        for node_label in self.node_specs.keys():
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
                        if not val in self.links.keys():
                            self.links[val] = []
                        self.links[val].append((n1, n2))
                        
                    
        #%% Third core sheet: Characteristics
        # This worksheet describes characteristics of the model (i.e. they can be considered as outputs).
        # A label, name and set of included compartments is required.
                    
        # Sweep through column headers to make sure the right tags exist. Basically checking spreadsheet format.
        cid_label = None
        cid_name = None
        cid_sheet = None
        cid_percentage = None
        cid_denom = None
        cid_order = None
        cid_default = None
        cid_entry = None
        cid_include_start = None
        cid_include_end = None
        for col_id in xrange(ws_characs.ncols):
            if ws_characs.cell_value(0, col_id) == 'Code Label': cid_label = col_id
            if ws_characs.cell_value(0, col_id) == 'Full Name': cid_name = col_id
            if ws_characs.cell_value(0, col_id) == 'Sheet Label': cid_sheet = col_id
            if ws_characs.cell_value(0, col_id) == 'Plot Percentage': cid_percentage = col_id
            if ws_characs.cell_value(0, col_id) == 'Denominator': cid_denom = col_id
            if ws_characs.cell_value(0, col_id) == 'Databook Order': cid_order = col_id
            if ws_characs.cell_value(0, col_id) == 'Default Value': cid_default = col_id
            if ws_characs.cell_value(0, col_id) == 'Entry Point': cid_entry = col_id
            if ws_characs.cell_value(0, col_id) == 'Includes': cid_include_start = col_id
        
        # Work out where the 'include' columns end when defining cascade characteristics.
        cid_list = np.array(sorted([cid_label, cid_name, cid_sheet, cid_percentage, cid_denom, cid_order, cid_default, cid_entry, ws_characs.ncols]))
        cid_include_end = cid_list[sum(cid_include_start > cid_list)] - 1
        
        if None in [cid_label, cid_name, cid_include_start, cid_include_end]:
            raise OptimaException('ERROR: Cascade characteristics worksheet does not have correct column headers.')
        
        # Actually append characteristic labels and full names to relevant odicts and lists. Labels are crucial.
        entry_dict = {}
        standard_sheet_count = ws_characs.nrows - 1      # All characteristics are printed to the standard databook sheet to begin with.
        for row_id in xrange(ws_characs.nrows):
            if row_id > 0 and ws_characs.cell_value(row_id, cid_label) not in ['']:
                charac_label = str(ws_characs.cell_value(row_id, cid_label))
                charac_name = str(ws_characs.cell_value(row_id, cid_name))
                
                self.charac_name_labels[charac_name] = charac_label
                self.charac_specs[charac_label] = odict()
                self.charac_specs[charac_label]['name'] = charac_name
                self.charac_specs[charac_label]['includes'] = []
                self.charac_specs[charac_label]['sheet_label'] = 'charac'    # Default label of databook sheet to print to.
                
                # Attribute a custom sheet label to this characteristic if available.
                if not cid_sheet is None:
                    val = str(ws_characs.cell_value(row_id, cid_sheet))
                    if val not in ['']:
                        if not val in self.databook['custom_sheet_names'].keys():
                            raise OptimaException('ERROR: Project databook sheet-label "%s" for characteristic "%s" has not been previously defined.' % (val, charac_label))
                        self.charac_specs[charac_label]['sheet_label'] = val
                        standard_sheet_count -= 1
                
                # Work out which compartments/characteristics this characteristic is a sum of.
                for col_id_inc in xrange(cid_include_end - cid_include_start + 1):
                    col_id = cid_include_start + col_id_inc
                    val = str(ws_characs.cell_value(row_id, col_id))
                    if val not in ['']:
                        if val not in self.node_specs.keys() + self.charac_specs.keys()[:-1]:
                            raise OptimaException('ERROR: Cascade characteristic "%s" is being defined with an inclusion reference to "%s", which has not been defined yet.' % (charac_label, val))
                        self.charac_specs[charac_label]['includes'].append(val)
                
                flat_list, dep_list = flattenDict(input_dict = self.charac_specs, base_key = charac_label, sub_keys = ['includes'], limit = self.recursion_limit)
                if len(flat_list) != len(set(flat_list)):
                    raise OptimaException('ERROR: Cascade characteristic "%s" contains duplicate references to a compartment (in its numerator) when all the recursion is flattened out.' % (charac_label))
                numerator_ref_list = dcp(flat_list)   # Useful for checking what compartments the numerator characteristic references.
                
                # Work out which compartment/characteristic this characteristic is normalised by.
                if not cid_denom is None:
                    val = str(ws_characs.cell_value(row_id, cid_denom))
                    if val not in ['']:
                        if val not in self.node_specs.keys() + self.charac_specs.keys()[:-1]:
                            raise OptimaException('ERROR: Cascade characteristic "%s" is being defined with reference to denominator "%s", which has not been defined yet.' % (charac_label, val))
                        self.charac_specs[charac_label]['denom'] = val
                
                # Store whether characteristic should be converted to percentages when plotting.
                if not cid_percentage is None:
                    val = str(ws_characs.cell_value(row_id, cid_percentage))
                    if val not in ['']:
                        self.charac_specs[charac_label]['plot_percentage'] = val
                
                # Store order that characteristics should be printed in project databook.
                # Defaults to high value (i.e. low priority), even if column does not exist.
                self.charac_specs[charac_label]['databook_order'] = ws_characs.nrows + 1   # Any value missing from a databook order column means their printouts are lowest priority.
                if not cid_order is None:
                    val = ws_characs.cell_value(row_id, cid_order)
                    if val not in ['']:
                        self.charac_specs[charac_label]['databook_order'] = int(val)
                        if int(val) < 0:
                            standard_sheet_count -= 1   # Unprinted characteristics do not end up on the standard sheet.
                        
                # Store characteristic default value if available.
                if not cid_default is None:
                    val = str(ws_characs.cell_value(row_id, cid_default))
                    if val not in ['']:
                        self.charac_specs[charac_label]['default'] = float(val)
                        
                # Store characteristic entry point if available.
                # Without this, the characteristic will not be converted into an initial value for the model.
                if not cid_entry is None:
                    val = str(ws_characs.cell_value(row_id, cid_entry))
                    if val not in ['']:
                        self.charac_specs[charac_label]['entry_point'] = val
                        
                        if self.charac_specs[charac_label]['databook_order'] < 0:
                            raise OptimaException('ERROR: Characteristic "%s" cannot be used to seed a model (i.e. have an entry point) if it does not feature in the databook (i.e. if it has a negative databook order).' % (charac_label))
                        if 'denom' in self.charac_specs[charac_label].keys():
                            denom_label = self.charac_specs[charac_label]['denom']
                            if denom_label not in self.charac_specs.keys()[:-1] or 'entry_point' not in self.charac_specs[denom_label].keys():
                                raise OptimaException('ERROR: At this time, characteristic "%s" cannot be used to seed a model (i.e. have an entry point) if its denominator "%s" is not a model-seeding characteristic.' % (charac_label, denom_label))
                        if val in entry_dict.keys():
                            raise OptimaException('ERROR: Cascade characteristic "%s" is not the first to use "%s" as an entry point.' % (charac_label, val))
                        if val not in self.node_specs.keys():
                            raise OptimaException('ERROR: There is no compartment named "%s" that can be used as an entry point by cascade characteristic "%s".' % (val, charac_label))                        
                        numerator_ref_list = set(numerator_ref_list)
                        try: numerator_ref_list.remove(val)
                        except: raise OptimaException('ERROR: Entry point "%s" for characteristic "%s" must be a compartment it includes (in its numerator), either directly or via characteristic reference.' % (val, charac_label))
                        entry_dict[val] = dcp(list(numerator_ref_list))
        
        # If all characteristics in this sheet are to be printed to custom databook sheets, no need for a default.
        if standard_sheet_count <= 0:
            self.make_sheet_characs = False
                            
        # Ensuring that no two characteristics with entry-points include each other's entry-points (or more complex referencing cycles).
        # This allows for disaggregated characteristic values to be partitioned from aggregate characteristics during model-seeding.
        for entry_label in entry_dict.keys():
            try:
                flattenDict(input_dict = entry_dict, base_key = entry_label, limit = self.recursion_limit)
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
        for col_id in xrange(ws_pars.ncols):
            if ws_pars.cell_value(0, col_id) == 'Transition Tag': cid_tag = col_id
            if ws_pars.cell_value(0, col_id) == 'Code Label': cid_label = col_id
            if ws_pars.cell_value(0, col_id) == 'Full Name': cid_name = col_id
            if ws_pars.cell_value(0, col_id) == 'Sheet Label': cid_sheet = col_id
            if ws_pars.cell_value(0, col_id) == 'Default Value': cid_default = col_id
            if ws_pars.cell_value(0, col_id) == 'Databook Order': cid_order = col_id
            if ws_pars.cell_value(0, col_id) == 'Function': cid_function = col_id
        if None in [cid_tag, cid_label, cid_name]:
            raise OptimaException('ERROR: Cascade transition-parameters worksheet does not have correct column headers.')
        
        # Store transition details in settings and make sure there is tag bijectivity between this sheet and the transition matrix.
        standard_sheet_count = ws_pars.nrows - 1      # All parameters are printed to the standard databook sheet to begin with.
        for row_id in xrange(ws_pars.nrows):
            tag = str(ws_pars.cell_value(row_id, cid_tag))
            label = str(ws_pars.cell_value(row_id, cid_label))
            name = str(ws_pars.cell_value(row_id, cid_name))
            if row_id > 0 and label not in ['']:
                self.linkpar_specs[label] = {'name':name, 'sheet_label':'linkpars'}
                self.linkpar_name_labels[name] = label
                if tag not in ['']:
                    if tag not in self.links:
                        raise OptimaException('ERROR: Cascade transition-parameter worksheet has a tag (%s) that is not in the transition matrix.' % tag)
                    self.linkpar_specs[label]['tag'] = tag
                else:
                    self.par_deps[label] = True     # Untagged parameters are dependencies for actual tagged transitions.
                
                # Attribute a custom sheet label to this parameter if available.
                if not cid_sheet is None:
                    val = str(ws_pars.cell_value(row_id, cid_sheet))
                    if val not in ['']:
                        if not val in self.databook['custom_sheet_names'].keys():
                            raise OptimaException('ERROR: Project databook sheet-label "%s" for parameter "%s" has not been previously defined.' % (val, label))
                        self.linkpar_specs[label]['sheet_label'] = val
                        standard_sheet_count -= 1
                
                # Store order that cascade parameters should be printed in project databook.
                self.linkpar_specs[label]['databook_order'] = ws_pars.nrows+1   # Any value missing from a databook order column means their printouts are lowest priority.
                if not cid_order is None:
                    val = ws_pars.cell_value(row_id, cid_order)
                    if val not in ['']:
                        self.linkpar_specs[label]['databook_order'] = int(val)
                        if int(val) < 0:    
                            standard_sheet_count -= 1   # Unprinted parameters do not end up on the standard sheet.
                
                # Store parameter default value if available.
                if not cid_default is None:
                    val = str(ws_pars.cell_value(row_id, cid_default))
                    if val not in ['']:
                        self.linkpar_specs[label]['default'] = float(val)
                        
                # Store the token stack corresponding to a custom function that defines this parameter, if available.
                if not cid_function is None:
                    val = str(ws_pars.cell_value(row_id, cid_function))
                    if val not in ['']:
                        if 'default' in self.linkpar_specs[label]:
                            raise OptimaException('ERROR: Parameter "%s" is a custom function of other parameters and characteristics. Specifying a default is thus restricted, so as to avoid user confusion.' % label)
                        if self.linkpar_specs[label]['databook_order'] >= 0:
                            raise OptimaException('ERROR: Parameter "%s" is a custom function of other parameters and characteristics. Tag "Databook Order" column with a negative number so that conflicts with user-provided values do not arise.' % label)
                        self.par_funcs[label] = True
                        expr_stack, var_dict = self.parser.produceStack(val)
                        self.linkpar_specs[label]['f_stack'] = expr_stack
                        self.linkpar_specs[label]['deps'] = var_dict
                        for var in var_dict.keys():
                            if not var in self.charac_specs.keys():
                                if not var in self.linkpar_specs.keys():
                                    raise OptimaException('ERROR: Dependency "%s" has not been defined by the time "%s" is loaded into settings.' % (var, label))
                            else:
                                if not var in self.charac_deps.keys():
                                    flat_list, dep_list = flattenDict(input_dict = self.charac_specs, base_key = var, sub_keys = ['includes','denom'], limit = self.recursion_limit)
                                    for dep in dep_list:
                                        self.charac_specs[dep]['par_dependency'] = True
                                        self.charac_deps[dep] = True
        
        # If all parameters in this sheet are to be printed to custom databook sheets, no need for a default.
        if standard_sheet_count <= 0:
            self.make_sheet_linkpars = False
        
        # Final validations.            
        for tag in self.links.keys():
            if tag not in [x['tag'] for x in self.linkpar_specs[:] if 'tag' in x]:
                raise OptimaException('ERROR: Transition matrix tag "%s" is not represented in transition-parameter worksheet.' % tag)
        
        test_labels = self.linkpar_specs.keys() + self.charac_specs.keys()       
        if len(test_labels) != len(set(test_labels)):
            raise OptimaException('ERROR: Cascade workbook appears to have duplicate characteristic/parameter code labels.')
        
        test_names = self.linkpar_name_labels.keys() + self.charac_name_labels.keys()
        if len(test_names) != len(set(test_names)):
            raise OptimaException('ERROR: Cascade workbook appears to have duplicate characteristic/parameter full names.')
            
            
    def plotCascade(self):
        
        import networkx as nx
        
        fig, ax = pl.subplots(figsize=(10,10))
        G = nx.DiGraph()
        plottable_nodes = [nid for nid in self.node_specs.keys() if 'tag_no_plot' not in self.node_specs[nid]]
        plottable_links = [link for lid in self.links for link in self.links[lid] if (link[0] in plottable_nodes and link[1] in plottable_nodes)]
        G.add_nodes_from(plottable_nodes)
        G.add_edges_from(plottable_links)

        # Use plot coordinates if stored and arrange the rest of the cascade out in a unit circle.
        pos = {}
        num_nodes = len(plottable_nodes)
        k = 0
        for node in plottable_nodes:
            try: pos[node] = (self.node_specs[node]['coords'][0], self.node_specs[node]['coords'][1])
            except: pos[node] = (np.sin(2.0*np.pi*k/num_nodes), np.cos(2.0*np.pi*k/num_nodes))
            k += 1
        
        # Generate edge label dictionary with tags from spreadsheet.
        el = {}
        for par_name in self.linkpar_specs.keys():
            if 'tag' in self.linkpar_specs[par_name]:
                for link in self.links[self.linkpar_specs[par_name]['tag']]:
                    el[link] = self.linkpar_specs[par_name]['tag']

        nx.draw_networkx_nodes(G, pos, node_shape = 'o', nodelist = [x for x in G.nodes() if not 'junction' in self.node_specs[x].keys()], node_size = 1250, node_color = 'w')
        nx.draw_networkx_nodes(G, pos, node_shape = 's', nodelist = [x for x in G.nodes() if 'junction' in self.node_specs[x].keys()], node_size = 750, node_color = 'w')
        ax.axis('tight')
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
#        nx.draw_networkx_edge_labels(G, pos, edge_labels = el, label_pos = 0.25, font_size = 14)
        
        [sp.set_visible(False) for sp in ax.spines.values()]
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Cascade Schematic')
        pl.show()


    def loadCustomDatabookFramework(self):
        
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
#         self.countrybook['sheet_names']['unitcost'] = 'Unit costs'
#         self.countrybook['sheet_names']['poptransitions'] = 'Population transitions'

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
        self.countrybook['disaggregations'] = odict() # ['smears','strains']
        # other values
        self.countrybook['disaggregations']['strains'] = ['DS-TB','MDR-TB','XDR-TB']
        self.countrybook['disaggregations']['smears'] = ['Smear-','Smear+'] # also potentially 'Extrapulmonary'
        self.countrybook['disaggregations']['populations'] = [] # determined dynamically at runtime
        self.countrybook['disaggregations']['regimens'] = ['DS-TB Regimen','MDR-TB Regimen','XDR-TB Regimen']
        self.countrybook['disaggregations']['programs'] = [] # determined dynamically at runtime
        
        # for univalue sheets, includes information on how data should be disaggregated 
        self.countrybook['sheet_classes'] = odict()
        self.countrybook['sheet_classes']['univalue'] = odict()
        self.countrybook['sheet_classes']['univalue']['population_sizes'] = ['populations']
        self.countrybook['sheet_classes']['univalue']['total_cases'] = ['populations','strains']
        self.countrybook['sheet_classes']['univalue']['incident_cases'] = ['populations','smears','strains']
        
        # sheet specific values
        self.countrybook['sheet_values'] = odict()
        self.countrybook['sheet_values']['other_epidemiology'] = odict()
        self.countrybook['sheet_values']['other_epidemiology']['Percentage of people vaccinated per year'] = ['populations']
        self.countrybook['sheet_values']['other_epidemiology']['Percentage of people who die from non-TB-related causes per year'] = ['populations']
        self.countrybook['sheet_values']['other_epidemiology']['Percentage of people who die from TB-related deaths per year'] = ['populations','smears','strains']
        self.countrybook['sheet_values']['other_epidemiology']['Birth rate (births per woman per year)'] = ['populations']
        
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


class PlottingSettings():
    
    
    def __init__(self):
        print("Loading plotting settings")
        self.defaultSettings()
        self.devSettings()
        
        
    def defaultSettings(self):
        
        pl.rcParams['font.size'] = 12
        pl.rcParams['font.family'] = 'sans-serif'
        
        pl.rcParams['savefig.dpi'] = 300
        
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