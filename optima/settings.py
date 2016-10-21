#%% Imports

from utils import odict, OptimaException

import xlrd
import pylab as pl
import numpy as np
 
    
#%% Settings class (for data that is effectively static per epidemic context)
#
#   This refers specifically to cascade metadata (loaded from a cascade workbook) and general defaults.
#   A cascade structure comprises of compartments/nodes and transitions/links that detail the flow of people through stages of a disease.
#   In general use, the attributes of this object are static regardless of project/simulation.
#
#   Notes:
#   Settings MUST be loaded into a project beforehand, as it details the framework for dealing with data and running the model.

class Settings(object):
    def __init__(self, cascade_path = './cascade.xlsx'):
        self.tvec_start = 2000.0     # Default start year for data input and simulations.
        self.tvec_end = 2030.0       # Default end year for data input and simulations.
        self.tvec_dt = 1.0/4          # Default timestep for simulations.
        
        self.plotSettings = PlottingSettings()
        
        
        self.startFresh()       # NOTE: Unnecessary as loading a cascade calls this anyway. But left here to be explicit. 
        self.loadCascadeSettings(cascade_path)
        
        # Project-data workbook metadata. Constant regardless of cascade structure.
        self.databook = odict()
        self.databook['sheet_names'] = odict()
        self.databook['sheet_names']['pops'] = 'Population Definitions'
        self.databook['sheet_names']['transmat'] = 'Transfer Definitions'
        self.databook['sheet_names']['transval'] = 'Transfer Details'
        self.databook['sheet_names']['charac'] = 'Epidemic Characteristics'
        self.databook['sheet_names']['linkpars'] = 'Cascade Parameters'
    
    def startFresh(self):
        ''' Resets all cascade contents and settings that are fundamental to how a project is structured. '''
            
        self.node_specs = odict()               # Relates compartment code-labels with special tags, if they exist. 
                                                # Key is a node label. Value is a dict including a 'dead' tag and networkx-related information.
        self.node_names = []                    # A corresponding list of full names for compartments.
        self.charac_specs = odict()             # Relates code-labels for defined characteristics (e.g. prevalence) with labels of compartments used in their definition.
                                                # Key is a characteristic label. Value is a dict containing characteristic name, a list of 'inclusions' and a normalising characteristic or compartment.
                                                # Be aware that inclusions/normalisations may refer to characteristics in the same odict.
        self.charac_name_labels = odict()       # Key is a characteristic name. Value is a characteristic label. (A partial reversed charac_specs.)
        self.links = odict()                    # Key is a tag. Value is a list of compartment-label tuple.
        self.linkpar_specs = odict()            # Key is a link-parameter label. Value is a dict including link tag, link-parameter name, default value.
        self.linkpar_name_labels = odict()      # Key is a link-parameter name. Value is a link-parameter label. (A partial reversed linkpar-specs.)
    
    def loadCascadeSettings(self, cascade_path):
        ''' Resets, then generates node and link settings based on cascade spreadsheet. '''
        
        self.startFresh()
        
        try: workbook = xlrd.open_workbook(cascade_path)
        except: raise OptimaException('ERROR: Cannot find cascade workbook from which to load model structure.')
        ws_nodes = workbook.sheet_by_name('Compartments')
        ws_charac = workbook.sheet_by_name('Cascade Characteristics')
        ws_links = workbook.sheet_by_name('Transitions')
        ws_pars = workbook.sheet_by_name('Transition Parameters')
        
        # First sheet: Compartments
        # Sweep through column headers to make sure the right tags exist. Basically checking spreadsheet format.
        cid_label = None
        cid_name = None
        cid_coords = None
        cid_no_plot = None
        cid_dead = None
        for col_id in xrange(ws_nodes.ncols):
            if ws_nodes.cell_value(0, col_id) == 'Code Label': cid_label = col_id
            if ws_nodes.cell_value(0, col_id) == 'Full Name': cid_name = col_id
            if ws_nodes.cell_value(0, col_id) == 'Plot Coordinates': cid_coords = col_id
            if ws_nodes.cell_value(0, col_id) == 'No Plot': cid_no_plot = col_id
            if ws_nodes.cell_value(0, col_id) == 'Dead Tag': cid_dead = col_id
        if None in [cid_label, cid_name]:
            raise OptimaException('ERROR: Cascade compartment worksheet does not have correct column headers.')
        
        # Actually append node labels and full names to relevant odicts and lists. Labels are crucial.
        for row_id in xrange(ws_nodes.nrows):
            if row_id > 0 and ws_nodes.cell_value(row_id, cid_label) not in ['']:
                node_label = str(ws_nodes.cell_value(row_id, cid_label))
                self.node_specs[node_label] = odict()
                self.node_names.append(str(ws_nodes.cell_value(row_id, cid_name)))
                
                # Not crucial, but store node plotting coordinates if available.
                try:
                    in_coords = str(ws_nodes.cell_value(row_id, cid_coords))
                    coords = in_coords.strip('()').split(',')
                    self.node_specs[node_label]['coords'] = (float(coords[0]),float(coords[1]))
                except: pass
                
                # Not crucial, but store information whether to plot if available.
                if not cid_no_plot is None:
                    val = ws_nodes.cell_value(row_id, cid_no_plot)
                    if val not in ['']:
                        self.node_specs[node_label]['tag_no_plot'] = val
                
                # Ditto with tags for dead compartments.
                if not cid_dead is None:
                    val = ws_nodes.cell_value(row_id, cid_dead)
                    if val not in ['']:
                        self.node_specs[node_label]['tag_dead'] = val
                        
        # Second sheet: Cascade Characteristics
        # Sweep through column headers to make sure the right tags exist. Basically checking spreadsheet format.
        cid_label = None
        cid_name = None
        cid_percentage = None
        cid_denom = None
        cid_include_start = None
        cid_include_end = None
        cid_order = None
        cid_default = None
        for col_id in xrange(ws_charac.ncols):
            if ws_charac.cell_value(0, col_id) == 'Code Label': cid_label = col_id
            if ws_charac.cell_value(0, col_id) == 'Full Name': cid_name = col_id
            if ws_charac.cell_value(0, col_id) == 'Plot Percentage': cid_percentage = col_id
            if ws_charac.cell_value(0, col_id) == 'Denominator': cid_denom = col_id
            if ws_charac.cell_value(0, col_id) == 'Includes': cid_include_start = col_id
            if ws_charac.cell_value(0, col_id) == 'Databook Order': cid_order = col_id
            if ws_charac.cell_value(0, col_id) == 'Default Value': cid_default = col_id
        
        # Work out where the 'include' columns end when defining cascade characteristics.
        cid_list = np.array(sorted([cid_label, cid_name, cid_percentage, cid_denom, cid_order, cid_default, ws_charac.ncols]))
        cid_include_end = cid_list[sum(cid_include_start > cid_list)] - 1
        
        if None in [cid_label, cid_name, cid_include_start, cid_include_end]:
            raise OptimaException('ERROR: Cascade characteristics worksheet does not have correct column headers.')
        
        # Actually append characteristic labels and full names to relevant odicts and lists. Labels are crucial.
        for row_id in xrange(ws_charac.nrows):
            if row_id > 0 and ws_charac.cell_value(row_id, cid_label) not in ['']:
                charac_label = str(ws_charac.cell_value(row_id, cid_label))
                charac_name = str(ws_charac.cell_value(row_id, cid_name))
                
                self.charac_name_labels[charac_name] = charac_label
                self.charac_specs[charac_label] = odict()
                self.charac_specs[charac_label]['name'] = charac_name
                self.charac_specs[charac_label]['includes'] = []
                
                for col_id_inc in xrange(cid_include_end - cid_include_start + 1):
                    col_id = cid_include_start + col_id_inc
                    val = str(ws_charac.cell_value(row_id, col_id))
                    if val not in ['']:
                        if val not in self.node_specs.keys() + self.charac_specs.keys()[:-1]:
                            raise OptimaException('ERROR: Cascade characteristic %s is being defined with an inclusion reference to %s, which has not been defined yet.' % (charac_label, val))
                        self.charac_specs[charac_label]['includes'].append(val)
        
                if not cid_denom is None:
                    val = str(ws_charac.cell_value(row_id, cid_denom))
                    if val not in ['']:
                        if val not in self.node_specs.keys() + self.charac_specs.keys()[:-1]:
                            raise OptimaException('ERROR: Cascade characteristic %s is being defined with reference to denominator %s, which has not been defined yet.' % (charac_label, val))
                        self.charac_specs[charac_label]['denom'] = val
                
                # Store whether characteristic should be converted to percentages when plotting.
                if not cid_percentage is None:
                    val = str(ws_charac.cell_value(row_id, cid_percentage))
                    if val not in ['']:
                        self.charac_specs[charac_label]['plot_percentage'] = val
                
                # Store order that characteristics should be printed in project databook.
                self.charac_specs[charac_label]['databook_order'] = ws_charac.nrows+1   # Any value missing from a databook order column means their printouts are lowest priority.
                if not cid_order is None:
                    val = ws_charac.cell_value(row_id, cid_order)
                    if val not in ['']:
                        self.charac_specs[charac_label]['databook_order'] = int(val)
                        
                # Store characteristic default value if available.
                if not cid_default is None:
                    def_val = str(ws_charac.cell_value(row_id, cid_default))
                    if def_val not in ['']:
                        self.charac_specs[charac_label]['default'] = float(def_val)
                    
        
        # Third sheet: Transitions
        # Quality-assurance test for the spreadsheet format.
        test = []
        for row_id in xrange(ws_links.nrows):
            val = ws_links.cell_value(row_id, 0)
            if row_id > 0 and val not in ['']:
                if val not in self.node_specs.keys():
                    raise OptimaException('ERROR: Cascade transitions worksheet has a row header (%s) that is not a known compartment code label.' % val)
                test.append(val)
        for label in self.node_specs.keys():
            if label not in test:
                raise OptimaException('ERROR: Compartment code label (%s) is not represented in row headers of transitions worksheet.' % label)
        test = []
        for col_id in xrange(ws_links.ncols):
            val = ws_links.cell_value(0, col_id)
            if col_id > 0 and val not in ['']:
                if val not in self.node_specs.keys():
                    raise OptimaException('ERROR: Cascade transitions worksheet has a column header (%s) that is not a known compartment code label.' % val)
                test.append(val)
        for label in self.node_specs.keys():
            if label not in test:
                raise OptimaException('ERROR: Compartment code label (%s) is not represented in column headers of transitions worksheet.' % label)
        
        # Store linked compartment tuples by tag.
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
            
        # Fourth sheet: Transition Parameters
        # Sweep through column headers to make sure the right tags exist. Basically checking spreadsheet format.
        cid_tag = None
        cid_label = None
        cid_name = None
        cid_default = None
        cid_order = None
        for col_id in xrange(ws_pars.ncols):
            if ws_pars.cell_value(0, col_id) == 'Tag': cid_tag = col_id
            if ws_pars.cell_value(0, col_id) == 'Code Label': cid_label = col_id
            if ws_pars.cell_value(0, col_id) == 'Full Name': cid_name = col_id
            if ws_pars.cell_value(0, col_id) == 'Default Value': cid_default = col_id
            if ws_pars.cell_value(0, col_id) == 'Databook Order': cid_order = col_id
        if None in [cid_tag, cid_label, cid_name]:
            raise OptimaException('ERROR: Cascade transition-parameters worksheet does not have correct column headers.')
        
        # Store transition details in settings and make sure there is tag bijectivity between this sheet and the transition matrix.
        for row_id in xrange(ws_pars.nrows):
            tag = str(ws_pars.cell_value(row_id, cid_tag))
            label = str(ws_pars.cell_value(row_id, cid_label))
            name = str(ws_pars.cell_value(row_id, cid_name))
            if row_id > 0 and tag not in [''] and label not in ['']:
                if tag not in self.links:
                    raise OptimaException('ERROR: Cascade transition-parameter worksheet has a tag (%s) that is not in the transition matrix.' % tag)
                self.linkpar_specs[label] = {'tag':tag, 'name':name}
                self.linkpar_name_labels[name] = label
                
                # Store order that cascade parameters should be printed in project databook.
                self.linkpar_specs[label]['databook_order'] = ws_pars.nrows+1   # Any value missing from a databook order column means their printouts are lowest priority.
                if not cid_order is None:
                    val = ws_pars.cell_value(row_id, cid_order)
                    if val not in ['']:
                        self.linkpar_specs[label]['databook_order'] = int(val)                
                
                # Store parameter default value if available.
                if not cid_default is None:
                    def_val = str(ws_pars.cell_value(row_id, cid_default))
                    if def_val not in ['']:
                        self.linkpar_specs[label]['default'] = float(def_val)
                    
        for tag in self.links.keys():
            if tag not in [x['tag'] for x in self.linkpar_specs[:]]:
                raise OptimaException('ERROR: Transition matrix tag (%s) is not represented in transition-parameter worksheet.' % tag)
        if len(self.linkpar_specs.keys()) != len(set(self.linkpar_specs.keys())):
            raise OptimaException('ERROR: Cascade transition-parameter worksheet appears to have duplicate parameter code labels.')
    




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
    
    