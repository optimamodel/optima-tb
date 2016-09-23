#%% Imports

from utils import odict, OptimaException

import xlrd
import networkx as nx
import pylab as pl
import numpy as np



#%% Settings class (for data that is effectively static per epidemic context)

class Settings(object):
    def __init__(self, cascade_path = './cascade.xlsx'):
        self.tvec_start = 2000.0     # Default start year for data input and simulations.
        self.tvec_end = 2030.0       # Default end year for data input and simulations.
        self.tvec_dt = 0.25          # Default timestep for simulations.
        
        self.startFresh()       # NOTE: Unnecessary as loading a cascade calls this anyway. But left here to be explicit. 
        self.loadCascadeSettings(cascade_path)
        
        # Project-data workbook metadata. Constant regardless of cascade structure.
        self.databook = odict()
        self.databook['sheet_names'] = odict()
        self.databook['sheet_names']['pops'] = 'Population Definitions'
        self.databook['sheet_names']['poptrans'] = 'Population Transitions'
        self.databook['sheet_names']['linkpars'] = 'Cascade Parameters'
    
    def startFresh(self):
        ''' Resets all cascade contents and settings that are fundamental to how a project is structured. '''
            
        self.node_labels = []
        self.node_names = []
        self.node_coords = []
        self.links = odict()                    # Key is a tag. Value is a compartment-label tuple.
        self.linkpar_specs = odict()            # Key is a link-parameter label. Value is a dict including link tag, link-parameter name, default value.
        self.linkpar_name_labels = odict()      # Key is a link-parameter name. Value is a link-parameter label. (A partial reversed linkpar-specs.)
    
    def loadCascadeSettings(self, cascade_path):
        ''' Resets, then generates node and link settings based on cascade spreadsheet. '''
        
        self.startFresh()
        
        try: workbook = xlrd.open_workbook(cascade_path)
        except: raise OptimaException('ERROR: Cannot find cascade workbook from which to load model structure.')
        ws_nodes = workbook.sheet_by_name('Compartments')
        ws_links = workbook.sheet_by_name('Transitions')
        ws_pars = workbook.sheet_by_name('Transition Parameters')
        
        # First sheet: Compartments
        # Sweep through column headers to make sure the right tags exist. Basically checking spreadsheet format.
        cid_label = None
        cid_name = None
        cid_coords = None
        for col_id in xrange(ws_nodes.ncols):
            if ws_nodes.cell_value(0, col_id) == 'Code Label': cid_label = col_id
            if ws_nodes.cell_value(0, col_id) == 'Full Name': cid_name = col_id
            if ws_nodes.cell_value(0, col_id) == 'Plot Coordinates': cid_coords = col_id
        if None in [cid_label, cid_name]:
            raise OptimaException('ERROR: Cascade compartment worksheet does not have correct column headers.')
        
        # Actually append node labels and full names to relevant lists. Labels are crucial.
        for row_id in xrange(ws_nodes.nrows):
            if row_id > 0 and ws_nodes.cell_value(row_id, cid_label) not in ['']:
                self.node_labels.append(str(ws_nodes.cell_value(row_id, cid_label)))
                self.node_names.append(str(ws_nodes.cell_value(row_id, cid_name)))
                
                # Not crucial, but store node plotting coordinates if available.
                try:
                    in_coords = str(ws_nodes.cell_value(row_id, cid_coords))
                    coords = in_coords.strip('()').split(',')
                    self.node_coords.append((float(coords[0]),float(coords[1])))
                except:
                    self.node_coords.append(None)
                    
                    
        
        # Second sheet: Transitions
        # Quality-assurance test for the spreadsheet format.
        test = []
        for row_id in xrange(ws_links.nrows):
            val = ws_links.cell_value(row_id, 0)
            if row_id > 0 and val not in ['']:
                if val not in self.node_labels:
                    raise OptimaException('ERROR: Cascade transitions worksheet has a row header (%s) that is not a known compartment code label.' % val)
                test.append(val)
        for label in self.node_labels:
            if label not in test:
                raise OptimaException('ERROR: Compartment code label (%s) is not represented in row headers of transitions worksheet.' % label)
        test = []
        for col_id in xrange(ws_links.ncols):
            val = ws_links.cell_value(0, col_id)
            if col_id > 0 and val not in ['']:
                if val not in self.node_labels:
                    raise OptimaException('ERROR: Cascade transitions worksheet has a column header (%s) that is not a known compartment code label.' % val)
                test.append(val)
        for label in self.node_labels:
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
                            self.links[val] = (n1, n2)
                        else:
                            raise OptimaException('ERROR: Cascade transition matrix appears to have duplicate labels for compartment links.')
            
        # Third sheet: Transition Parameters
        # Sweep through column headers to make sure the right tags exist. Basically checking spreadsheet format.
        cid_tag = None
        cid_label = None
        cid_name = None
        cid_default = None
        for col_id in xrange(ws_pars.ncols):
            if ws_pars.cell_value(0, col_id) == 'Tag': cid_tag = col_id
            if ws_pars.cell_value(0, col_id) == 'Code Label': cid_label = col_id
            if ws_pars.cell_value(0, col_id) == 'Full Name': cid_name = col_id
            if ws_pars.cell_value(0, col_id) == 'Default Value': cid_default = col_id
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
                
                # Not crucial, but store parameter default value if available.
                try:
                    def_val = str(ws_pars.cell_value(row_id, cid_default))
                    if def_val not in ['']:
                        self.linkpar_specs[label]['default'] = float(def_val)
                except:
                    pass
                    
        for tag in self.links.keys():
            if tag not in [x['tag'] for x in self.linkpar_specs[:]]:
                raise OptimaException('ERROR: Transition matrix tag (%s) is not represented in transition-parameter worksheet.' % tag)
        if len(self.linkpar_specs.keys()) != len(set(self.linkpar_specs.keys())):
            raise OptimaException('ERROR: Cascade transition-parameter worksheet appears to have duplicate parameter code labels.')
    
    def plotCascade(self):
        fig, ax = pl.subplots(figsize=(10,10))
        G = nx.DiGraph()
        G.add_nodes_from(self.node_labels)
        G.add_edges_from(self.links[:])

        # Use plot coordinates if stored and arrange the rest of the cascade out in a unit circle.
        pos = {}
        num_nodes = len(self.node_labels)
        k = 0
        for node in self.node_labels:
            try: pos[node] = (self.node_coords[k][0], self.node_coords[k][1])
            except: pos[node] = (np.sin(2.0*np.pi*k/num_nodes), np.cos(2.0*np.pi*k/num_nodes))
            k += 1
        
        # Generate edge label dictionary with tags from spreadsheet.
        el = {}
        for par_name in self.linkpar_specs.keys():
            el[self.links[self.linkpar_specs[par_name]['tag']]] = self.linkpar_specs[par_name]['tag']

        nx.draw_networkx_nodes(G, pos, node_size = 1250, node_color = 'w')
        ax.axis('tight')
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
#        nx.draw_networkx_edge_labels(G, pos, edge_labels = el, label_pos = 0.25, font_size = 14)
        
        [sp.set_visible(False) for sp in ax.spines.values()]
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Cascade Schematic')
        pl.show()