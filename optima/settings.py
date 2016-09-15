#%% Imports

from utils import odict, OptimaException
from xlrd import open_workbook



#%% Settings class (for data that is effectively static per epidemic context)

class Settings(object):
    def __init__(self, cascade_path = './cascade.xlsx'):
#        self.dt = 0.2           # Timestep
#        self.start = 2000.0     # Default start year
#        self.end = 2030.0       # Default end year
        
        self.node_labels = []
        self.node_names = []
        self.links = odict()
        
        self.loadCascadeSettings(cascade_path)
    
    # Resets, then generates node and link settings based on cascade spreadsheet.
    def loadCascadeSettings(self, cascade_path):
        self.node_labels = []
        self.node_names = []
        self.links = odict()
        
        try: workbook = open_workbook(cascade_path)
        except: raise OptimaException('ERROR: Cannot find cascade sheet from which to load model structure.')
        ws_nodes = workbook.sheet_by_name('Compartments')
        ws_links = workbook.sheet_by_name('Transitions')
        
        # First sheet: Compartments
        # Sweep through column headers to make sure the right tags exist. Basically checking spreadsheet format.
        cid_label = None
        cid_name = None
        for col_id in xrange(ws_nodes.ncols):
            if ws_nodes.cell_value(0, col_id) == 'Code Label': cid_label = col_id
            if ws_nodes.cell_value(0, col_id) == 'Full Name': cid_name = col_id
        if None in [cid_label, cid_name]:
            raise OptimaException('ERROR: Cascade compartment worksheet does not have correct column headers.')
        
        # Actually append node labels and full names to relevant lists. Labels are crucial.
        for row_id in xrange(ws_nodes.nrows):
            if row_id > 0 and ws_nodes.cell_value(row_id, cid_label) not in ['']:
                self.node_labels.append(str(ws_nodes.cell_value(row_id, cid_label)))
                self.node_names.append(str(ws_nodes.cell_value(row_id, cid_name)))
        
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
        
        # Store linked compartment tuples as transitions.
        test = [] 
        for row_id in xrange(ws_links.nrows):
            for col_id in xrange(ws_links.ncols):
                if row_id > 0 and col_id > 0:
                    n1 = str(ws_links.cell_value(row_id, 0))
                    n2 = str(ws_links.cell_value(0, col_id))
                    val = str(ws_links.cell_value(row_id, col_id))
                    if not '' in [n1,n2,val]:
                        self.links[(n1, n2)] = val
        
        # Make sure every transition has a unique spreadsheet label.
        if len(set(self.links[:])) != len(self.links[:]):
            raise OptimaException('ERROR: Cascade transitions worksheet appears to have duplicate labels for compartment links.')