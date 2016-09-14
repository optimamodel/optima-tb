#%% Imports

from utils import OptimaException
from xlrd import open_workbook



#%% Settings class (for data that is effectively static per epidemic context)

class Settings(object):
    def __init__(self, spreadsheet_path = './cascade.xlsx'):
#        self.dt = 0.2           # Timestep
#        self.start = 2000.0     # Default start year
#        self.end = 2030.0       # Default end year
        
        self.node_labels = []
        self.node_names = []
        
        self.loadCascadeSettings(spreadsheet_path)
    
    # Generates node and link settings based on cascade spreadsheet.
    def loadCascadeSettings(self, spreadsheet_path):
        self.node_labels = []
        
        try: workbook = open_workbook(spreadsheet_path)
        except: raise OptimaException('ERROR: Cannot find cascade sheet from which to load model structure.')
        ws_nodes = workbook.sheet_by_name('Compartments')
        
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