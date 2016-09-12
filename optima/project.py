#%% Imports

from utils import odict
from model import model
from uuid import uuid4 as uuid



#%% Project class (i.e. one self-contained geographical unit)

class Project(object):
    ''' The main Optima project class. Almost all Optima functionality is provided by this class. '''

    def __init__(self, name='default'):
        ''' Initialize project. '''

        self.name = name
        self.uid = uuid()
        self.data = {}

        self.parsets = odict()
        self.results = odict()

        return None
        
        
    def runSim(self):
        ''' Run model using a selected parset and store/return results. '''

        results = model()        
        
        return results