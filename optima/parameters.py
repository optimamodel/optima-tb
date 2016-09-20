#%% Imports

from utils import odict



#%% Parameter class that stores one array of values converted from raw project data

class Parameter(object):
    ''' Class to hold one group of parameter values. '''
    
    def __init__(self, label, t = None, y = None):
        self.label = label
        if t is None: t = odict()
        if y is None: y = odict()
        self.t = t      # Time data.
        self.y = y      # Value data.


#%% Parset class that contains one set of parameters converted from raw project data

class ParameterSet(object):
    ''' Class to hold all parameters. '''
    
    def __init__(self, name='default'):
        self.name = name 
        self.pop_names = []         # List of population names.
        self.pop_labels = []        # List of population labels.
        self.pars = []
        self.par_ids = {}
    
    def makePars(self, data):
        self.pop_names = data['pops']['name_labels'].keys()
        for name in self.pop_names:
            self.pop_labels.append(data['pops']['name_labels'][name])
        for l, label in enumerate(data['linkpars']):
            self.par_ids[label] = l
            self.pars.append(Parameter(label))
            for pop_id in data['linkpars'][label]:
                self.pars[-1].t[pop_id] = data['linkpars'][label][pop_id]['t']
                self.pars[-1].y[pop_id] = data['linkpars'][label][pop_id]['y']
        