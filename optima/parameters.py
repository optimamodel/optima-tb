#%% Imports

from utils import odict, OptimaException
import numpy as np



#%% Parameter class that stores one array of values converted from raw project data

class Parameter(object):
    ''' Class to hold one group of parameter values. '''
    
    def __init__(self, label, t = None, y = None):
        self.label = label
        if t is None: t = odict()
        if y is None: y = odict()
        self.t = t      # Time data.
        self.y = y      # Value data.
        
    def interpolate(self, tvec = None, pop_label = None):
        ''' Take parameter values and construct an array matching input time vector. '''
        
        # Validate input.
        if pop_label is None: raise OptimaException('ERROR: Cannot interpolate parameter %s without knowing which population to do it for.' % pop_label)
        if tvec is None: raise OptimaException('ERROR: Cannot interpolate parameter %s without providing a time vector.' % self.label)
        
#        output = np.ones()
            
        
        if self.t[pop_label] is None:
            pass
        
        return output


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
        