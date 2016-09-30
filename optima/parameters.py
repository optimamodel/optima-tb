#%% Imports

from utils import odict, OptimaException
from interpolation import interpolateFunc

from copy import deepcopy as dcp
import numpy as np



#%% Parameter class that stores one array of values converted from raw project data

class Parameter(object):
    ''' Class to hold one group of parameter values. '''
    
    def __init__(self, label, t = None, y = None, y_format = None):
        self.label = label
        if t is None: t = odict()
        if y is None: y = odict()
        if y_format is None: y_format = odict()
        self.t = t                      # Time data.
        self.y = y                      # Value data.
        self.y_format = y_format        # Value format data (e.g. Probability, Fraction or Number).
        
    def interpolate(self, tvec = None, pop_label = None):
        ''' Take parameter values and construct an array matching input time vector. '''
        
        # Validate input.
        if pop_label not in self.t.keys(): raise OptimaException('ERROR: Cannot interpolate parameter %s without referring to a proper population label.' % pop_label)
        if tvec is None: raise OptimaException('ERROR: Cannot interpolate parameter %s without providing a time vector.' % self.label)
        if not len(self.t[pop_label]) > 0: raise OptimaException('ERROR: There are no timepoint values for parameter %s, population %s.' % (self.label, pop_label))
        if not len(self.t[pop_label]) == len(self.y[pop_label]): raise OptimaException('ERROR: Parameter %s, population %s, does not have corresponding values and timepoints.' % (self.label, pop_label))

        if len(self.t[pop_label]) == 1:
            output = np.ones(len(tvec))*self.y[pop_label][0]    # Don't bother running interpolation loops if constant.
        else:
            input_t = dcp(self.t[pop_label])
            input_y = dcp(self.y[pop_label])
            
            # Pad the input vectors for interpolation with minimum and maximum timepoint values, to avoid extrapolation blowing up.
            ind_min, t_min = min(enumerate(self.t[pop_label]), key = lambda p: p[1])
            ind_max, t_max = max(enumerate(self.t[pop_label]), key = lambda p: p[1])
            y_at_t_min = self.y[pop_label][ind_min]
            y_at_t_max = self.y[pop_label][ind_max]
            
    #        print t_min
    #        print t_max
    #        print tvec        
            
            if tvec[0] < t_min:
                input_t = np.append(tvec[0], input_t)
                input_y = np.append(y_at_t_min, input_y)
            if tvec[-1] > t_max:
                input_t = np.append(input_t, tvec[-1])
                input_y = np.append(input_y, y_at_t_max)
                
    #        print input_t
    #        print input_y
            
            output = interpolateFunc(input_t, input_y, tvec)
        
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
        
        self.transfers = odict()   # List of inter-population transitions.
        self.transfers['age'] = odict()
    
    def makePars(self, data):
        self.pop_names = data['pops']['name_labels'].keys()
        
        for name in self.pop_names:
            self.pop_labels.append(data['pops']['name_labels'][name])
            
        # Cascade parameters.
        for l, label in enumerate(data['linkpars']):
            self.par_ids[label] = l
            self.pars.append(Parameter(label = label))
            for pop_id in data['linkpars'][label]:
                self.pars[-1].t[pop_id] = data['linkpars'][label][pop_id]['t']
                self.pars[-1].y[pop_id] = data['linkpars'][label][pop_id]['y']
                self.pars[-1].y_format[pop_id] = data['linkpars'][label][pop_id]['y_format']
        
        # Age migrations.
        for source in data['transfers']['aging'].keys():
            for sink in data['transfers']['aging'][source].keys():
                self.transfers['age'][source] = {'target':sink, 'value':float(1/data['pops']['ages'][source]['range'])}