#%% Imports

from utils import odict, OptimaException
from interpolation import interpolateFunc

import logging
logger = logging.getLogger(__name__)

from copy import deepcopy as dcp
import numpy as np



#%% Parameter class that stores one array of values converted from raw project data

class Parameter(object):
    ''' Class to hold one set of parameter values disaggregated by populations. '''
    
    def __init__(self, label, t = None, y = None, y_format = None, y_factor = None):
        self.label = label
        
        # These ordered dictionaries have population labels as keys.
        if t is None: t = odict()
        if y is None: y = odict()
        if y_format is None: y_format = odict()
        if y_factor is None: y_factor = odict()
        self.t = t                      # Time data.
        self.y = y                      # Value data.
        self.y_format = y_format        # Value format data (e.g. Probability, Fraction or Number).
        self.y_factor = y_factor        # Scale factor of data (i.e. 1., or None indicating that scaling should not occur during automated calibration
        
    def interpolate(self, tvec = None, pop_label = None):
        ''' Take parameter values and construct an interpolated array corresponding to the input time vector. '''
        
        # Validate input.
        if pop_label not in self.t.keys(): raise OptimaException('ERROR: Cannot interpolate parameter "%s" without referring to a proper population label.' % pop_label)
        if tvec is None: raise OptimaException('ERROR: Cannot interpolate parameter "%s" without providing a time vector.' % self.label)
        if not len(self.t[pop_label]) > 0: raise OptimaException('ERROR: There are no timepoint values for parameter "%s", population "%s".' % (self.label, pop_label))
        if not len(self.t[pop_label]) == len(self.y[pop_label]): raise OptimaException('ERROR: Parameter "%s", population "%s", does not have corresponding values and timepoints.' % (self.label, pop_label))

        if len(self.t[pop_label]) == 1:
            output = np.ones(len(tvec))*self.y[pop_label][0]    # Don't bother running interpolation loops if constant. Good for performance.
        else:
            input_t = dcp(self.t[pop_label])
            input_y = dcp(self.y[pop_label])
            
            # Pad the input vectors for interpolation with minimum and maximum timepoint values, to avoid extrapolated values blowing up.
            ind_min, t_min = min(enumerate(self.t[pop_label]), key = lambda p: p[1])
            ind_max, t_max = max(enumerate(self.t[pop_label]), key = lambda p: p[1])
            y_at_t_min = self.y[pop_label][ind_min]
            y_at_t_max = self.y[pop_label][ind_max]
            
            # This padding effectively keeps edge values constant for desired time ranges larger than data-provided time ranges.
            if tvec[0] < t_min:
                input_t = np.append(tvec[0], input_t)
                input_y = np.append(y_at_t_min, input_y)
            if tvec[-1] > t_max:
                input_t = np.append(input_t, tvec[-1])
                input_y = np.append(input_y, y_at_t_max)
            
            output = interpolateFunc(input_t, input_y, tvec)
        
        return output
    
    def __repr__(self, *args, **kwargs):
        return "Parameter: %s \n t       : %s \n y       : %s \n y_format: %s\n  y_factor: %s\n"%(self.label,self.t,self.y,self.y_format,self.y_factor)
        


#%% Parset class that contains one set of parameters converted from raw project data

class ParameterSet(object):
    ''' Class to hold all parameters. '''
    
    def __init__(self, name='default'):
        self.name = name 
        # TODO: for DK, define what difference is between pop_names and pop_labels, and when 
        # each should be used
        self.pop_names = []         # List of population names.
        self.pop_labels = []        # List of population labels.
        self.pars = odict()
        self.pars['cascade'] = []
        self.pars['characs'] = []
        self.par_ids = {'cascade':{}, 'characs':{}}
        
        self.transfers = odict()    # Dictionary of inter-population transitions.
        
        logging.info("Created ParameterSet: %s"%self.name)
    
    def makePars(self, data):
        self.pop_names = data['pops']['name_labels'].keys()
        
        for name in self.pop_names:
            self.pop_labels.append(data['pops']['name_labels'][name])
            
        # Cascade parameters.
        for l, label in enumerate(data['linkpars']):
            self.par_ids['cascade'][label] = l
            self.pars['cascade'].append(Parameter(label = label))
            for pop_id in data['linkpars'][label]:
                self.pars['cascade'][-1].t[pop_id] = data['linkpars'][label][pop_id]['t']
                self.pars['cascade'][-1].y[pop_id] = data['linkpars'][label][pop_id]['y']
                self.pars['cascade'][-1].y_format[pop_id] = data['linkpars'][label][pop_id]['y_format']
                self.pars['cascade'][-1].y_factor[pop_id] = data['linkpars'][label][pop_id]['y_factor']
                
        # Characteristic parameters (e.g. popsize/prevalence).
        # Despite being mostly data to calibrate against, is still stored in full so as to interpolate initial value.
        for l, label in enumerate(data['characs']):
            self.par_ids['characs'][label] = l
            self.pars['characs'].append(Parameter(label = label))
            for pop_id in data['characs'][label]:
                self.pars['characs'][-1].t[pop_id] = data['characs'][label][pop_id]['t']
                self.pars['characs'][-1].y[pop_id] = data['characs'][label][pop_id]['y']
                self.pars['characs'][-1].y_format[pop_id] = data['characs'][label][pop_id]['y_format']
                self.pars['characs'][-1].y_factor[pop_id] = data['characs'][label][pop_id]['y_factor']
        
        # Migrations, including aging.
        for trans_type in data['transfers'].keys():
            if trans_type not in self.transfers: self.transfers[trans_type] = odict()
            
            for source in data['transfers'][trans_type].keys():
                if source not in self.transfers[trans_type]: self.transfers[trans_type][source] = Parameter(label = trans_type + '_from_' + source)
                for target in data['transfers'][trans_type][source].keys():
                    self.transfers[trans_type][source].t[target] = data['transfers'][trans_type][source][target]['t']
                    self.transfers[trans_type][source].y[target] = data['transfers'][trans_type][source][target]['y']
                    self.transfers[trans_type][source].y_format[target] = data['transfers'][trans_type][source][target]['y_format']
                    self.transfers[trans_type][source].y_factor[target] = data['transfers'][trans_type][source][target]['y_factor']
    
    def __getMinMax(self,y_format):
        if y_format.lower() == 'fraction':
            return (0.,1.)
        elif y_format.lower() in ['number','proportion']:
            return (0.,None)
        else:
            raise OptimaException("Unknown y_format '%s' encountered while returning min-max bounds"%y_format)
    
    
    def extract(self,getMinMax=False):
        """
        Extract parameters values into a list
        
        Note that this method is expecting one y-value per parameter, and does 
        not allow for time-varying parameters. 
        
        To avoid this, users should use the assumption value to specify values, or 
        mark the 'Calibrate?' column in the cascade spreadsheet with "-1" (== settings.DO_NOT_SCALE value)
        
        If getMinMax=True, additionally returns the min and max values for each of the parameters returned.
        This depends on their format.
        """    
        import settings 
        
        paramvec = [] # Not efficient - would prefer: np.zeros(len(self.pop_labels)*len(self.pars['cascade']))
        minmax = []
        casc_labels = []
        index= 0
        for pop_id in self.pop_labels:
            for (j,casc_id) in enumerate(self.par_ids['cascade']): 
                
                if self.pars['cascade'][j].y_factor[pop_id] == settings.DO_NOT_SCALE:
                    continue
                #paramvec[index] = [self.pars['cascade'][j].y[pop_id]]
                paramvec.append(self.pars['cascade'][j].y[pop_id])
                if getMinMax:
                    minmax.append(self.__getMinMax(self.pars['cascade'][j].y_format[pop_id]))
                    casc_labels.append(casc_id)
                index+=1
                
        if getMinMax:
            # Hmm, not a big fan of different return signatures for the one method ... 
            return paramvec,minmax,casc_labels
        return paramvec#[:index]
    
    
    def update(self,paramvec,yearvec=None,y_format_vec=None,y_factor_vec=None):
        """
        Update parameters from a list of values
        
        TODO: added yearvec,y_format_vec,y_factor_vec with intention that they can be used to fully update
        the parameter set cascade parameter
        """
        import settings
        
        index = 0
        for (i,pop_id) in enumerate(self.pop_labels):
            for (j,casc_id) in enumerate(self.par_ids['cascade']): 
                if self.pars['cascade'][j].y_factor[pop_id] == settings.DO_NOT_SCALE:
                    continue
                if len(self.pars['cascade'][j].y[pop_id]) != len(paramvec[index]):
                    raise OptimaException("Could not update parameter set '%s' for pop=%s,cascade=%s as updated parameter has different length."%(self.name,pop_id,casc_id))
                self.pars['cascade'][j].y[pop_id] = paramvec[index]
                index += 1
                
        logger.info("Updated ParameterSet %s with new values"%self.name)
    
    
    def __repr__(self, *args, **kwargs):
        return "ParameterSet: %s \npars: \n%s"%(self.name, self.pars) 
    
    
    