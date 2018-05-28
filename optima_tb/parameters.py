#%% Imports

from optima_tb.utils import odict, OptimaException
from optima_tb.interpolation import interpolateFunc
from optima_tb.databook import getEmptyData
from optima_tb.settings import DO_NOT_SCALE, DEFAULT_YFACTOR

import logging
logger = logging.getLogger(__name__)

from copy import deepcopy as dcp
import numpy as np
from uuid import uuid4 as uuid
from csv import reader, writer
import operator



#%% Parameter class that stores one array of values converted from raw project data

class Parameter(object):
    ''' Class to hold one set of parameter values disaggregated by populations. '''
    
    def __init__(self, label, t = None, y = None, y_format = None, y_factor = None, autocalibrate = None):
        self.label = label
        
        # These ordered dictionaries have population labels as keys.
        if t is None: t = odict()
        if y is None: y = odict()
        if y_format is None: y_format = odict()
        if y_factor is None: y_factor = odict()
        if autocalibrate is None: autocalibrate = odict()
        self.t = t # Time data 
        self.y = y # Value data
        self.y_format = y_format                # Value format data (e.g. Probability, Fraction or Number).
        self.y_factor = y_factor                # Scaling factor of data. Corresponds to different transformations whether format is fraction or number.
        self.autocalibrate = autocalibrate      # A set of boolean flags corresponding to y_factor that denote whether this parameter can be autocalibrated.

    @property
    def pops(self):
        return sorted(self.t.keys())

    def has_values(self, pop_name):
        # Returns True if this Parameter has values specified for the given population
        # Essentially, if this function returns True, then the `interpolate()` method will return usable values
        return np.any(np.isfinite(self.y[pop_name]))

    def clear(self):
        ''' Remove all existing values '''
        for pop in self.pops:
            self.t[pop] = np.zeros((0,))
            self.y[pop] = np.zeros((0,))

    def insertValuePair(self, t, y, pop_label):
        ''' Check if the inserted t value already exists for the population parameter. If not, append y value. If so, overwrite y value. '''
        # Make sure it stays sorted
        if t in self.t[pop_label]:
            self.y[pop_label][self.t[pop_label]==t] = y
        else:
            idx = np.searchsorted(self.t[pop_label],t)
            self.t[pop_label] = np.insert(self.t[pop_label],idx,t)
            self.y[pop_label] = np.insert(self.y[pop_label],idx,y)

    def removeValueAt(self, t, pop_label):
        '''
        # Remove t value if at least one other time point exists    
        '''

        if t in self.t[pop_label]:
            if len(self.t[pop_label]) <= 1:
                return False
            else:
                idx = (self.t[pop_label]==t).nonzero()[0][0]
                self.t[pop_label] = np.delete(self.t[pop_label], idx)
                self.y[pop_label] = np.delete(self.y[pop_label], idx)
                return True
        else:
            return True

    def removeBetween(self,t_remove,pop_label):
        # t is a two element vector [min,max] such that
        # times > min and < max are removed
        # Note that the endpoints are not included!
        original_t = self.t[pop_label]
        for tval in original_t:
            if tval > t_remove[0] and tval < t_remove[1]:
                self.removeValueAt(tval,pop_label)

    def getValueAt(self,t,pop_label):
        # Return parameter value without any scaling
        return self.interpolate(np.array([t]),pop_label)/np.abs(self.y_factor[pop_label])

    def interpolate(self, tvec = None, pop_label = None, extrapolate_nan = False):
        ''' Take parameter values and construct an interpolated array corresponding to the input time vector. '''
        
        # Validate input.
        if pop_label not in self.t.keys(): raise OptimaException('ERROR: Cannot interpolate parameter "%s" without referring to a proper population label.' % pop_label)
        if tvec is None: raise OptimaException('ERROR: Cannot interpolate parameter "%s" without providing a time vector.' % self.label)
        if not len(self.t[pop_label]) > 0: raise OptimaException('ERROR: There are no timepoint values for parameter "%s", population "%s".' % (self.label, pop_label))
        if not len(self.t[pop_label]) == len(self.y[pop_label]): raise OptimaException('ERROR: Parameter "%s", population "%s", does not have corresponding values and timepoints.' % (self.label, pop_label))

        if len(self.t[pop_label]) == 1 and not extrapolate_nan:
            output = np.ones(len(tvec))*(self.y[pop_label][0]*np.abs(self.y_factor[pop_label]))   # Don't bother running interpolation loops if constant. Good for performance.
        else:
            input_t = dcp(self.t[pop_label])
            input_y = dcp(self.y[pop_label])*np.abs(self.y_factor[pop_label])
            
            if not extrapolate_nan:
                # Pad the input vectors for interpolation with minimum and maximum timepoint values, to avoid extrapolated values blowing up.
                ind_min, t_min = min(enumerate(self.t[pop_label]), key = lambda p: p[1])
                ind_max, t_max = max(enumerate(self.t[pop_label]), key = lambda p: p[1])
                y_at_t_min = self.y[pop_label][ind_min]*np.abs(self.y_factor[pop_label])
                y_at_t_max = self.y[pop_label][ind_max]*np.abs(self.y_factor[pop_label])
                
                # This padding effectively keeps edge values constant for desired time ranges larger than data-provided time ranges.
                if tvec[0] < t_min:
                    input_t = np.append(tvec[0], input_t)
                    input_y = np.append(y_at_t_min, input_y)
                if tvec[-1] > t_max:
                    input_t = np.append(input_t, tvec[-1])
                    input_y = np.append(input_y, y_at_t_max)
            
            elif len(input_t) == 1:       # The interpolation function currently complains about single data points, so this is the simplest hack for nan extrapolation.
                input_t = np.append(input_t, input_t[-1] + 1e-12)
                input_y = np.append(input_y, input_y[-1])
            
            output = interpolateFunc(x = input_t, y = input_y, xnew = tvec, extrapolate_nan = extrapolate_nan)
        
        return output
    
    def __repr__(self, *args, **kwargs):
        return "Parameter: %s\n\nt\n%s\ny\n%s\ny_format\n%s\ny_factor\n%s\n"%(self.label,self.t,self.y,self.y_format,self.y_factor)
        

#%% Parset class that contains one set of parameters converted from raw project data

class ParameterSet(object):
    ''' Class to hold all parameters. '''
    
    def __init__(self, name='default'):
        self.name = name 
        self.uid = uuid()

        self.pop_names = []         # List of population names.
                                    # Names are used only for user interface and can be more elaborate than simple labels.
        self.pop_labels = []        # List of population labels.
                                    # Labels are used throughout the codebase as variable names (technically dict keys).
        self.pars = odict()
        self.pars['cascade'] = []
        self.pars['characs'] = []
        self.par_ids = {'cascade':{}, 'characs':{}}
        
        self.transfers = odict()    # Dictionary of inter-population transitions.
        self.contacts = odict()     # Dictionary of inter-population interaction weights.
        
        logging.info("Created ParameterSet: %s"%self.name)
        
    def getPar(self, label):
        for par_type in ['cascade','characs']:
            if label in self.par_ids[par_type].keys():
                return self.pars[par_type][self.par_ids[par_type][label]]
        raise OptimaException('ERROR: Label "%s" cannot be found in parameter set "%s" as either a cascade parameter or characteristic.' % (label, self.name))
    
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
                if data['linkpars'][label][pop_id]['y_factor'] == DO_NOT_SCALE:
                    self.pars['cascade'][-1].y_factor[pop_id] = DEFAULT_YFACTOR
                    self.pars['cascade'][-1].autocalibrate[pop_id] = False
                else:
                    self.pars['cascade'][-1].y_factor[pop_id] = data['linkpars'][label][pop_id]['y_factor']
                    self.pars['cascade'][-1].autocalibrate[pop_id] = True
                
        # Characteristic parameters (e.g. popsize/prevalence).
        # Despite being mostly data to calibrate against, is still stored in full so as to interpolate initial value.
        for l, label in enumerate(data['characs']):
            self.par_ids['characs'][label] = l
            self.pars['characs'].append(Parameter(label = label))
            for pop_id in data['characs'][label]:
                self.pars['characs'][l].t[pop_id] = data['characs'][label][pop_id]['t']
                self.pars['characs'][l].y[pop_id] = data['characs'][label][pop_id]['y']
                self.pars['characs'][l].y_format[pop_id] = data['characs'][label][pop_id]['y_format']
                if data['characs'][label][pop_id]['y_factor'] == DO_NOT_SCALE:
                    self.pars['characs'][-1].y_factor[pop_id] = DEFAULT_YFACTOR
                    self.pars['characs'][-1].autocalibrate[pop_id] = False
                else:
                    self.pars['characs'][-1].y_factor[pop_id] = data['characs'][label][pop_id]['y_factor']
                    self.pars['characs'][-1].autocalibrate[pop_id] = True
            
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
                    
        self.contacts = dcp(data['contacts'])   # Simple copying of the contacts structure into data. No need to be an object.
    
    
    def inflate(self,tvec):
        """
        Inflates cascade parameters only
        
        """
        for (id,par) in enumerate(self.pars['cascade']):
           
            for j,pop_label in enumerate(self.pop_labels):        
                self.pars['cascade'][id].y[j] = par.interpolate(tvec = tvec, pop_label = pop_label)
                self.pars['cascade'][id].t[j] = tvec
           
    
    def __getMinMax(self,y_format):
        if y_format.lower() == 'fraction':
            return (0.,1.)
        elif y_format.lower() in ['number','proportion']:
            return (0.,np.inf)
        else:
            raise OptimaException("Unknown y_format '%s' encountered while returning min-max bounds"%y_format)
    
    
    def extract(self,settings=None,getMinMax=False,getYFactor=False):
        """
        Extract parameters values into a list
        
        Note that this method is expecting one y-value per parameter, and does 
        not allow for time-varying parameters. If time-varying parameters are expected, then 
        getYFactor will return the scaling factor (y_factor) that can be used to manipulate
        
        To avoid values from being extracted, users should use the assumption value to specify values, or 
        mark the 'Autocalibrate' column in the cascade spreadsheet with "-1" (== settings.DO_NOT_SCALE value)
        
        Note that parameters defined as functions will not be extracted.
        
        If getMinMax=True, this function additionally returns the min and max values for each of the parameters returned,
        depending on y_format. Each minmax value is a tuple of form (min,max). Note that min/max values can be either a float, int, or None.
        
        Params:
            getMinMax     boolean. Additionally return min and max values for each of the parameters returned in format: 
                            @TODO: document format for this and return statements
            getYFactor    boolean: 
                            
        Return:
            paramvec
            minmax        minmax values. [] if getMinMax was False, else a list of tuples: [(min,max)]
            casc_labels    
        """    
#        import optima_tb.settings as settings 
        
        paramvec = [] # Not efficient - would prefer: np.zeros(len(self.pop_labels)*len(self.pars['cascade']))
        minmax = []
        casc_labels = []
        index= 0
        for pop_id in self.pop_labels:
            #for (j,casc_id) in enumerate(self.par_ids['cascade']): 
            for (casc_id,j) in sorted(self.par_ids['cascade'].items(), key=operator.itemgetter(1)):
#                print casc_id
#                print j

                # Do not extract parameters that are functions (if settings are available) or otherwise marked not to be calibrated.
                if self.pars['cascade'][j].autocalibrate[pop_id] == False or (not settings is None and casc_id in settings.par_funcs):
                    continue
                #paramvec[index] = [self.pars['cascade'][j].y[pop_id]]
                if getYFactor:
                    paramvec.append(self.pars['cascade'][j].y_factor[pop_id])
                else:
                    paramvec.append(self.pars['cascade'][j].y[pop_id])
                    
                if getMinMax:
                    minmax_bounds = self.__getMinMax(self.pars['cascade'][j].y_format[pop_id])
                    if getYFactor and minmax_bounds[1] is not None:
                        # use the bounds to calculate the corresponding minmax values for the bounds: 
#                        # possibilities are (0.,1.) for fraction, and (0,None) for number and proportion.
#                        # - yfactor-min has to remain 0, so leave
#                        yval_max= np.max(self.pars['cascade'][j].y[pop_id])
#                        if yval_max == 0:
#                            logger.info("ParameterSet: max value of 0 observed for fraction for casc_id=%s"%casc_id)
#                            yval_max = settings.TOLERANCE
#                        tmp_upper = minmax_bounds[1]
#                        tmp_upper /= yval_max
#                        minmax_bounds = (minmax_bounds[0], tmp_upper)
                        minmax_bounds = (0., np.inf)      # There is probably no reason to enforce bounds on y_factor in this layer as fraction ranges are handled in model.py.
                    # if we're grabbing the minmax values for the y-values directly, then use the values
                    minmax.append(minmax_bounds)
                    # also grab the cascade labels for debugging and logging purposes
                    casc_labels.append((casc_id,pop_id))
                index+=1
                
        if getMinMax:
            return paramvec,minmax,casc_labels
        
        return paramvec,None,casc_labels 
    
    def extractEntryPoints(self,proj_settings,useInitCompartments=False):
        """
        Extract initial compartments: 
        """
        import optima_tb.settings as settings
        
        init_compartments = []
        charac_labels = []
        if useInitCompartments:
            for pop_id in self.pop_labels:
                for (charac_id,j) in sorted(self.par_ids['characs'].items(), key=operator.itemgetter(1)):
                    #if 'entry_point' in proj_settings.charac_specs[charac_id].keys() and self.pars['characs'][j].y_factor[pop_id] != settings.DO_NOT_SCALE:
#                    print charac_id
#                    print j
                    if 'entry_point' in proj_settings.charac_specs[charac_id].keys() and self.pars['characs'][j].autocalibrate[pop_id] == True:
                        init_compartments.append(self.pars['characs'][j].y_factor[pop_id])
                        charac_labels.append((charac_id,pop_id))                        
        return init_compartments,charac_labels
         
    
    def update(self, values, par_pop_labels, isYFactor=False):
        """
        Update parameters from a list of values, given a corresponding list of parameter labels and population labels, arranged in pairs.

        values - list of values to insert
        par_pop_labels - list of corresponding tuples (par_label,pop_label)

        TODO: extend function so that can also add years or format values
        TODO: remove index ...?
        """
        import optima_tb.settings as settings
        
        index = 0
        for par_pop_label in par_pop_labels:
            par_label = par_pop_label[0]
            pop_label = par_pop_label[1]
            par = self.getPar(par_label)    # Returns a parameter or characteristic as required.        
                
            if isYFactor:
                
#                if par.autocalibrate[pop_label]:
#                    logger.info("ParameterSet %s is acknowledging that it can calibrate the y-factor for parameter %s, population %s." % (self.name, par_label, pop_label))
#                else:
#                    logger.info("ParameterSet %s has not modifed its y-factor value for parameter %s, population %s, due to y-factor constraints." % (self.name, par_label, pop_label))
#                    continue
                
                if par.y_factor[pop_label] != values[index]:
                    logger.info("ParameterSet %s is updating its y-factor for parameter %s, population %s from %f to %f." % (self.name, par_label, pop_label, par.y_factor[pop_label], values[index]))
                
                par.y_factor[pop_label] = values[index]
            else:
                par.y[pop_label] = values[index]
            index += 1
            
        
#        index = 0
#        for (i,pop_id) in enumerate(self.pop_labels):
#            #for (j,casc_id) in enumerate(self.par_ids['cascade']): 
#            for (casc_id,j) in sorted(self.par_ids['cascade'].items(), key=operator.itemgetter(1)):
#                # perform checks
#                if self.pars['cascade'][j].y_factor[pop_id] == settings.DO_NOT_SCALE:
#                    continue
#                if not isYFactor and len(self.pars['cascade'][j].y[pop_id]) != len(paramvec[index]):
#                    raise OptimaException("Could not update parameter set '%s' for pop=%s,cascade=%s as updated parameter has different length."%(self.name,pop_id,casc_id))
#                # update y or y_factor, based on parameters
#                if isYFactor:
#                    self.pars['cascade'][j].y_factor[pop_id] = paramvec[index]
#                else:
#                    self.pars['cascade'][j].y[pop_id] = paramvec[index]
#                # finally, update index count
#                index += 1
                
        # logger.info('Updated ParameterSet "%s" with new values.' % self.name)
    
#    def updateEntryPoints(self,proj_settings,compartment_t0,charac_labels):
#        """
#        
#        
#        """
#        import optima_tb.settings as settings 
#        
#        index = 0 
#        for pop_id in self.pop_labels:
#            for (charac_id,j) in sorted(self.par_ids['characs'].items(), key=operator.itemgetter(1)):
#                if 'entry_point' in proj_settings.charac_specs[charac_id].keys() and self.pars['characs'][j].y_factor[pop_id] != settings.DO_NOT_SCALE:
#                    if charac_labels[index] == charac_id:
#                        self.pars['characs'][j].y_factor[pop_id]= compartment_t0[index]
#                        #self.pars['characs'][j].y[pop_id][0] *= compartment_t0[index]
#                        #print "initial for %s is %g"%(charac_id,self.pars['characs'][j].y[pop_id][0] )
#                        index += 1
#                    else:
#                        logger.debug("Updating entry points: characteristic label for entry point doesn't match updated value [%s,%s]"%(charac_labels[index],charac_id))
                        
        
    #TODO: Fix as the DO_NOT_SCALE setting no longer represents an absence of rescaling; attribute autocalibrate does.
    def _updateFromYFactor(self):
        """
        Goes through each of the cascade parameters and updates the y-values to be
        y-value*y_factor. This is used for example in autocalibration (which can scale y_factor),
        which may have set the y_factor to a non-unary value but then needs to update y values at the end. 
        
        The purpose of doing this is to ensure that y_values are correct, if inspected, to avoid errors if a process
        only checks y_value rather than y_value*y_factor.
        
        """
        import optima_tb.settings as settings
        
        for (i,pop_id) in enumerate(self.pop_labels):
            for (j,casc_id) in enumerate(self.par_ids['cascade']): 
                if self.pars['cascade'][j].y_factor[pop_id] == settings.DO_NOT_SCALE:
                    continue
                self.pars['cascade'][j].y[pop_id] *= self.pars['cascade'][j].y_factor[pop_id] 
                self.pars['cascade'][j].y_factor[pop_id] = 1.0
                
        for pop_id in self.pop_labels:
            for (charac_id,j) in sorted(self.par_ids['characs'].items(), key=operator.itemgetter(1)):
                if self.pars['characs'][j].y_factor[pop_id] == settings.DO_NOT_SCALE:
                    continue
                self.pars['characs'][j].y[pop_id] *= self.pars['characs'][j].y_factor[pop_id] 
                self.pars['characs'][j].y_factor[pop_id] = 1.0
        
        
        
    
    def __repr__(self, *args, **kwargs):
        return "ParameterSet: %s \npars: \n%s"%(self.name, self.pars) 
    
    
    def __add__(a,b):
        """
        Add two parameter sets together, with the value of any parameter
        cascade in ParameterSet 'b' being added to the corresponding value 
        in ParameterSet 'a'. If no corresponding value in a exists, it is inserted.
        
        As only the cascade values in 'b' are copied over, it is mandatory at this
        stage that 'a' should be the ParameterSet kept for it's characteristics and
        transfer definitions.
        
        Throws an OptimaException if parameters with matching names have different
        y_formats
        
        Usage:
        default  = ParameterSet()
        largeVal = ParameterSet()
        c = default + largeVal
        print type(c)
        > "ParameterSet"
        
        """
        logger.debug("Adding two parameter sets together: %s + %s"%(a.name,b.name))
        c = dcp(a)
        for par_name, b_index in b.par_ids['cascade'].iteritems():
            
            # find corresponding par_id in c
            c_index = c.par_ids['cascade'][par_name]
            # for each population referenced:
            for pop in b.pars['cascade'][b_index].t.keys():     
                # check that the y_format matches: if not, throw an error
                if b.pars['cascade'][b_index].y_format[pop] != c.pars['cascade'][c_index].y_format[pop]:
                    raise OptimaException("ERROR: trying to combine two Parameters with different y_formats: ")
                # add or insert value of b into c
                for i,t_val in enumerate(b.pars['cascade'][b_index].t[pop]):
                    
                    tmp =  c.pars['cascade'][c_index].t[pop]
                    if t_val in c.pars['cascade'][c_index].t[pop]: 
                        mask = tmp == t_val
                        c.pars['cascade'][c_index].y[pop][mask] += b.pars['cascade'][b_index].y[pop][i]
                    else:
                        c.pars['cascade'][c_index].y[pop] = np.append(c.pars['cascade'][c_index].y[pop],[b.pars['cascade'][b_index].y[pop][i]])
                        c.pars['cascade'][c_index].t[pop] = np.append(c.pars['cascade'][c_index].t[pop],[t_val])
                
                # correct for min/max, based on format type: as presumably 'a' was already correct, 
                # we only need to check that the min max wasn't violated when we're adding values together, and therefore
                # only have to check when we've added something from b. 
                minmax = c.__getMinMax(a.pars['cascade'][c_index].y_format[pop])
                if minmax[0] is not None and sum(c.pars['cascade'][c_index].y[pop]<minmax[0]) > 0: 
                    logger.debug("ParameterSet.__add__ : Observed min that is less than accepted min value for parameter type '%s': cascade label=%s for population=%s"%(c.pars['cascade'][c_index].y_format[pop],par_name,pop))
                    c.pars['cascade'][c_index].y[pop][c.pars['cascade'][c_index].y[pop]<minmax[0]] = minmax[0]
                if minmax[1] is not None and sum(c.pars['cascade'][c_index].y[pop]>minmax[1]) > 0: 
                    logger.debug("ParameterSet.__add__ : Observed max that is less than accepted max value for parameter type '%s': cascade label=%s for population=%s"%(c.pars['cascade'][c_index].y_format[pop],par_name,pop))
                    c.pars['cascade'][c_index].y[pop][c.pars['cascade'][c_index].y[pop]>minmax[1]] = minmax[1]
        return c
        

def __read_block(csvdata,pops,index_start,index_end):
    import itertools

    datadict = odict()
    
    c_label = None
    p_label = None
    y_format = None
    y_factor = None
    ts = []
    ys = []
    for i in range(index_start, index_end):
        row = csvdata[i]
        if len(row)!=2: # either a cascade label or population label
            # clear out buffer from previous rows for the previous c_label and p_labels
            if c_label is not None and p_label is not None:
                datadict[c_label][p_label]['t'] = np.array(ts)
                datadict[c_label][p_label]['y'] = np.array(ys)
                datadict[c_label][p_label]['y_factor'] = y_factor
                datadict[c_label][p_label]['y_format'] = y_format
                p_label = None
            # set up for next c_label or p_label
            if len(row)==1: # and not row[0] in pops:
                c_label = row[0]
                datadict[c_label] = odict()
            elif len(row)==3:
                p_label = row[0]
                y_format, y_factor = row[1:]
                datadict[c_label][p_label] = odict()
                ys = []
                ts = []
                
        else:
            ts.append(float(row[0]))
            ys.append(float(row[1]))
    # finish off and empty buffer: 
    datadict[c_label][p_label]['t'] = np.array(ts)
    datadict[c_label][p_label]['y'] = np.array(ys)
    datadict[c_label][p_label]['y_factor'] = y_factor
    datadict[c_label][p_label]['y_format'] = y_format
    return datadict
                

def __write_block(csvwriter,datadict):
    for cid in datadict:
        csvwriter.writerow([cid.label])
        for (p,tdata) in cid.t.iteritems():
            csvwriter.writerow([p,cid.y_format[p],cid.y_factor[p]])
            for (i,ti) in enumerate(tdata):
                csvwriter.writerow([ti,cid.y[p][i]])
         
def export_paramset(parset):
    """
    Saves to file: <parset.name>.csv
    
    Format:
        Populations
        label1,name1
        label2,name2
        Cascade
        <casc_label1>
        <Pop1>,<y_format>,<y_factor>
        t1,y1
        t2,y2
        t3,y3
        <Pop2>,<y_format>,<y_factor>
        t1,y1
        t2,y2
        t3,y3
        <casc_label2>
        <Pop1>,<y_format>,<y_factor>
        t1,y1
        <Pop2>,<y_format>,<y_factor>
        t1,y1
        t2,y2
        ...
        Characteristic
        <char_label1>
        <Pop1>,<y_format>,<y_factor>
        t1,y1
        <Pop2>,<y_format>,<y_factor>
        t1,y1
        <char_label2>
        <Pop1>,<y_format>,<y_factor>
        t1,y1
        <Pop2>,<y_format>,<y_factor>
        t1,y1
        ...
    """
    # setup writer
    filename = parset.name + ".csv"
    with open(filename, 'wb') as csvfile:
        pswriter = writer(csvfile, delimiter=',')
        # write populations
        pswriter.writerow(['Populations'])
        for (i,p) in enumerate(parset.pop_labels):
            pswriter.writerow([p,parset.pop_names[i]])
        # write cascade
        pswriter.writerow(['Cascade'])
        __write_block(pswriter, parset.pars['cascade'])
        # write characteristics
        pswriter.writerow(['Characteristics'])
        __write_block(pswriter, parset.pars['characs'])
        
    # finish
    logger.info("Exported Parameterset '%s' to %s"%(parset.name,filename))
    


def load_paramset(parset_filename):
    """
    Loads and builds paramset of a csv file, assuming structure as specified in export_paramset(). 
    In the new parset, the parset.name is set as the filename (without the csv extension) and a new uid is assigned.
    
    Params:
        parset_filename    filename as string
        
        
    TODO: add 
    """
    import os
    import itertools
    parset_filename = os.path.abspath(parset_filename)
    
    # Setup data structure
    data = getEmptyData()
    index = 0
    indices = {}
    # Fill in data from csv
    with open(parset_filename, 'rb') as csvfile:
        psreader = reader(csvfile, delimiter=',')
        csvdata = list(psreader)
    row_count = len(csvdata)
    # obtain row indices for the three sections: 'Populations','Cascade' and 'Characteristics'
    for row in csvdata:
        if len(row) == 1 and row[0] in ['Populations','Cascade','Characteristics']:
            indices[row[0]] = index
        index += 1
       
    with open(parset_filename, 'rb') as csvfile:
        # Read in populations
        psreader = reader(csvfile, delimiter=',')
        for i in range(indices['Populations']+1,indices['Cascade']):
            row = csvdata[i]
            data['pops']['name_labels'][row[1]] = row[0]
            data['pops']['label_names'][row[0]] = row[1]
        pops = data['pops']['label_names'].keys()

        # Read in cascade
        data['linkpars'] = __read_block(csvdata, pops, indices['Cascade']+1, indices['Characteristics'])
        
        # Read in characteristics
        data['characs'] = __read_block(csvdata, pops, indices['Characteristics']+1,row_count)
        
        
    # Use ps.makePars method to automagically construct the parset from the data
    parset_name = os.path.splitext(parset_filename)[0]
    ps = ParameterSet(name=parset_name)
    ps.makePars(data)
    return ps
    
