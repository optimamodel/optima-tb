#%% Imports

from utils import odict, OptimaException
from interpolation import interpolateFunc
from databook import getEmptyData

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
        self.uid   = uuid()
    
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
                self.pars['characs'][l].t[pop_id] = data['characs'][label][pop_id]['t']
                self.pars['characs'][l].y[pop_id] = data['characs'][label][pop_id]['y']
                self.pars['characs'][l].y_format[pop_id] = data['characs'][label][pop_id]['y_format']
                self.pars['characs'][l].y_factor[pop_id] = data['characs'][label][pop_id]['y_factor']
            
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
    
    
    def extract(self,getMinMax=False,getYFactor=False):
        """
        Extract parameters values into a list
        
        Note that this method is expecting one y-value per parameter, and does 
        not allow for time-varying parameters. 
        
        To avoid this, users should use the assumption value to specify values, or 
        mark the 'Calibrate?' column in the cascade spreadsheet with "-1" (== settings.DO_NOT_SCALE value)
        
        If getMinMax=True, additionally returns the min and max values for each of the parameters returned.
        This depends on their format.
        
        Params:
            getMinMax     boolean. Additionally return min and max values for each of the parameters returned in format: 
                            @TODO: document format for this and return statements
            getYFactor    boolean: 
                            
        Return:
            paramvec
            minmax        minmax values. None if getMinMax was False
            casc_labels    
        """    
        import settings 
        
        paramvec = [] # Not efficient - would prefer: np.zeros(len(self.pop_labels)*len(self.pars['cascade']))
        minmax = []
        casc_labels = []
        index= 0
        for pop_id in self.pop_labels:
            #for (j,casc_id) in enumerate(self.par_ids['cascade']): 
            for (casc_id,j) in sorted(self.par_ids['cascade'].items(), key=operator.itemgetter(1)):
                
                if self.pars['cascade'][j].y_factor[pop_id] == settings.DO_NOT_SCALE:
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
                        # possibilities are (0.,1.) for fraction, and (0,None) for number and proportion.
                        # - yfactor-min has to remain 0, so leave
                        yval_max= np.max(self.pars['cascade'][j].y[pop_id])
                        if yval_max == 0:
                            logger.info("ParameterSet: max value of 0 observed for fraction for casc_id=%s"%casc_id)
                            yval_max = setting.TOLERANCE
                        tmp_upper = minmax_bounds[1]
                        tmp_upper /= yval_max
                        minmax_bounds = (minmax_bounds[0], tmp_upper)
                    # if we're grabbing the minmax values for the y-values directly, then use the values
                    minmax.append(minmax_bounds)
                    # also grab the cascade labels for debugging and logging purposes
                    casc_labels.append(casc_id)
                index+=1
                
        if getMinMax:
            return paramvec,minmax,casc_labels
        
        return paramvec,None,casc_labels 
    
    def extractEntryPoints(self,proj_settings,useInitCompartments=False):
        """
        Extract initial compartments: 
        """
        import settings
        
        init_compartments = []
        charac_labels = []
        if useInitCompartments:
            for pop_id in self.pop_labels:
                for (charac_id,j) in sorted(self.par_ids['characs'].items(), key=operator.itemgetter(1)):
                    if 'entry_point' in proj_settings.charac_specs[charac_id].keys() and self.pars['characs'][j].y_factor[pop_id] != settings.DO_NOT_SCALE:
                        init_compartments.append(self.pars['characs'][j].y[pop_id][0])
                        charac_labels.append(charac_id)                        
        return init_compartments,charac_labels
         
    
    def update(self,paramvec,isYFactor=False):
        """
        Update parameters from a list of values
        
        Params:
            paramvec    new values
            y_factor_vec
        
        TODO: improve this function to future proof it and make it more robust i.e. reference only by population and compartment id.
        TODO: extend function so that can also add years or format values
        TODO: remove index ...?
        """
        import settings
        
        index = 0
        for (i,pop_id) in enumerate(self.pop_labels):
            #for (j,casc_id) in enumerate(self.par_ids['cascade']): 
            for (casc_id,j) in sorted(self.par_ids['cascade'].items(), key=operator.itemgetter(1)):
                # perform checks
                if self.pars['cascade'][j].y_factor[pop_id] == settings.DO_NOT_SCALE:
                    continue
                if not isYFactor and len(self.pars['cascade'][j].y[pop_id]) != len(paramvec[index]):
                    raise OptimaException("Could not update parameter set '%s' for pop=%s,cascade=%s as updated parameter has different length."%(self.name,pop_id,casc_id))
                # update y or y_factor, based on parameters
                if isYFactor:
                    self.pars['cascade'][j].y_factor[pop_id] = paramvec[index]
                else:
                    self.pars['cascade'][j].y[pop_id] = paramvec[index]
                # finally, update index count
                index += 1
                
        logger.info("Updated ParameterSet %s with new values"%self.name)
    
    def updateEntryPoints(self,settings,compartment_t0,charac_labels):
        """
        
        
        """
        index = 0 
        for pop_id in self.pop_labels:
            for (charac_id,j) in sorted(self.par_ids['characs'].items(), key=operator.itemgetter(1)):
                if 'entry_point' in settings.charac_specs[charac_id].keys():# and self.pars['characs'][j].y_factor[pop_id] == settings.DO_NOT_SCALE:
                    
                    if charac_labels[index] == charac_id:
                        self.pars['characs'][j].y[pop_id][0] = compartment_t0[index]
                        index += 1
                    else:
                        pass #logger.debug("Updating entry points: characteristic label for entry point doesn't match updated value [%s,%s]"%(charac_labels[index],charac_id))
                        
        
    
    def _updateFromYFactor(self):
        """
        Goes through each of the cascade parameters and updates the y-values to be
        y-value*y_factor. This is used for example in autocalibration (which uses y_factor),
        which may have set the y_factor to a non-unary value but then needs to update y values at the end. 
        The benefit of doing this is that most values will inspect y_value and 
        not y_value*y_factor; updating it thus avoids any instances where checks don't take both into account.
    
        
        """
        import settings
        
        for (i,pop_id) in enumerate(self.pop_labels):
            for (j,casc_id) in enumerate(self.par_ids['cascade']): 
                if self.pars['cascade'][j].y_factor[pop_id] == settings.DO_NOT_SCALE:
                    continue
                self.pars['cascade'][j].y[pop_id] *= self.pars['cascade'][j].y_factor[pop_id] 
                self.pars['cascade'][j].y_factor[pop_id] = 1.0
        
        
        
    
    def __repr__(self, *args, **kwargs):
        return "ParameterSet: %s \npars: \n%s"%(self.name, self.pars) 


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
    