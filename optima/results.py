import logging
logger = logging.getLogger(__name__)

from utils import OptimaException, odict, defaultrepr, objrepr
import numpy as np


class Result(object):
    '''
    Class to contain one set of results
    '''
    def __init__(self, name = None, ispercentage = False, pops = None, tot = None):
        self.name = name                 # Name of this parameter
        self.ispercentage = ispercentage # Whether or not the result is a percentage
        self.pops = pops                 # The model result by population, if available
        self.tot = tot                   # The model result total, if available
            
    def __repr__(self):
        ''' Print out useful information when called '''
        output = defaultrepr(self)
        return output
        

#%% Resultset class that contains one set of results
class ResultSet(object):
    """
    Class to hold one set of Results
    Fields: 
            name               name of result
            parset_name        name (not index) of parset used to create this result set
            parset_id          uuid of relevant parset (in case of duplication or missing parset_name)
            sim_settings       settings for simulation
            pop_labels         population tags for which siumulation was run (e.g. [])
            char_label         characteristics tags for which simulation was run
            ind_observed_data  indices for t that correspond to times at which observed data was collected
                               so that the observed data points can be compared against simulated data
                               using data_observed - data_sim[indices_observed_data]
            t_observed_data    t values for indices_observed_data points
                               popdata_sim        simulated compartment data points (format: ?)
            t_step             t_steps in real time 
            dt                 dt used to create this resultset
            outputs            simulated characteristic data points (format: ?)
            m_pops             totals per pop per compartment
            
        Optional: ------------------------------------------------------
            data_observed      datapoints supplied (format: dataset)
            calibration_fit    calibration fit score
            calibration_metric calibration metric used
    
        Example
            if dt = 0.25, and observed data was taken annually
            then 
                indices_observed_data = [0,4,8,..]. This can be determined directly from dt
                t_observed_data       = [2000.,2001.,2002., ...]
    """
    
    def __init__(self, model, parset, settings):
        
        self.name = 'results - ' + parset.name
        self.parset_name = parset.name
        self.parset_id  = parset.uid
        
        self.dt = settings.tvec_dt
        self.t_step = model.sim_settings['tvec']
        self.indices_observed_data = np.where(self.t_step%1.0==0)
        self.t_observed_data = self.t_step[self.indices_observed_data]
        
        self.outputs = model.calculateOutputs(settings = settings)
        
        # Set up for future use
        self.calibration_fit = None
        self.calibration_score = None
        
        """
        # remaining fields to be set:
        self.popdata_sim
        self.chardata_sim
        self.data_observed
        """
        
        # work-in-progress: in time, these sections should be removed and only the data
        # points we are interested should be listed
        self.m_pops = model.pops
        self.sim_settings = model.sim_settings
        self.pop_labels = self.outputs[0].keys()
        self.comp_labels = settings.node_names
        self.comp_specs = settings.node_specs
        self.comp_label_names = self.__generateLabelNames(self.comp_specs.keys(),self.comp_labels)
        self.char_labels = self.outputs.keys() # definitely need a better way of determining these
        # /work-in-progress
    
    def __repr__(self):
        ''' Print out useful information when called'''
        #output = '============================================================\n'
        #output += '      Project name: %s\n'    % (self.parset_name if self.parset_name is not None else None)
        #output += '      Date created: %s\n'    % getdate(self.created)
        #output += '               UID: %s\n'    % self.uid
        output = '====================================================================================\n'
        output += objrepr(self)
        return output
        
        
    def extractSimulationData(self):
        """
        
        Currently, only extract the simulated data for characteristics, as we should only use
        this in fitting for calibration.
        """
        pass
        
    def __generateLabelNames(self,labels,names):
        """
        
        Note: this could potentially go in utils or another location, as it's not really 
        specific to results.py. 
        """
        label_names = odict()
        for (i,label) in enumerate(labels):
            label_names[label] = names[i]
        return label_names
        
        
        
    def getPopDatapoints(self, pop_labels = None, comp_labels = None):
        """
        Returns the data points for the simulation, for each compartment, 
        that should correspond to times of the observed data. This method is intended for
        use with calibration. 
        
        [Comment to DK: will we need compartments at all during calibration, scenarios or optimization?
        Possibly for post-hoc analysis, but perhaps it would be worthwhile to assume no, and only
        for some cases should a user - perhaps via a flag in settings - save all data.]
        
        [Intended]
        If pop_id or comp_id are specified, then the method returns the simulated data points for
        just those populations or compartments.
        Otherwise, the simulated data points for all popuations and compartments are returned.

        @param pop_id: population id, either as single id or as list
        @param comp_id: compartment id, either as single id or as list
                
        """
        if pop_labels is not None:
            if isinstance(pop_labels, list):
                pops = pop_labels
            else:
                pops = [pop_labels]
        else:
            pops = self.pop_labels
            
        if comp_labels is not None:
            if isinstance(comp_labels, list):
                comp = comp_labels
            else:
                comp = [comp_labels]
        else:
            comp = self.comp_labels
        
        # this currently uses odict for the output ... Hmm, not sure whether this is the best
        datapoints = odict() 
        for cj in comp:
            datapoints[cj] = odict() 
            for pi in pops:
                datapoints[cj][pi] = self.outputs[cj][pi][self.indices_observed_data]
        return datapoints
        
    
    
    def getCharacteristicDatapoints(self,pop_label = None, char_label = None):
        """
        Returns the data points for the simulation, for each characteristics, 
        that should correspond to times of the observed data. This method is intended for
        use with calibration. 
        
        [Intended]
        If pop_id or char_id are specified, then the method returns the simulated data points for
        just those populations or compartments.
        Otherwise, the simulated data points for all popuations and compartments are returned.
        
        @param pop_id: population id, either as single id or as list
        @param char_id: characteristic id, either as single label or as list of labels
        
        """
        if pop_label is not None:
            if isinstance(pop_label, list):
                pops = pop_label
            else:
                pops = [pop_label]
        else:
            pops = self.pop_labels
            
        if char_label is not None:
            if isinstance(char_label, list):
                chars = char_label
            else:
                chars = [char_label]
        else:
            chars = self.char_labels
        
        # this currently uses odict for the output ... Hmm, not sure whether this is the best
        datapoints = odict() 
        for cj in chars:
            datapoints[cj] = odict() 
            for pi in pops:
                datapoints[cj][pi] = self.outputs[cj][pi][self.indices_observed_data]
        return datapoints
                
    def export(self, filestem=None, bypop=False, sep=',', writetofile=True):
        """
        Export method for characteristics results obtained from a simulation that should correspond 
        to times of the observed data (i.e. annually). This method is intended for use with runSim 
        currently and will be extended to include optimization and scenario results.
        """
        if filestem is None:  # Doesn't include extension, hence filestem
            filestem = self.name
        filename = filestem + '.csv'
        npts = len(self.t_observed_data)
        keys = self.char_labels
        
        output = sep.join(['Indicator','Population'] + ['%i'%t for t in self.t_observed_data]) # Create header and years
        for key in keys:
            if bypop: output += '\n' # Add a line break between different indicators
            if bypop: popkeys = ['tot']+self.pop_labels # include total even for bypop
            else:     popkeys = ['tot']
            for pk, popkey in enumerate(popkeys):
                output += '\n'
                if bypop==False and popkey!='tot': data = self.outputs[key][popkey][:]
                elif bypop and popkey!='tot': data = self.outputs[key][popkey][:]
                else: data = self.outputs[key][0][:]
                output += key+sep+popkey+sep
                for t in range(npts):
                    output += ('%i'+sep) % data[t]
            
        if writetofile: 
            with open(filename, 'w') as f: f.write(output)
            print 'Results exported to "%s"' % (filename)
            return None
        else:
            return output
        
        