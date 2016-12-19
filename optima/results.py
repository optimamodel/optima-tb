import logging
logger = logging.getLogger(__name__)

from utils import OptimaException, odict

import numpy as np

"""




"""

class ResultSet():
    """
    
    A ResultSet contains information:
    
    
    resultset fields: 
            parset_name        name (not index) of parset used to create this result set
            population_labels  minimal data required for populations
            compartment_label  minimal data required for compartments        
            popdata_sim          simulated compartment data points (format: ?)
            outputs              simulated characteristic data points (format: ?)
                            
            indices_observed_data    indices for t that correspond to times at which observed data collected 
                                     so that the observed data points can be compared against simulated data
                                     using data_observed - data_sim[indices_observed_data]
            t_observed_data    t values for indices_observed_data points
            t_step             t_steps in real time 
            dt                 dt used to create this resultset
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
    
    def __init__(self,model,parset,settings):
        
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
        self.population_label
        self.compartment_labels
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
        
        
        
    def getPopDatapoints(self,pop_id = None, comp_id = None):
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
        pass
        
    
    
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
                
    
    
    