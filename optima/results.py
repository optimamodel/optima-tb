import logging
logger = logging.getLogger(__name__)

from utils import OptimaException

"""




"""

class ResultSet():
    """
    
    A ResultSet contains information:
    
    
    resultset fields: 
            parset_name        name (not index) of parset used to create this result set
            population_labels  minimal data required for populations
            compartment_label  minimal data required for compartments        
            data_sim           simulated data points (format: dataset)
            indices_observed_data    indices for t that correspond to times at which observed data collected 
                                     so that the observed data points can be compared against simulated data
                                     using data_observed - data_sim[indices_observed_data]
            t_step             t_steps in real time 
            dt                 dt used to create this resultset
        Optional: ------------------------------------------------------
            data_observed      datapoints supplied (format: dataset)
            calibration_fit    calibration fit score
            calibration_metric calibration metric used
    
    Example
        if dt = 0.25, and observed data was taken annually
        then indices_observed_data = [0,4,8,..]. This can be determined directly from dt
    
    
    """
    
    def __init__(self,model,parset):
        
        self.parset_name = parset.name

        self.t_step = model.sim_settings['tvec']
        self.indices_observed_data = np.where(self.t_step%1.0==0)
        
        
    def getSimulatedDatapoints(self,pop_id = None, comp_id = None):
        """
        Returns the data points for the simulation that should correspond to times
        of the observed data. This method is intended for 
        
        [Intended]
        If pop_id or comp_id are specified, then the method returns the simulated data points for
        just those populations or compartments.
        Otherwise, the simulated data points for all popuations and compartments are returned.
        
        """
        pass
        
    