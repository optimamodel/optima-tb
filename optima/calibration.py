import logging
logger = logging.getLogger(__name__)

from utils import OptimaException, tic, toc

import numpy as np

"""

Calibration and sensitivity analysis

"""


def performSensitivityAnalysis(project,transmission,epsilon=0.05):
    """
    
    
    """
    pass

def calculateFit(model,parset_name):
    """
    
    
    TODO implement
    """
    pass
    
def _calculateFitscore(y_obs, y_fit,metric="minsquare"):
    """
    
    TODO implement calculateFit
    """
    availfns = globals().copy()
    availfns.update(locals())
    try:
        availfns.get('_calc_%s'%metric)(y_obs,y_fit)
    except:
        raise NotImplementedError("No method associated with _calc_%s (calibration.py)"%metric)


def _calc_minsquare(y_obs,y_fit):
    """
    
    """
    pass

def _calc_R2(y_obs,y_fit):
    """
    
    """
    pass


def makeManualCalibration(paramset,rate_dict,plot=False):
    """
    
    """

    for pop_label in rate_dict.keys():
        logging.info("Updating parameter values for population=%s"%(pop_label))
        _setRateCalibration(paramset, rate_dict[pop_label], pop_label)
    

def _setRateCalibration(parset,value_dictionary,pop_label):
    """
    Creats parset and applies a transition/rate change that is applicable for the entire duration of the simulation 
    i.e. individual timepoints are not able to be specified.
    
    Params:
        parset
        valueDictionary
        pop_label
        
        
    Examples
        rate_dict_adults = {'mortality_rate': [0.02],
                            'num_births'    : [0] }
        rate_dict_sac = {'mortality_rate': [0.2],
                         'num_births'    : [5000] }
        runRateCalibration(project,'myNewParset',rate_dict,'Adults')
        
    
    """
    
    
    # update the values for a given popuation
    for k,v in value_dictionary.iteritems():
        par_index = parset.par_ids['cascade'][k]
        if isinstance(v,list):
            parset.pars['cascade'][par_index].y[pop_label] = np.array(v)
        elif isinstance(v, float) or isinstance(v, int):
            parset.pars['cascade'][par_index].y[pop_label] = np.array([v])
        if len(parset.pars['cascade'][par_index].t[pop_label]) > 1: # i.e. there are multiple years that we've just ignored
            # then take the first
            parset.pars['cascade'][par_index].t[pop_label] = np.array([parset.pars['cascade'][par_index].t[pop_label][0]])





class Calibration():
    """
    
    
    """
    
    def performAutofit(project,name=None,maxiters=1000,maxtime=None):
        """
        Run an autofit and save resulting parameterset
        
        Params:
            project
            name        name of resulting parameterset
            maxiters    max number of maximum iterations
            maxtime
        
        Currently, autofit is separated from loadAutofit to allow later interfacing to 
        other FEs.  
        
        """
        self.maxiters = maxiters
        
        loadAutofit(project)
        autofit()
    
    def loadAutofit(project):
        """
        
        TODO impement loading of autofit
        """
        pass
    
    
    def autofit():
        """
        
        TODO implement autofit
        """
        pass 
        