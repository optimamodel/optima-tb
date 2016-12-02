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

def calculateFitFunc(sim_data,sim_tvec,obs_data,metric):
    """
    
    """
    score = []
    char_labels = sim_data.keys()
    pop_labels = sim_data[0].keys()
    
    for char in char_labels:
        if char not in obs_data.keys():
            logger.info("Results: could not extract characteristic datapoint values for characteristic '%s' as this does not appear in databook"%char)
            continue
        logger.debug("calculating fit for char=%s"%char)
        for pop in pop_labels:
            y_obs = obs_data[char][pop]['y']
            t_obs = obs_data[char][pop]['t']
            # only grab sim_data and sim_tvec where there are corresponding values of obs_data
            t_indices = np.nonzero(np.in1d(sim_tvec, t_obs))[0]
            y_fit = sim_data[char][pop][t_indices]
            # calc and add to scores
            s = _calculateFitscore(y_obs, y_fit, metric)
            logger.debug("--- calc fit score = %g"%s)
            score.append(s)
    
    return score
    
def _calculateFitscore(y_obs, y_fit,metric="meansquare"):
    """
    

    """
    availfns = globals().copy()
    availfns.update(locals())
    try:
        return availfns.get('_calc_%s'%metric)(y_obs,y_fit)
    except:
        raise NotImplementedError("No method associated with _calc_%s (calibration.py)"%metric)


def _calc_meansquare(y_obs,y_fit):
    """
    Calcs the RMS error. 
    
    Note: could also use implementation from sklearn in future ... 
    """
    return np.sqrt(((y_fit - y_obs) ** 2).mean())


def _calc_R2(y_obs,y_fit):
    """
    
    """
    raise NotImplementedError


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
        