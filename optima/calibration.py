import logging
logger = logging.getLogger(__name__)

from utils import OptimaException, tic, toc, odict
import settings
import asd
from parameters import ParameterSet

import numpy as np
from copy import deepcopy as dcp

"""

Calibration and sensitivity analysis

"""


def performSensitivityAnalysis(project,transmission,epsilon=0.05):
    """
    
    TODO: implement performSensitivityAnalysis
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
            logger.debug("Results: could not extract characteristic datapoint values for characteristic '%s' as this does not appear in databook"%char)
            continue
        #logger.debug("calculating fit for char=%s"%char)
        for pop in pop_labels:
            y_obs = obs_data[char][pop]['y']
            t_obs = obs_data[char][pop]['t']
            # only grab sim_data and sim_tvec where there are corresponding values of obs_data
            t_indices = np.nonzero(np.in1d(sim_tvec, t_obs))[0]
            y_fit = sim_data[char][pop][t_indices]
            # calc and add to scores
            s = _calculateFitscore(y_obs, y_fit, metric)
            #logger.debug("---- for values obs = "+' '.join("%.2f"%ii for ii in y_obs)+" and yfit = "+' '.join("%.2f"%ii for ii in y_fit))
            #logger.debug("-------- calc fit score = %s"%' '.join("%.2f"%ii for ii in s))
            score.append(s)
    
    return np.concatenate(score).ravel()
    
def _getFitscoreFunc(metric):
    """
    
    
    """
    availfns = globals().copy()
    availfns.update(locals())
    try:
        return availfns.get('_calc_%s'%metric)
    except:
        raise NotImplementedError("No method associated with _calc_%s (calibration.py)"%metric)
        
    
def _calculateFitscore(y_obs, y_fit,metric="meansquare"):
    """
    

    """
    return _getFitscoreFunc(metric)(y_obs,y_fit)
    

def _calc_meansquare(y_obs,y_fit):
    """
    Calcs the RMS error. 
    
    Note: could also use implementation from sklearn in future ... 
    """
    return np.sqrt(((y_fit - y_obs) ** 2).mean())


def _calc_wape(y_obs,y_fit):
    """
    Calculates the weighted absolute percentage error 
    """
    return abs(y_fit - y_obs) / (y_obs.mean() + settings.TOLERANCE)


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
    
    # update the values for a given population
    for k,v in value_dictionary.iteritems():
        par_index = parset.par_ids['cascade'][k]
        if isinstance(v,list):
            parset.pars['cascade'][par_index].y[pop_label] = np.array(v)
        elif isinstance(v, float) or isinstance(v, int):
            parset.pars['cascade'][par_index].y[pop_label] = np.array([v])
        if len(parset.pars['cascade'][par_index].t[pop_label]) > 1: # i.e. there are multiple years that we've just ignored
            # then take the first
            parset.pars['cascade'][par_index].t[pop_label] = np.array([parset.pars['cascade'][par_index].t[pop_label][0]])


def performAutofit(project,paramset,new_parset_name,target_characs=None,**calibration_settings):
    """
    Run an autofit and save resulting parameterset
    
    Params:
        project
        paramset
        new_parset_name     name of resulting parameterset
        calibration_settings
   
    """
    # setup:
    useYFactor = calibration_settings['useYFactor']
    logger.debug("Autofit: useYFactor == %s"%useYFactor)
    metric = project.settings.fit_metric
    paramvec,minmax,casc_labels = paramset.extract(getMinMax=True,getYFactor=useYFactor)  # array representation of initial values for p0, with bounds
    if len(paramvec) == 0:
        raise OptimaException("No available cascade parameters to calibrate during autofitting. Please set at least one 'Calibrate?' value to be not equal to %g"%settings.DO_NOT_SCALE)
    
    mins, maxs = zip(*minmax)
    sample_param = dcp(paramset)   # ParameterSet created just to be overwritten
    sample_param.name = "calibrating"
    
    if target_characs is None: 
        # if no targets characteristics are supplied, then we autofit to all characteristics 
        target_data_characs = project.data['characs']
        logger.info("Autofit: fitting to all characteristics")
    else:
        target_data_characs = odict()
        for k in target_characs:
            target_data_characs[k] = project.data['characs'][k]
        logger.info("Autofit: fitting to the following target characteristics =[%s]"%(",".join(target_characs)))
    
    print paramvec
    print minmax
    
    def objective_calc(p_est):
        ''' Function used by ASD algorithm to run and evaluate fit of parameter set'''    
        print p_est
        sample_param.update(p_est,isYFactor=useYFactor)
        results = project.runSim(parameterset = sample_param)
        datapoints = results.getCharacteristicDatapoints()
        score = calculateFitFunc(datapoints,results.t_observed_data,target_data_characs,metric)
        return score
    
    
    parvecnew, fval, exitflag, output = asd.asd(objective_calc, paramvec, xmin=mins,xmax=maxs,xnames=casc_labels,**calibration_settings)
    
    print "pre ---------------------------"
    print sample_param.pars['cascade']
    print "/pre ---------------------------"
    sample_param.update(parvecnew,isYFactor=useYFactor)
    sample_param._updateFromYFactor()
    sample_param.name = new_parset_name
    print "post ---------------------------"
    print sample_param.pars['cascade']
    print "/post ---------------------------"
  
    return sample_param
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        