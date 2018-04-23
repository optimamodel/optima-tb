import logging
logger = logging.getLogger(__name__)

from optima_tb.utils import OptimaException, tic, toc, odict
import optima_tb.settings as settings
import optima_tb.asd as asd
from optima_tb.parameters import ParameterSet

from optima_tb.model import Compartment, Characteristic, Link, Parameter

import numpy as np
from copy import deepcopy as dcp

from optima_tb.interpolation import interpolateFunc

"""

Calibration and sensitivity analysis

"""

# CALIBRATION LAYOUT
# Auto calibration proceeds using y-factors. The input should thus be a list
# of parameters to calibrate. The output metric would be a list of characteristics
# or flow rates that need to be matched. So calibration requires
# - A list of parameters to adjust
# - A list of output characteristics to match
#
# For example, the input parameters could be a list of birth rates and transfer rates
# and HIV infection rates, and the output parameters could be a list of characteristics
# In theory these could be population-specific?

def calculateObjective(y_factors,pars_to_adjust,output_quantities,parset,project):
    # y-factors, array of y-factors to apply to specified output_quantities
    # pars_to_adjust - list of tuples (par_label,pop_label) recognized by parset.update()
    # output_quantnties - a tuple like (pop,var,weight,metric) understood by model.getPop[pop].getVar(var)
    # TODO - support Link flow rates and parameters, need to auto adjust

    parset.update(y_factors, pars_to_adjust, isYFactor=True)

    results = project.runSim(parset = parset, store_results = False)

    objective = 0.0

    for pop_label,var_label,weight,metric in output_quantities
        var = result.model.getPop(pop_label).getVariable(var_label)
        if isinstance(var,Characteristic):
            target = project.data['characs'][var_label][pop_label]
        else:
            raise OptimaException('Not yet implemented')

        y = target['y']
        y2 = interpolateFunc(var[0].t,var[0].vals,target['t'])

        objective += weight* _calculateFitscore(y, y2, metric)

    return objective

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


def performAutofit(project,parset,pars_to_adjust,output_quantities):
    """
    Run an autofit and save resulting parameterset
    
    Params:
        project
        paramset
        new_parset_name     name of resulting parameterset
        target_characs      a list of characteristic and population label pairs
        calibration_settings
        useYFactor            boolean of flag whether we should use yvalues directly, or y_factor
   
    """

    args = {
        'project': project,
        'parset': dcp(parset),
        'pars_to_adjust': pars_to_adjust,
        'output_quantities': output_quantities,
    }

    x0 = parset.getPar(adjust[0],adjust[1]).
    optim_args = {
                 'stepsize': proj.settings.autofit_params['stepsize'],
                 'maxiters': proj.settings.autofit_params['maxiters'],
                 'maxtime': proj.settings.autofit_params['maxtime'],
                 'sinc': proj.settings.autofit_params['sinc'],
                 'sdec': proj.settings.autofit_params['sdec'],
                 'fulloutput': False,
                 'reltol': None,
                 'xmin': xmin,
                 'xmax': xmax,
                 }

    if not max_time is None:
        optim_args['maxtime'] = max_time

    best_attribute_list, _, _ = asd(_objective, x0, args, **optim_args)

    # setup:
    
    # setup for cascade parameters
    paramvec,minmax,par_pop_labels = paramset.extract(settings=project.settings,getMinMax=True,getYFactor=useYFactor)  # array representation of initial values for p0, with bounds
    # setup for characteristics
    compartment_init,charac_pop_labels = paramset.extractEntryPoints(project.settings,useInitCompartments=useInitCompartments)
    # min maxes for compartments are always (0,np.inf):
    charac_minmax = [(0,np.inf) for i in charac_pop_labels]
    minmax += charac_minmax
    
#    print par_pop_labels
#    print charac_pop_labels
    
    
    if len(paramvec)+len(compartment_init) == 0:
        raise OptimaException("No available cascade parameters or initial characteristic sizes to calibrate during autofitting. At least one 'Autocalibrate' cascade value must be listed something other than 'n' or '-1'.")
    
    mins, maxs = zip(*minmax)
    sample_param = dcp(paramset)   # ParameterSet created just to be overwritten
    sample_param.name = "calibrating"
    
#    print mins
#    print maxs
    
    if target_characs is None: 
        # if no targets characteristics are supplied, then we autofit to all characteristics 
        target_data_characs = dcp(project.data['characs'])
        for label in target_data_characs.keys():
            target_data_characs[label]['pops_to_fit'] = {pop_label:True for pop_label in paramset.pop_labels}
        logger.info("Autofit: fitting to all characteristics")
    else:
        target_data_characs = odict()
        for pair in target_characs:
            if not pair[0] in target_data_characs:
                target_data_characs[pair[0]] = dcp(project.data['characs'][pair[0]])
            if 'pops_to_fit' not in target_data_characs[pair[0]]:
                target_data_characs[pair[0]]['pops_to_fit'] = {}
            target_data_characs[pair[0]]['pops_to_fit'][pair[1]] = True
        logger.info("Autofit: fi
    parvecnew, fval, details = asd.asd(calculateObjective, paramvec+compartment_init, args={}, xmin=mins, xmax=maxs, **calibration_settings)
    
#    # Compare old and new values 
#    print paramvec
#    print parvecnew[:len(paramvec)]
#    print compartment_init
#    print parvecnew[len(paramvec):]
    
    sample_param.update(parvecnew, par_pop_labels+charac_pop_labels, isYFactor=useYFactor)
#    sample_param._updateFromYFactor()
    sample_param.name = new_parset_name
    
    return sample_param
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        