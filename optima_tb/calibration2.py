import logging
logger = logging.getLogger(__name__)

from optima_tb.utils import OptimaException, tic, toc, odict
import optima_tb.settings as settings
from optima_tb.asd import asd
import optima_tb.plotting2 as oplt

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

def update_parset(parset,y_factors,pars_to_adjust):
    for i,x in enumerate(pars_to_adjust):
        par_label = x[0]
        pop_label = x[1]

        if par_label in parset.par_ids['cascade'] or par_label in parset.par_ids['characs']:
            parset.update([y_factors[i]], [(par_label,pop_label)], isYFactor=True)
        else: # For now, must be in there...
            tokens = par_label.split('_from_')
            par = parset.transfers[tokens[0]][tokens[1]]
            par.y_factor[pop_label] = y_factors[i]

def calculateObjective(y_factors,pars_to_adjust,output_quantities,parset,project):
    # y-factors, array of y-factors to apply to specified output_quantities
    # pars_to_adjust - list of tuples (par_label,pop_label) recognized by parset.update()
    # output_quantities - a tuple like (pop,var,weight,metric) understood by model.getPop[pop].getVar(var)
    # TODO - support Link flow rates and parameters, need to auto adjust

    update_parset(parset,y_factors, pars_to_adjust)

    result = project.runSim(parset = parset, store_results = False)

    objective = 0.0

    for pop_label,var_label,weight,metric in output_quantities:
        var = result.model.getPop(pop_label).getVariable(var_label)
        if isinstance(var[0],Characteristic):
            target = project.data['characs'][var_label][pop_label]
        else:
            raise OptimaException('Not yet implemented')

        y = target['y']
        y2 = interpolateFunc(var[0].t,var[0].vals,target['t'])

        objective += weight* sum(_calculateFitscore(y, y2, metric))

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

def performAutofit(proj,parset,pars_to_adjust,output_quantities,max_time=60):
    """
    Run an autofit and save resulting parameterset

    pars_to_adjust - list of tuples, (par,pop,scale) where allowed Y-factor range is 1-scale to 1+scale
    Params:
        project
        paramset
        new_parset_name     name of resulting parameterset
        target_characs      a list of characteristic and population label pairs
        calibration_settings
        useYFactor            boolean of flag whether we should use yvalues directly, or y_factor
   
    """

    args = {
        'project': proj,
        'parset': dcp(parset),
        'pars_to_adjust': pars_to_adjust,
        'output_quantities': output_quantities,
    }

    x0 = []
    xmin = []
    xmax = []
    for i,x in enumerate(pars_to_adjust):
        par_label, pop_label, scale = x
        if par_label in parset.par_ids['cascade'] or par_label in parset.par_ids['characs']:
            par = parset.getPar(par_label)
            x0.append(par.y_factor[pop_label])
        else:
            tokens = par_label.split('_from_')
            par = parset.transfers[tokens[0]][tokens[1]]
            x0.append(par.y_factor[pop_label])
        xmin.append(1-scale)
        xmax.append(1+scale)

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

    x1, _, _ = asd(calculateObjective, x0, args, **optim_args)

    update_parset(args['parset'],x1,pars_to_adjust)

    return args['parset']



def calibrate_demographics(project,parset,max_time=60):

    # Adjust birth rate
    birth_pars = [('b_rate',parset.getPar('b_rate').pops[0],0.9)]

    # Adjust all transfer parameters
    transfer_pars = []
    for x in parset.transfers.values(): # for transfer type
        for y in x.values(): # for from_pop
            for pop in y.pops:
                transfer_pars.append((y.label,pop,0.9))

    pars_to_adjust = birth_pars + transfer_pars

    # Collate the output demographic quantities (just alive for all pops)
    output_quantities = []
    for pop in parset.pop_labels:
        output_quantities.append((pop,'alive',1.0,"wape"))

    calibrated_parset = performAutofit(project, parset, pars_to_adjust, output_quantities,max_time=max_time)

    return calibrated_parset

# Adjust death rates


