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

def update_parset(parset,y_factors,pars_to_adjust):
    # Insert updated y-values into the parset
    # - parset : a ParameterSet object
    # - y_factors : Array with as many elements as pars_to_adjust
    # - pars_to_adjust : Array of tuples (par_name,pop_label,...) with special value pop='all' supported
    #                    Must have as many elements as y_factors. pop=None is not allowed - it must be converted
    #                    to a full list of pops previously (in performAutofit)

    for i,x in enumerate(pars_to_adjust):
        par_name = x[0]
        pop_label = x[1]

        if par_name in parset.par_ids['cascade'] or par_name in parset.par_ids['characs']:
            par = parset.getPar(par_name)
            if pop_label == 'all':
                for pop in par.pops:
                    par.y_factor[pop] = y_factors[i]
            else:
                par.y_factor[pop_label] = y_factors[i]
        else:
            tokens = par_name.split('_from_')
            par = parset.transfers[tokens[0]][tokens[1]]
            par.y_factor[pop_label] = y_factors[i]


def calculateObjective(y_factors,pars_to_adjust,output_quantities,parset,project):
    # y-factors, array of y-factors to apply to specified output_quantities
    # pars_to_adjust - list of tuples (par_label,pop_label) recognized by parset.update()
    # output_quantities - a tuple like (pop,var,weight,metric) understood by model.getPop[pop].getVar(var)

    update_parset(parset,y_factors, pars_to_adjust)

    result = project.runSim(parset = parset, store_results = False)

    objective = 0.0

    for var_label,pop_label,weight,metric in output_quantities:
        var = result.model.getPop(pop_label).getVariable(var_label)
        if var_label in project.data['characs'].keys():
            target = project.data['characs'][var_label][pop_label]
        elif var_label in project.data['linkpars'].keys():
            target = project.data['linkpars'][var_label][pop_label]
        else:
            raise NotImplementedError

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

def _calc_fractional(y_obs,y_fit):
    return np.abs((y_fit-y_obs)/np.clip(y_obs,1,None)) # Use clipping to handle cases where data has value 0

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

    # Expand out pop=None in pars_to_adjust
    p2 = []
    for par_tuple in pars_to_adjust:
        if par_tuple[1] is None:  # If the pop name is None
            par = parset.getPar(par_tuple[0])
            for pop_label in par.pops:
                p2.append((par_tuple[0], pop_label, par_tuple[2], par_tuple[3]))
        else:
            p2.append(par_tuple)
    pars_to_adjust = p2

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
        par_label, pop_label, scale_min, scale_max = x
        if par_label in parset.par_ids['cascade'] or par_label in parset.par_ids['characs']:
            par = parset.getPar(par_label)
            if pop_label == 'all':
                x0.append(np.mean([par.y_factor[p] for p in par.pops]))
            else:
                x0.append(par.y_factor[pop_label])
        else:
            tokens = par_label.split('_from_')
            par = parset.transfers[tokens[0]][tokens[1]]
            x0.append(par.y_factor[pop_label])
        xmin.append(scale_min)
        xmax.append(scale_max)

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

    for i,x in enumerate(pars_to_adjust):
        par_label = x[0]
        pop_label = x[1]

        if par_label in parset.par_ids['cascade'] or par_label in parset.par_ids['characs']:
            par = args['parset'].getPar(par_label)

            if pop_label is None or pop_label == 'all':
                for pop in par.pops:
                    print '%s - %s, scale=%.2f' % (par_label, pop, par.y_factor[pop])
            else:
                print '%s - %s, scale=%.2f' % (par_label, pop_label, par.y_factor[pop_label])

        else: # For now, must be in there...
            tokens = par_label.split('_from_')
            par = args['parset'].transfers[tokens[0]][tokens[1]]
            print '%s - %s, scale=%.2f' % (tokens[0], tokens[1], par.y_factor[pop_label])

    for i,x in enumerate(pars_to_adjust):
        par_label = x[0]
        pop_label = x[1]

        if par_label in parset.par_ids['cascade'] or par_label in parset.par_ids['characs']:
            par = args['parset'].getPar(par_label)

            if pop_label is None or pop_label == 'all':
                for pop in par.pops:
                    # parset.getPar('b_rate').y_factor['0-2']=1.06
                    # parset.transfers['aging']['0-2'].y_factor['3-14'] = 1.01
                    #
                    print "parset.getPar('%s').y_factor['%s']=%.2f"  % (par_label, pop, par.y_factor[pop])
            else:
                print "parset.getPar('%s').y_factor['%s']=%.2f" % (par_label, pop_label, par.y_factor[pop_label])

        else: # For now, must be in there...
            tokens = par_label.split('_from_')
            par = args['parset'].transfers[tokens[0]][tokens[1]]
            print "parset.transfers['%s']['%s'].y_factor['%s']=%.2f" % (tokens[0], tokens[1], pop_label,par.y_factor[pop_label])

    return args['parset']



def calibrate_demographics(project,parset,max_time=60):

    # Adjust birth rates
    birth_pars = []
    for pop in parset.getPar('b_rate').pops:
        birth_pars.append(('b_rate',pop,0.5,4.0))
    birth_pars = [birth_pars[0]] # Don't touch migration actually

    # Adjust all transfer parameters
    transfer_pars = []
    for x in parset.transfers.values(): # for transfer type
        for y in x.values(): # for from_pop
            for pop in y.pops:
                transfer_pars.append((y.label,pop,0.1,10.0))

    death_pars = []
    # for pop in parset.getPar('doth_rate').pops:
    #     birth_pars.append(('doth_rate',pop,0.6,1.4))
    death_pars.append(('doth_rate','15-64 (HIV+)',0.1,10.0))
    death_pars.append(('doth_rate','65+ (HIV+)',0.1,10.0))

    emi_pars = []
    emi_pars.append(('emi_rate','all',0.1,2.0))

    pars_to_adjust = birth_pars + transfer_pars + death_pars + emi_pars

    # Collate the output demographic quantities (just alive for all pops)
    output_quantities = []
    for pop in parset.pop_labels:
        output_quantities.append(('alive',pop,1.0,"fractional"))

    calibrated_parset = performAutofit(project, parset, pars_to_adjust, output_quantities,max_time=max_time)

    return calibrated_parset


def calibrate_foi(project,parset,max_time=60):

    # Adjust force of infection
    foi_pars = []
    for pop in parset.getPar('spd_infxness').pops:
        foi_pars.append(('spd_infxness',pop,0.01,100))

    pars_to_adjust = foi_pars

    # Collate the output demographic quantities (just alive for all pops)
    output_quantities = []
    for pop in parset.pop_labels:
        output_quantities.append((pop,'l_inf',1.0,"fractional"))
        output_quantities.append((pop,'ac_inf',1.0,"fractional"))

    calibrated_parset = performAutofit(project, parset, pars_to_adjust, output_quantities,max_time=max_time)

    return calibrated_parset



