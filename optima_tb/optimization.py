from optima_tb.asd import asd
from optima_tb.model import runModel
from optima_tb.utils import odict

import logging
import logging.config
logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

import numpy as np
from copy import deepcopy as dcp


def calculateObjective(alloc, settings, parset, progset, options):
    
    options_iter = dcp(options)
    for k in xrange(len(alloc)):
        options_iter['init_alloc'][k] = alloc[k]
    
    results = runModel(settings = settings, parset = parset, progset = progset, options = options_iter)
    
    index_start = np.where(results.sim_settings['tvec'] >= options['progs_start'])[0][0]
#    print index_start
#    print results.sim_settings['tvec'][index_start:]
    charac_label = options['objective_weight'].keys()[0]
    objective = 0
    for pop_label in results.outputs[charac_label].keys():
        objective += sum(results.outputs[charac_label][pop_label][index_start:])*results.dt*options['objective_weight'][charac_label]
#    print objective
    
    return objective


def optimizeFunc(settings, parset, progset, options = None):
    
    if options is None:
        logger.warn("An options dictionary was not supplied for optimisation. A default one will be constructed.")
        options = dict()
        
    if 'progs_start' not in options:
        options['progs_start'] = 2015.0
        logger.warn("Programs will start overwriting calibrated parameter values in the year: %f" % options['progs_start'])

    if 'init_alloc' not in options:
        options['init_alloc'] = {}
        
    total_budget = 0
    for prog in progset.progs:
        if prog.label not in options['init_alloc']:
            options['init_alloc'][prog.label] = prog.getDefaultBudget()
        total_budget += options['init_alloc'][prog.label]
    
    # Convert alloc into an ordered dictionary.
    options['init_alloc'] = odict(options['init_alloc'])
        
    if 'objective_weight' not in options:
        options['objective_weight'] = {settings.charac_pop_count : -1.0}
            
#    print options
#    print total_budget

#    results = runModel(settings = settings, parset = parset, progset = progset, options = options)
#    
#    index_start = np.where(results.sim_settings['tvec'] >= options['progs_start'])[0][0]
#    print index_start
#    print results.sim_settings['tvec'][index_start:]
#    charac_label = options['objective_weight'].keys()[0]
#    objective = 0
#    for pop_label in results.outputs[charac_label].keys():
#        objective += sum(results.outputs[charac_label][pop_label][index_start:])*results.dt*options['outcome_weight'][charac_label]
#    print objective
    
    args = {'settings':settings,
            'parset':parset, 
            'progset':progset,
            'options':options}
    
    alloc = dcp(options['init_alloc'].values())
    alloc_new, obj_vals, exit_reason = asd(calculateObjective, alloc, args=args, maxiters=3)#, xmin=xmin, maxtime=maxtime, maxiters=maxiters, verbose=verbose, randseed=randseed, label=thislabel, **kwargs)
    
#    print alloc_new
#    print obj_vals
#    print exit_reason
    
    results = (alloc_new, obj_vals, exit_reason)
    
    return results
    