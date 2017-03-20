from optima_tb.asd import asd
from optima_tb.model import runModel
from optima_tb.utils import odict, OptimaException, tic, toc

import logging
import logging.config
logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger()

import numpy as np
from copy import deepcopy as dcp

def constrainAllocation(alloc, options):

    alloc = np.array(alloc)     # Converting to np array just in case.

    # Convert negative allocation values to zeros.
    alloc[alloc < 0.0] = 0.0    
    
    if 'total' in options['constraints']:
        sum_current = sum(alloc)
        if sum_current == 0: raise OptimaException('ERROR: Allocation was constrained to have a sum of zero during optimization.')
        alloc *= options['constraints']['total']/sum_current
        
    
    return dcp(alloc)

def calculateObjective(alloc, settings, parset, progset, options, algorithm_refs = None):
    '''
    Calculates the objective function value for a certain allocation.
    The algorithm_refs dictionary is for storing algorithm-useful values during the process, e.g. the previous allocation tested.
    '''
    
    print 'Unconstrained alloc...'
    print alloc
    print sum(alloc)
    
    alloc = constrainAllocation(alloc = alloc, options = options)   
    
    print 'Constrained alloc...'
    print alloc
    print sum(alloc)
    
    # Makes sure that time is not wasted on an allocation that was run in the previous step.
    # Avoids useless iteration tests where a program is given negative funding, then constrained to be zero once more.
    # This paragraph is for checking identical allocations.    
    identical_alloc = False
    if 'previous_alloc' in algorithm_refs:
        if not algorithm_refs['previous_alloc'] is None:
            identical_alloc = True
            for k in xrange(len(alloc)):
                if not alloc[k] == algorithm_refs['previous_alloc'][k]:
                    identical_alloc = False
                    break
        algorithm_refs['previous_alloc'] = dcp(alloc)   # Store the new allocation for the next iteration.
    
    # This if-else clause prevents the actual expensive model run.
    if identical_alloc and 'previous_results' in algorithm_refs and not algorithm_refs['previous_results'] is None:
        results = algorithm_refs['previous_results']
        logger.warn("Optimization iteration is testing the same allocation. Skipping model run.")
    else:
        options_iter = dcp(options)
        for k in xrange(len(alloc)):
            options_iter['init_alloc'][k] = alloc[k]
        t = tic()
        results = runModel(settings = settings, parset = parset, progset = progset, options = options_iter)
        print toc(t)

    if 'previous_results' in algorithm_refs:
        algorithm_refs['previous_results'] = dcp(results)   # Store the new results for the next iteration.
    
    index_start = np.where(results.sim_settings['tvec'] >= options['progs_start'])[0][0]
#    print index_start
#    print results.sim_settings['tvec'][index_start:]
    objective = 0.0
    
    for objective_label in options['objectives']:
#        charac_label = options['objective_weight'].keys()[0]
        weight = 1.0
        if 'weight' in options['objectives'][objective_label]: weight = options['objectives'][objective_label]['weight']
        
        # Handle an objective that wants to know the value in a particular year, e.g. a snapshot of a characteristic.
        if 'year' in options['objectives'][objective_label]:
            index_start = np.where(results.sim_settings['tvec'] >= options['objectives'][objective_label]['year'])[0][0]
            for pop_label in results.outputs[objective_label].keys():
                objective += results.outputs[objective_label][pop_label][index_start] * weight
            index_start = np.where(results.sim_settings['tvec'] >= options['progs_start'])[0][0]
        # Handle the default objective that wants to know a cumulative value. Timestep scaling required.
        else:
            for pop_label in results.outputs[objective_label].keys():
                objective += sum(results.outputs[objective_label][pop_label][index_start:]) * results.dt * weight
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
        
    for prog in progset.progs:
        if prog.label not in options['init_alloc']:
            options['init_alloc'][prog.label] = prog.getDefaultBudget()
    
    # Convert alloc into an ordered dictionary so that keys and values are ordered during optimisation.
    options['init_alloc'] = odict(options['init_alloc'])
    
    # If user has not supplied constraints, then make the optimisation a redistribution of default budgets.
    if 'constraints' not in options:
        options['constraints'] = {}
    if 'total' not in options:
        options['constraints']['total'] = sum(options['init_alloc'].values())
        
    if 'objectives' not in options:
        options['objectives'] = {settings.charac_pop_count : {'weight':-1,'year':2030.0}}
            
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
            'options':options,
            'algorithm_refs':{'previous_alloc':None, 'previous_results':None}}
    
    alloc = dcp(np.array(options['init_alloc'].values()))
    alloc_new, obj_vals, exit_reason = asd(calculateObjective, alloc, args=args, maxiters=7)#, xmin=xmin, maxtime=maxtime, maxiters=maxiters, verbose=verbose, randseed=randseed, label=thislabel, **kwargs)
    
#    print alloc_new
#    print obj_vals
#    print exit_reason
    
    results = (alloc_new, obj_vals, exit_reason)
    
    return results
    